"""Extend reaction_analysis's decomposition_to_elements case from the
6-compound hand-picked validation batch (validate_reaction_analysis_case1.py)
to the full dataset, at the same scope as reaction_icohp.py's case 1
(compute_reaction_icohp_case1.py, mission #5): every non-single-element
compound in mp_dataset/structures/ whose elements all have an entry in
ELEMENT_REFERENCE (reused directly from compute_reaction_icohp_case1.py,
not re-derived) and whose elemental references are themselves computed.

Populates reactions_dataset/entries/ (one entry.json per compound actually
used -- elements are shared/cached across every target that needs them,
never re-parsed) and reactions_dataset/reactions/ (one
decomposition_to_elements Reaction per target), then runs the whole batch
through reaction_analysis.delta.compute_delta and cross-checks every
result's delta_per_atom_eV against reaction_icohp_case1.csv's
delta_icohp_per_atom (sign-flipped, see validate_reaction_analysis_case1.py
for why these should agree to floating-point precision) -- a much stronger
validation than the 3-compound spot check, since it covers every compound
both pipelines can compute.

Raw ICOHPLIST.lobster/CONTCAR files are deliberately NOT copied into
reactions_dataset/entries/<id>/ at this scale (unlike the small validation
batch) -- they already live in mp_dataset/structures/<id>/, which
CompoundEntry.source_path points back to; duplicating ~260 compounds'
worth of these files would be pure redundant storage for no provenance
gain, since the originals are already tracked in this repo.

Does NOT re-run the mission #5 statistical/correlation analysis
(stats_analysis_reaction_icohp.py / REPORT_reaction_icohp.md) -- the
per-atom cross-check below confirms this pipeline reproduces the same
numbers, so re-deriving the same correlations would be redundant, not a
new result.
"""

from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).parent))

from pymatgen.core import Composition  # noqa: E402

import compute_reaction_icohp_case1 as ri_case1  # noqa: E402
from reaction_analysis.balance import derive_element_coefficients  # noqa: E402
from reaction_analysis.delta import compute_delta  # noqa: E402
from reaction_analysis.parse_lobster import parse_compound_entry  # noqa: E402
from reaction_analysis.schema import CompoundEntry, Reaction, ReactionMember  # noqa: E402

STRUCTURES_ROOT = ri_case1.STRUCTURES_ROOT
ELEMENT_REFERENCE = ri_case1.ELEMENT_REFERENCE
is_computed = ri_case1.is_computed

DATASET_ROOT = REPO_ROOT / "reactions_dataset"
OLD_CSV = REPO_ROOT / "analysis" / "reaction_icohp_case1.csv"
OUT_CSV = REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv"
OUT_JSON = REPO_ROOT / "analysis" / "reaction_analysis_case1_full.json"
BONDTYPE_CSV = REPO_ROOT / "analysis" / "icohp_icobi_bondtype.csv"

_entry_cache: dict[str, CompoundEntry] = {}


def _load_bondtype_and_ismetal_maps() -> tuple[dict[str, str], dict[str, bool]]:
    """compute_icohp_icobi_bondtype.py's is_metal-first, first-shell-ICOBI
    classification (commit 60fe81a) -- used here instead of the raw
    mp_metadata.json bond_type, which is classify()'s composition-only
    heuristic (fetch_candidates.py) and leaves ~2/3 of this dataset's
    compounds (mostly is_metal=True compounds containing an anion-like
    element, deliberately excluded from classify()'s "metallic" bucket)
    with bond_type=None. Also returns its is_metal column (build_is_metal_map()
    in compute_icobi_antibonding_all.py: main-campaign percolation CSV,
    mp_metadata.json fallback otherwise -- 0 NaN across the whole dataset,
    vs. 188/591 structures dirs where meta["is_metal"] itself is None,
    mostly main-campaign compounds whose metadata predates is_metal being
    added directly). Both fall back to the raw meta[...] value per-row (see
    main()) for any compound not in this CSV, e.g. one added after the
    classifier was last run."""
    import pandas as pd

    if not BONDTYPE_CSV.exists():
        return {}, {}
    df = pd.read_csv(BONDTYPE_CSV)
    return (
        dict(zip(df["compound_id"], df["icobi_label"])),
        dict(zip(df["compound_id"], df["is_metal"])),
    )


def _get_entry(compound_id: str, role: str) -> CompoundEntry:
    if compound_id in _entry_cache:
        return _entry_cache[compound_id]
    entry_dir = DATASET_ROOT / "entries" / compound_id
    entry_dir.mkdir(parents=True, exist_ok=True)
    entry = parse_compound_entry(STRUCTURES_ROOT / compound_id, role=role, compound_id=compound_id)
    (entry_dir / "entry.json").write_text(entry.model_dump_json(indent=2))
    _entry_cache[compound_id] = entry
    return entry


def main() -> None:
    (DATASET_ROOT / "entries").mkdir(parents=True, exist_ok=True)
    (DATASET_ROOT / "reactions").mkdir(parents=True, exist_ok=True)

    n_ok = n_skipped_single = n_skipped_element = n_not_ready = n_failed = 0
    rows: list[dict] = []
    detail: dict[str, dict] = {}
    bond_type_map, is_metal_map = _load_bondtype_and_ismetal_maps()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
            meta_path = compound_dir / "mp_metadata.json"
            if not meta_path.exists() or not is_computed(compound_dir):
                continue
            compound_id = compound_dir.name
            meta = json.loads(meta_path.read_text())
            formula = meta.get("formula")
            if not formula:
                continue
            try:
                comp = Composition(formula)
            except Exception:
                continue
            elements = [str(e) for e in comp.elements]

            if len(elements) == 1:
                n_skipped_single += 1
                continue

            missing_refs = [e for e in elements if e not in ELEMENT_REFERENCE]
            if missing_refs:
                detail[compound_id] = {"error": f"no elemental reference for {missing_refs}"}
                n_skipped_element += 1
                continue

            ref_ids = [ELEMENT_REFERENCE[e] for e in elements]
            not_ready = [rid for rid in ref_ids if not is_computed(STRUCTURES_ROOT / rid)]
            if not_ready:
                detail[compound_id] = {"error": f"reference(s) not yet computed: {not_ready}"}
                n_not_ready += 1
                continue

            try:
                target = _get_entry(compound_id, role="target")
                element_entries = {}
                for el, ref_id in zip(elements, ref_ids):
                    ref_entry = _get_entry(ref_id, role="element")
                    ref_comp = ref_entry.composition_per_formula_unit()
                    assert set(ref_comp) == {el}, f"{ref_id} not pure {el}: {ref_comp}"
                    element_entries[el] = ref_entry

                products = derive_element_coefficients(target, element_entries)
                reaction = Reaction(
                    reaction_id=f"{compound_id}__decomposition",
                    type="decomposition_to_elements",
                    reactants=[ReactionMember(compound_id=compound_id, coefficient=1.0)],
                    products=products,
                )
                (DATASET_ROOT / "reactions" / f"{reaction.reaction_id}.json").write_text(
                    reaction.model_dump_json(indent=2)
                )

                all_entries = {compound_id: target, **{ref_id: element_entries[el] for el, ref_id in zip(elements, ref_ids)}}
                result = compute_delta(reaction, all_entries)[0]
                detail[compound_id] = result.model_dump()
                if result.error:
                    n_failed += 1
                    print(f"FAILED {compound_id}: {result.error}")
                    continue

                n_ok += 1
                rows.append({
                    "compound_id": compound_id,
                    "mp_id": meta.get("mp_id"),
                    "formula": meta.get("formula"),
                    "family": meta.get("family"),
                    "bond_type": bond_type_map.get(compound_id, meta.get("bond_type")),
                    "is_metal": is_metal_map.get(compound_id, meta.get("is_metal")),
                    "theoretical": meta.get("theoretical"),
                    "energy_above_hull_eV_per_atom": meta.get("energy_above_hull_eV_per_atom"),
                    "delta_per_formula_unit_eV": result.delta_per_formula_unit_eV,
                    "delta_per_atom_eV": result.delta_per_atom_eV,
                    "delta_per_bond_eV": result.delta_per_bond_eV,
                })
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                detail[compound_id] = {"error": f"{type(exc).__name__}: {exc}"}
                n_failed += 1
                print(f"FAILED {compound_id}: {type(exc).__name__}: {exc}")

    OUT_JSON.write_text(json.dumps(detail, indent=2, default=str))
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"\nWrote {len(rows)} rows to {OUT_CSV} "
        f"(ok={n_ok}, skipped_single_element={n_skipped_single}, "
        f"skipped_no_reference={n_skipped_element}, reference_not_ready={n_not_ready}, "
        f"failed={n_failed}); {len(_entry_cache)} distinct CompoundEntry objects parsed"
    )

    # Cross-check against the existing reaction_icohp.py case-1 numbers
    # (opposite sign convention, see module docstring)
    old_by_id = {}
    with OLD_CSV.open() as f:
        for row in csv.DictReader(f):
            old_by_id[row["compound_id"]] = row

    n_compared = n_match = n_mismatch = n_no_old = 0
    mismatches = []
    for row in rows:
        old_row = old_by_id.get(row["compound_id"])
        if old_row is None:
            n_no_old += 1
            continue
        n_compared += 1
        old_per_atom_flipped = -float(old_row["delta_icohp_per_atom"])
        if abs(row["delta_per_atom_eV"] - old_per_atom_flipped) < 1e-4:
            n_match += 1
        else:
            n_mismatch += 1
            mismatches.append((row["compound_id"], row["delta_per_atom_eV"], old_per_atom_flipped))

    print(
        f"\nCross-check vs reaction_icohp_case1.csv: {n_compared} compared "
        f"({n_match} match, {n_mismatch} mismatch, {n_no_old} in new but not in old CSV)"
    )
    for cid, new_v, old_v in mismatches[:20]:
        print(f"  MISMATCH {cid}: new={new_v:.6f} old(flipped)={old_v:.6f}")


if __name__ == "__main__":
    main()
