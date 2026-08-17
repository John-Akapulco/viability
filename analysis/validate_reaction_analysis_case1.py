"""First real-data run of reaction_analysis (as opposed to the synthetic
Na/K/NaK/Na2K fixtures tests/test_delta.py already covers) -- validates the
new schema-driven pipeline against the existing ad hoc reaction_icohp.py
(mission #5, analysis/reaction_icohp_case1.csv) on 3 already-computed real
compounds spanning different bond_type families, before deciding whether to
populate reactions_dataset/ at full scale.

Cross-check logic: reaction_icohp.py's delta_icohp_per_atom uses the
OPPOSITE sign convention (reactant - products) from reaction_analysis's
delta_per_atom_eV (products - reactants, see delta.py docstring), and is
computed on the reactant's own actual cell/atom-count. reaction_analysis's
delta_per_atom_eV is per-formula-unit-normalized (Z divided out) instead.
Both are intensive per-atom quantities of the same underlying reaction, so
after flipping the sign they should agree to floating-point precision --
this is not expected to be approximate, since both read the exact same
ICOHPLIST.lobster files and both sum unfiltered over every bond label.

Populates reactions_dataset/entries/ and reactions_dataset/reactions/ with
these 5 compounds + 2 reactions as the first real (non-synthetic) content
in that directory -- previously empty by design (see its README).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from reaction_analysis.balance import derive_element_coefficients
from reaction_analysis.delta import compute_delta
from reaction_analysis.parse_lobster import parse_compound_entry
from reaction_analysis.schema import Reaction, ReactionMember

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
DATASET_ROOT = REPO_ROOT / "reactions_dataset"

# (target compound_id, [element compound_ids], reaction_icohp_case1.csv compound_id)
CASES = [
    ("extension_CaO_mp-2605", ["extension_Ca_mp-21", "extension_O2_mp-1524462"], "extension_CaO_mp-2605"),
    ("extension_Ca3N2_mp-844", ["extension_Ca_mp-21", "gasref_N2_dimerbox"], "extension_Ca3N2_mp-844"),
    ("exp_metastable_AsPd2_mp-1080831", ["extension_As_mp-158", "extension_Pd_mp-2"], "exp_metastable_AsPd2_mp-1080831"),
]


def _write_entry(compound_id: str, role: str) -> None:
    entry_dir = DATASET_ROOT / "entries" / compound_id
    entry_dir.mkdir(parents=True, exist_ok=True)
    src = STRUCTURES_ROOT / compound_id
    entry = parse_compound_entry(src, role=role, compound_id=compound_id)
    (entry_dir / "entry.json").write_text(entry.model_dump_json(indent=2))
    for fname in ("ICOHPLIST.lobster", "ICOBILIST.lobster", "CONTCAR"):
        fsrc = src / fname
        if fsrc.exists():
            (entry_dir / fname).write_bytes(fsrc.read_bytes())


def main() -> None:
    (DATASET_ROOT / "entries").mkdir(parents=True, exist_ok=True)
    (DATASET_ROOT / "reactions").mkdir(parents=True, exist_ok=True)

    old_by_id = {}
    with open(REPO_ROOT / "analysis" / "reaction_icohp_case1.csv") as f:
        import csv
        for row in csv.DictReader(f):
            old_by_id[row["compound_id"]] = row

    entries = {}
    reactions = []
    for target_id, element_ids, old_id in CASES:
        all_ids = [target_id, *element_ids]
        for cid in all_ids:
            if cid not in entries:
                role = "element" if cid in element_ids and cid not in (t for t, _, _ in CASES) else "target"
                # an id used as both a target elsewhere and an element here never
                # happens in this fixed case list, so this role pick is unambiguous
                _write_entry(cid, role)
                from reaction_analysis.schema import CompoundEntry
                entries[cid] = CompoundEntry.model_validate_json(
                    (DATASET_ROOT / "entries" / cid / "entry.json").read_text()
                )

        target = entries[target_id]
        element_entries = {}
        for el_id in element_ids:
            comp = entries[el_id].composition_per_formula_unit()
            assert len(comp) == 1, f"{el_id} is not a pure element: {comp}"
            (element,) = comp.keys()
            element_entries[element] = entries[el_id]

        products = derive_element_coefficients(target, element_entries)
        reaction = Reaction(
            reaction_id=f"{target_id}__decomposition",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id=target_id, coefficient=1.0)],
            products=products,
        )
        (DATASET_ROOT / "reactions" / f"{reaction.reaction_id}.json").write_text(
            reaction.model_dump_json(indent=2)
        )
        reactions.append((reaction, old_id))

    print(f"{'compound':<32} {'new delta/atom (prod-react)':>28} {'old -delta/atom (react-prod flipped)':>38} {'match':>7}")
    all_ok = True
    for reaction, old_id in reactions:
        results = compute_delta(reaction, entries)
        result = results[0]
        old_row = old_by_id[old_id]
        old_per_atom_flipped = -float(old_row["delta_icohp_per_atom"])
        ok = abs(result.delta_per_atom_eV - old_per_atom_flipped) < 1e-4
        all_ok &= ok
        print(
            f"{old_id:<32} {result.delta_per_atom_eV:>28.6f} {old_per_atom_flipped:>38.6f} "
            f"{'OK' if ok else 'MISMATCH'}"
        )
        if result.error:
            print(f"  ERROR: {result.error}")
            all_ok = False
        if result.warnings:
            print(f"  warnings: {result.warnings}")

    print()
    print("ALL MATCH" if all_ok else "SOME MISMATCHES -- do not trust reaction_analysis until resolved")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
