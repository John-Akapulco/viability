"""Case 1 (decomposition into elements) of reaction_icohp.py, run across
every compound in the dataset (186-compound main campaign + extension
batches), not just the 8 extension reactions already computable when
METRIC_DEFINITION_reaction_icohp.md was validated.

Needs one elemental reference compound_dir per element that appears
anywhere in the dataset (62 distinct elements). 12 already existed;
50 were added in bulk by download_elements_reference.py, plus a
dedicated Zn fix (extension_Zn_mp-79 replaces the pre-existing
theo_metastable_Zn_mp-2647117, which is theoretical=True bcc Zn --
wrong structure type, chosen by the main campaign's stratification
criteria for an unrelated purpose, not a valid elemental reference; using
it here would have been a real error, caught before use).

Pure single-element compounds (the reference compounds themselves, plus
any single-element compound_id elsewhere in the dataset) are skipped --
decomposition of an element into itself is degenerate (see
METRIC_DEFINITION_reaction_icohp.md section 5 "IDENTITY TEST": pymatgen's
Reaction balancer collapses reactant==product into a null reaction,
delta is not well-defined). Those belong to case 2 (polymorph comparison)
instead, not case 1.

Skipped, not guessed, if any element's reference is missing ICOHPLIST.lobster/
CONTCAR yet (job not finished) -- per-compound try/except, batch does not
die on one bad/incomplete compound, same convention as
compute_antibonding_all.py and compute_antibonding_extension.py.

Writes analysis/reaction_icohp_case1.json (full detail per compound) and
analysis/reaction_icohp_case1.csv (delta_icohp_per_atom next to each
compound's own metadata, for the follow-up stats script).
"""

import json
import sys
import warnings
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
import reaction_icohp as ri  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT_JSON = Path(__file__).parent / "reaction_icohp_case1.json"
OUT_CSV = Path(__file__).parent / "reaction_icohp_case1.csv"
BONDTYPE_CSV = Path(__file__).parent / "icohp_icobi_bondtype.csv"


def _load_bondtype_and_ismetal_maps() -> tuple[dict, dict]:
    """analysis/compute_icohp_icobi_bondtype.py's bond_type (is_metal-
    first, first-shell-ICOBI classification, commit 60fe81a) and its
    is_metal column (build_is_metal_map() in compute_icobi_antibonding_all.py:
    main-campaign percolation CSV, mp_metadata.json fallback otherwise --
    already covers the whole dataset with 0 NaN in icohp_icobi_bondtype.csv).
    Both preferred over the raw mp_metadata.json fields directly: bond_type
    is classify()'s composition-only heuristic, None for most is_metal=True
    compounds containing an anion-like element; is_metal itself is None for
    188/591 structures dirs (mostly main-campaign compounds whose metadata
    predates is_metal being added directly -- see build_dataset.py's "present
    for the 6 pilots, None for campaign compounds" comment). See
    populate_reaction_analysis_case1_full.py's identically-named helper for
    the same fix."""
    if not BONDTYPE_CSV.exists():
        return {}, {}
    df = pd.read_csv(BONDTYPE_CSV)
    return (
        dict(zip(df["compound_id"], df["icobi_label"])),
        dict(zip(df["compound_id"], df["is_metal"])),
    )

# element symbol -> reference compound_dir name
# K has no low-energy non-theoretical bcc entry available in MP (the real
# standard state for this alkali metal) -- exp_metastable_K_mp-10157
# (fcc, EAH=6.9 meV/at, non-theoretical) is the best available pick, same
# "not individually verified, best available" caveat as several of the
# download_elements_reference.py entries. Left as-is rather than guessed.
# N was extension_N2_mp-1059834 until 2026-08-17: that MP entry is a
# polymeric N-N solid (1.296 A, two contacts/atom), not molecular N2 (gas
# N-N is 1.10 A). Replaced by gasref_N2_dimerbox, an isolated-dimer
# relaxation (1.113 A, ICOHP -22.998 eV, matching Reitz & Dronskowski's
# -23.161 eV to 0.7%) built specifically to fix this reference.
ELEMENT_REFERENCE = {
    "Ca": "extension_Ca_mp-21",
    "C": "extension_Cgraphite_mp-48",
    "K": "exp_metastable_K_mp-10157",
    "Mn": "extension_Mn_mp-35",
    "N": "gasref_N2_dimerbox",
    "Na": "exp_stable_Na_mp-10172",
    "O": "extension_O2_mp-1524462",
    "Pb": "extension_Pb_mp-20483",
    "S": "extension_S_mp-77",
    "Si": "hull_covalent_Si_mp-149",
    "Ti": "extension_Ti_mp-46",
    "Zn": "extension_Zn_mp-79",
    "Ag": "extension_Ag_mp-124", "Al": "extension_Al_mp-134", "As": "extension_As_mp-158",
    "Au": "extension_Au_mp-81", "B": "extension_B_mp-160", "Ba": "extension_Ba_mp-122",
    "Be": "extension_Be_mp-87", "Br": "extension_Br_mp-23154", "Cd": "extension_Cd_mp-94",
    "Cl": "extension_Cl_mp-22848", "Co": "extension_Co_mp-102", "Cr": "extension_Cr_mp-90",
    "Cs": "extension_Cs_mp-1", "Cu": "extension_Cu_mp-30", "F": "extension_F_mp-561203",
    "Fe": "extension_Fe_mp-13", "Ga": "extension_Ga_mp-142", "Ge": "extension_Ge_mp-32",
    "H": "extension_H_mp-730101", "Hf": "extension_Hf_mp-103", "Hg": "extension_Hg_mp-10861",
    "I": "extension_I_mp-23153", "In": "extension_In_mp-1055994", "Ir": "extension_Ir_mp-101",
    "Li": "extension_Li_mp-1018134", "Mg": "extension_Mg_mp-153", "Mo": "extension_Mo_mp-129",
    "Nb": "extension_Nb_mp-75", "Ni": "extension_Ni_mp-23", "Os": "extension_Os_mp-49",
    "P": "extension_P_mp-568348", "Pd": "extension_Pd_mp-2", "Pt": "extension_Pt_mp-126",
    "Rb": "extension_Rb_mp-70", "Re": "extension_Re_mp-8", "Rh": "extension_Rh_mp-74",
    "Ru": "extension_Ru_mp-33", "Sb": "extension_Sb_mp-104", "Sc": "extension_Sc_mp-67",
    "Se": "extension_Se_mp-14", "Sn": "extension_Sn_mp-623511", "Sr": "extension_Sr_mp-139",
    "Ta": "extension_Ta_mp-50", "Tc": "extension_Tc_mp-113", "Te": "extension_Te_mp-19",
    "Tl": "extension_Tl_mp-82", "V": "extension_V_mp-146", "W": "extension_W_mp-91",
    "Y": "extension_Y_mp-112", "Zr": "extension_Zr_mp-131",
}


def is_computed(compound_dir: Path) -> bool:
    return (compound_dir / "ICOHPLIST.lobster").exists() and (compound_dir / "CONTCAR").exists()


def main():
    results = {}
    rows = []
    n_ok = n_skipped_element = n_skipped_single = n_not_ready = n_failed = 0
    bond_type_map, is_metal_map = _load_bondtype_and_ismetal_maps()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
            meta_path = compound_dir / "mp_metadata.json"
            if not meta_path.exists() or not is_computed(compound_dir):
                continue
            compound_id = compound_dir.name
            meta = json.loads(meta_path.read_text())
            if meta.get("quality_excluded"):
                results[compound_id] = {"error": f"quality_excluded: {meta.get('quality_excluded_reason', 'see mp_metadata.json')}"}
                continue
            formula = meta.get("formula")
            if not formula:
                continue
            try:
                comp = Composition(formula)
            except Exception:
                continue
            elements = [str(e) for e in comp.elements]

            if len(elements) == 1:
                results[compound_id] = {"error": "single-element compound -- case 2 (polymorph), not case 1"}
                n_skipped_single += 1
                continue

            missing_refs = [e for e in elements if e not in ELEMENT_REFERENCE]
            if missing_refs:
                results[compound_id] = {"error": f"no elemental reference for {missing_refs}"}
                n_skipped_element += 1
                continue

            ref_dirs = [STRUCTURES_ROOT / ELEMENT_REFERENCE[e] for e in elements]
            not_ready = [d.name for d in ref_dirs if not is_computed(d)]
            if not_ready:
                results[compound_id] = {"error": f"reference(s) not yet computed: {not_ready}"}
                n_not_ready += 1
                continue

            try:
                r = ri.reaction_delta_icohp(compound_dir, ref_dirs)
                results[compound_id] = r
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
                    "reaction_string": r["reaction_string"],
                    "delta_icohp_total": r["delta_icohp_total"],
                    "delta_icohp_per_atom": r["delta_icohp_per_atom"],
                })
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                results[compound_id] = {"error": f"{type(exc).__name__}: {exc}"}
                n_failed += 1
                print(f"FAILED {compound_id}: {type(exc).__name__}: {exc}")

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))
    print(
        f"\nWrote {len(results)} entries to {OUT_JSON} "
        f"(ok={n_ok}, skipped_single_element={n_skipped_single}, "
        f"skipped_no_reference={n_skipped_element}, reference_not_ready={n_not_ready}, "
        f"failed={n_failed})"
    )

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
