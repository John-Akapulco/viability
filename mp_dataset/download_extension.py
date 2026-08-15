"""Download a small, hand-picked set of compounds the user asked to add to
the study, on top of the 186-compound campaign (select_campaign.py) and the
6-compound pilot (fetch_candidates.py). Unlike those two, this list was not
selected by any automated criterion -- it is explicitly requested, and
several entries break assumptions the rest of the pipeline relies on
(missing from Materials Project entirely, or physically magnetic where the
rest of the campaign was filtered to non-magnetic). Both exceptions are
handled explicitly here and in prepare_extension_vasp_lobster.py, not
silently.

Two sources:
  - Materials Project (10 compounds): batched `material_ids` lookup, same
    pattern as download_campaign.py.
  - Crystallography Open Database (2 compounds, S4N2 and S4N4): absent from
    MP entirely (checked across the full N-S chemsys, no entry). Verified
    real, literature-backed COD entries used instead (DOIs recorded in
    metadata) -- CIFs already fetched to the scratchpad and parsed
    successfully before this script was written.

Writes one directory per compound under mp_dataset/structures/ (same
layout/tooling as every other compound in this project -- percolation_path.py,
prepare_vasp_lobster.py, etc. all just iterate that directory), each with
POSCAR + mp_metadata.json, family="extension".
"""

import json
import os
import warnings
from pathlib import Path

from mp_api.client import MPRester
from pymatgen.io.cif import CifParser

import sys

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from fetch_candidates import classify  # noqa: E402  (reused unmodified)

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()
OUT_ROOT = Path(__file__).parent / "structures"
COD_CIF_DIR = Path(
    "/tmp/claude-1002/-home-gilles-viability/254d15db-9c00-4d72-ae7b-fb85a1389cff/scratchpad/cod_check"
)

# (label, mp_id, requested_spacegroup, note)
MP_SOURCED = [
    ("CaO", "mp-2605", "Fm-3m", "rock-salt, on-hull"),
    ("CaN", "mp-1058549", "Fm-3m", "rock-salt-type CaN; 0.398 eV/at above hull"),
    ("Ca3N2", "mp-844", "Ia-3", "anti-bixbyite, on-hull, 40 sites -- larger than the usual <=20-site cutoff, kept because explicitly requested"),
    ("Ca", "mp-21", None, "elemental, on-hull"),
    ("Pb", "mp-20483", None, "elemental, on-hull"),
    ("N2", "mp-1059834", None, "periodic proxy for the molecule; no MP entry for solid N2 is on-hull (~1.1 eV/at above -- expected, real N2 is a gas)"),
    ("O2", "mp-1524462", None, "on-hull; O2 is a real triplet-ground-state (paramagnetic) molecule -- needs a dedicated spin-polarized VASP run, see prepare_extension_vasp_lobster.py"),
    ("PbN6", "mp-667338", None, "lead azide, 84 sites -- much larger than the rest of the dataset, kept because explicitly requested"),
    ("Mn2O7", "mp-28338", None, "Mn(VII), d0 -- expected non-magnetic, standard (non-spin-polarized) run should be fine"),
    ("MnO2", "mp-510408", "P4_2/mnm", "rutile/pyrolusite-type; Mn4+ is d3, physically magnetic (MnO2 is a known antiferromagnet) -- needs the same spin-polarized exception as O2; MP itself reports this structure as is_metal=True (gap=0), a plausible symptom of the well-known GGA self-interaction-error failure mode for correlated Mn oxides -- flagged, not corrected, no +U used here"),
    ("Mn", "mp-35", None, "alpha-Mn, the real room-temperature standard state (on-hull), elemental decomposition reference for Mn2O7/MnO2 -- 29 sites/cell (4 crystallographically distinct Mn sites), a genuinely complex non-collinear antiferromagnet in reality; the spin-polarized exception used here is a simple collinear FM initial guess, a known-crude approximation for this specific element, documented rather than fixed"),
    ("S", "mp-77", None, "alpha-S8 (Fddd), the standard room-temperature reference state for sulfur, elemental decomposition reference for S4N2/S4N4 -- 32 sites/cell, non-magnetic (closed-shell), no exception needed"),
]

# (label, cod_file_id, cif_path, doi, journal, year, note)
COD_SOURCED = [
    (
        "S4N4", "7017102", COD_CIF_DIR / "s4n4_7017102.cif",
        "10.1039/c1dt11418b", "Dalton Transactions", 2011,
        "tetrasulfur tetranitride, P21/n, Z=4 -- no Materials Project entry exists for this composition (verified: zero hits across the entire N-S chemsys)",
    ),
    (
        "S4N2", "4031496", COD_CIF_DIR / "s4n2_4031496.cif",
        None, "J. Chem. Soc., Dalton Trans.", 1981,
        "tetrasulfur dinitride, P4_2nm, Z=4 -- no Materials Project entry exists for this composition",
    ),
]


def download_mp_sourced() -> int:
    mp_ids = [mp_id for _, mp_id, _, _ in MP_SOURCED]
    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=mp_ids,
            fields=[
                "material_id", "structure", "formula_pretty", "energy_above_hull",
                "theoretical", "nsites", "symmetry", "is_metal", "band_gap",
            ],
        )
    by_id = {str(d.material_id): d for d in docs}

    n_ok = 0
    for label, mp_id, requested_sg, note in MP_SOURCED:
        d = by_id.get(mp_id)
        if d is None:
            print(f"MISSING {label} ({mp_id})")
            continue
        actual_sg = d.symmetry.symbol if d.symmetry else None
        if requested_sg and actual_sg != requested_sg:
            print(f"WARNING {label} ({mp_id}): requested sg={requested_sg} but MP reports {actual_sg}")

        compound_dir = OUT_ROOT / f"extension_{label}_{mp_id}"
        compound_dir.mkdir(exist_ok=True)
        d.structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")

        elements = {str(e) for e in d.structure.composition.elements}
        bond_type = classify(elements, d.is_metal)

        meta = {
            "label": label,
            "mp_id": mp_id,
            "formula": d.formula_pretty,
            "family": "extension",
            "source": "materials_project",
            "bond_type": bond_type,
            "is_metal": d.is_metal,
            "band_gap_eV": d.band_gap,
            "energy_above_hull_eV_per_atom": d.energy_above_hull,
            "theoretical": d.theoretical,
            "nsites": d.nsites,
            "spacegroup": actual_sg,
            "requested_spacegroup": requested_sg,
            "note": note,
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {compound_dir.name:35s} sg={actual_sg} nsites={d.nsites} is_metal={d.is_metal}")
        n_ok += 1
    return n_ok


def download_cod_sourced() -> int:
    n_ok = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for label, cod_id, cif_path, doi, journal, year, note in COD_SOURCED:
            parser = CifParser(cif_path)
            structure = parser.parse_structures(primitive=False)[0]
            sg_symbol, sg_number = structure.get_space_group_info()

            compound_dir = OUT_ROOT / f"extension_{label}_cod{cod_id}"
            compound_dir.mkdir(exist_ok=True)
            structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")
            (compound_dir / f"cod_{cod_id}.cif").write_text(Path(cif_path).read_text())

            elements = {str(e) for e in structure.composition.elements}
            # is_metal unknown -- no MP entry to draw on for this composition;
            # left None rather than guessed, to be filled in later from this
            # project's own converged VASP band gap if/when that's computed
            # (still needs the same locally-coarse-k-mesh caution documented
            # in analysis/METRIC_DEFINITION_antibonding.md before being trusted).
            bond_type = classify(elements, is_metal=False)

            meta = {
                "label": label,
                "cod_id": cod_id,
                "formula": structure.composition.reduced_formula,
                "family": "extension",
                "source": "COD",
                "doi": doi,
                "journal": journal,
                "year": year,
                "bond_type": bond_type,
                "is_metal": None,
                "nsites": len(structure),
                "spacegroup": sg_symbol,
                "spacegroup_number": sg_number,
                "note": note,
            }
            (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
            print(f"OK  {compound_dir.name:35s} sg={sg_symbol} nsites={len(structure)} (source: COD {cod_id}, doi={doi})")
            n_ok += 1
    return n_ok


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    n_mp = download_mp_sourced()
    n_cod = download_cod_sourced()
    print(f"\n{n_mp + n_cod}/{len(MP_SOURCED) + len(COD_SOURCED)} extension compounds downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
