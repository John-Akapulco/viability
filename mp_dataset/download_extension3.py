"""Third extension batch: a single compound, elemental Ti (mp-46), the
missing decomposition-into-elements reference for the 3 TiO2 high-pressure
polymorphs already in the extension set (extension_TiO2II_mp-1439,
extension_TiO2baddeleyite_mp-430, extension_TiO2OII_mp-9173) -- needed for
analysis/reaction_icohp.py case 1 (compound -> elements).

mp-72 is MP's own on-hull Ti entry (energy_above_hull=0.0) but is the
omega-Ti phase (P6/mmm, 3 sites/cell, AlB2-type) -- a well-documented
PBE-GGA artifact where the omega phase is spuriously favored over the real
room-temperature ground state by a small margin. mp-46 (hcp alpha-Ti,
P6_3/mmc, 2 sites/cell, 0.0152 eV/at above MP's own hull) is the real
experimental standard state, used here instead -- same "real standard
state over literal DFT on-hull" convention already used for Mn
(alpha-Mn, mp-35, download_extension.py) and Sulfur (alpha-S8, mp-77).
Flagged, not corrected, consistent with how this project treats other
known GGA artifacts (e.g. the MnO2 is_metal self-interaction-error note).

Writes one directory under mp_dataset/structures/, same layout as
download_extension.py / download_extension2.py, family="extension".
"""

import json
import os
import sys
from pathlib import Path

from mp_api.client import MPRester

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from fetch_candidates import classify  # noqa: E402  (reused unmodified)

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()
OUT_ROOT = Path(__file__).parent / "structures"

LABEL = "Ti"
MP_ID = "mp-46"
EXPECTED_SG = "P6_3/mmc"
NOTE = (
    "hcp alpha-Ti, the real room-temperature standard state; 0.0152 eV/at "
    "above MP's own hull -- MP's literal on-hull Ti entry (mp-72, "
    "energy_above_hull=0.0) is the omega-Ti phase (P6/mmm, AlB2-type), a "
    "well-documented PBE-GGA artifact that spuriously favors omega over "
    "alpha by a small margin; alpha-Ti used instead as the elemental "
    "decomposition reference for the TiO2 high-pressure polymorphs "
    "(extension_TiO2II/baddeleyite/OII), same convention as the alpha-Mn "
    "and alpha-S8 references (flagged, not corrected)."
)


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=[MP_ID],
            fields=[
                "material_id", "structure", "formula_pretty", "energy_above_hull",
                "theoretical", "nsites", "symmetry", "is_metal", "band_gap",
            ],
        )
    if not docs:
        print(f"MISSING {LABEL} ({MP_ID})")
        return
    d = docs[0]
    actual_sg = d.symmetry.symbol if d.symmetry else None
    if actual_sg != EXPECTED_SG:
        print(f"WARNING {LABEL} ({MP_ID}): expected sg={EXPECTED_SG} but MP reports {actual_sg}")

    compound_dir = OUT_ROOT / f"extension_{LABEL}_{MP_ID}"
    compound_dir.mkdir(exist_ok=True)
    d.structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")

    elements = {str(e) for e in d.structure.composition.elements}
    bond_type = classify(elements, d.is_metal)

    meta = {
        "label": LABEL,
        "mp_id": MP_ID,
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
        "expected_spacegroup": EXPECTED_SG,
        "note": NOTE,
    }
    (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
    print(f"OK  {compound_dir.name:35s} sg={actual_sg} nsites={d.nsites} "
          f"EAH={d.energy_above_hull:.4f} is_metal={d.is_metal}")


if __name__ == "__main__":
    main()
