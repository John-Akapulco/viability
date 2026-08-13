"""Download the 6 curated MP structures (2 families x 3 bonding types) as
POSCAR files, one directory per compound, ready for a LOBSTER-compatible
VASP static run.
"""

import os
from pathlib import Path

from mp_api.client import MPRester

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

SELECTION = [
    # (family, bond_type, mp_id, label)
    ("hull", "ionic", "mp-22862", "NaCl"),
    ("hull", "covalent", "mp-149", "Si"),
    ("hull", "metallic", "mp-1487", "AlNi"),
    ("metastable", "ionic", "mp-23259", "LiBr"),
    ("metastable", "covalent", "mp-169", "C_rhombohedral"),
    ("metastable", "metallic", "mp-2323", "BeCu"),
]

OUT_ROOT = Path(__file__).parent / "structures"


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    with MPRester(API_KEY) as mpr:
        for family, bond_type, mp_id, label in SELECTION:
            doc = mpr.materials.summary.search(
                material_ids=[mp_id],
                fields=["material_id", "structure", "energy_above_hull", "formula_pretty"],
            )[0]
            structure = doc.structure
            compound_dir = OUT_ROOT / f"{family}_{bond_type}_{label}_{mp_id}"
            compound_dir.mkdir(exist_ok=True)
            structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")
            meta = {
                "material_id": mp_id,
                "formula": doc.formula_pretty,
                "family": family,
                "bond_type": bond_type,
                "energy_above_hull_eV_per_atom": doc.energy_above_hull,
                "num_sites": len(structure),
            }
            with open(compound_dir / "mp_metadata.json", "w") as f:
                import json

                json.dump(meta, f, indent=2)
            print(f"OK  {compound_dir.name:45s} ({len(structure)} sites)")


if __name__ == "__main__":
    main()
