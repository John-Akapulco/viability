"""Generate the conventional standard cell for each of the 6 pilot
compounds (from their relaxed CONTCAR, not the raw MP structure), and
prepare VASP+LOBSTER inputs for it via prepare_vasp_lobster.py's existing
prepare_one() -- reused unchanged, not reimplemented.

Writes mp_dataset/structures_conventional_pilot/{compound_id}/. Does not
touch mp_dataset/structures/ (the primitive-cell results stay intact for
comparison).
"""

import json
import sys
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sys.path.insert(0, str(Path(__file__).parent))
import prepare_vasp_lobster as pvl  # noqa: E402  (reused unchanged)

PILOT_COMPOUNDS = [
    "hull_ionic_NaCl_mp-22862",
    "hull_covalent_Si_mp-149",
    "hull_metallic_AlNi_mp-1487",
    "metastable_ionic_LiBr_mp-23259",
    "metastable_covalent_C_rhombohedral_mp-169",
    "metastable_metallic_BeCu_mp-2323",
]

SRC_ROOT = Path(__file__).parent / "structures"
OUT_ROOT = Path(__file__).parent / "structures_conventional_pilot"


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    for compound_id in PILOT_COMPOUNDS:
        src_dir = SRC_ROOT / compound_id
        contcar = src_dir / "CONTCAR"
        if not contcar.exists():
            print(f"SKIP {compound_id}: no CONTCAR found at {contcar}")
            continue

        primitive = Structure.from_file(str(contcar))
        sga = SpacegroupAnalyzer(primitive)
        conventional = sga.get_conventional_standard_structure()

        out_dir = OUT_ROOT / compound_id
        out_dir.mkdir(exist_ok=True)
        conventional.to(filename=str(out_dir / "POSCAR"), fmt="poscar")

        src_meta = json.loads((src_dir / "mp_metadata.json").read_text())
        meta = {
            **src_meta,
            "cell_type": "conventional_standard",
            "n_sites_primitive": len(primitive),
            "n_sites_conventional": len(conventional),
            "expansion_factor": len(conventional) / len(primitive),
            "spacegroup_symbol_sga": sga.get_space_group_symbol(),
            "spacegroup_number_sga": sga.get_space_group_number(),
        }
        (out_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))

        pvl.prepare_one(out_dir)
        print(
            f"{compound_id}: {len(primitive)} -> {len(conventional)} sites "
            f"(x{meta['expansion_factor']:.2f}), SG {meta['spacegroup_symbol_sga']}"
        )


if __name__ == "__main__":
    main()
