"""Prepare VASP+LOBSTER inputs for the mp_dataset/structures/extension_*
compounds (download_extension.py), reusing prepare_vasp_lobster.py's POTCAR
resolution, base INCAR settings, and SLURM template UNCHANGED -- the only
difference is a spin-polarized (ISPIN=2) exception for the compounds
flagged in EXTENSION_SPIN_OVERRIDES below, since the rest of this project's
pipeline (INCAR_SETTINGS) assumes non-magnetic compounds and several of the
extension compounds are physically magnetic (O2 triplet; MnO2 and Mn have
partially-filled Mn d-shells; Fe/Co/Ni are ferromagnetic and Cr is an
antiferromagnetic/SDW metal among the download_elements_reference.py
elemental references).

Initial MAGMOM guesses are simple collinear ferromagnetic starting points
(VASP relaxes them self-consistently) -- explicitly NOT an attempt at the
correct magnetic ordering (e.g. MnO2's real antiferromagnetism, alpha-Mn's
genuinely non-collinear structure, or Cr's real SDW incommensurate order).
Documented as a known approximation in mp_metadata.json already; repeated
here at the point where it actually takes effect.

Does not submit anything.
"""

import sys
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.lobster.inputs import Lobsterin
from pymatgen.io.vasp import Incar, Kpoints

sys.path.insert(0, str(Path(__file__).parent))
from prepare_vasp_lobster import INCAR_SETTINGS, SLURM_TEMPLATE, build_potcar  # noqa: E402  (reused unmodified)

STRUCTURES_ROOT = Path(__file__).parent / "structures"

# element -> initial MAGMOM (mu_B/atom); elements not listed get 0.0
EXTENSION_SPIN_OVERRIDES = {
    "extension_O2_mp-1524462": {"O": 1.5},
    "extension_MnO2_mp-510408": {"Mn": 3.0, "O": 0.0},
    "extension_Mn_mp-35": {"Mn": 3.0},
    "extension_Fe_mp-13": {"Fe": 2.2},
    "extension_Co_mp-102": {"Co": 1.7},
    "extension_Ni_mp-23": {"Ni": 0.6},
    "extension_Cr_mp-90": {"Cr": 1.0},
    # download_extension4.py: alkali superoxides (open-shell O2^- radical,
    # same physical class as the O2 molecule override above) -- SrO2 in
    # the same batch is the alkaline-earth PEROXIDE (closed-shell O2^2-)
    # and deliberately does NOT get an exception here.
    "extension_CsO2_exp_mp-1441": {"O": 1.0},
    "extension_CsO2_theo_mp-1096936": {"O": 1.0},
}


def magmom_for(structure: Structure, per_element: dict) -> list[float]:
    return [per_element.get(str(site.specie), 0.0) for site in structure]


def prepare_one(compound_dir: Path) -> None:
    poscar_path = compound_dir / "POSCAR"
    structure = Structure.from_file(poscar_path)
    symbols_in_order = list(dict.fromkeys(str(s) for s in structure.species))

    potcar_path = compound_dir / "POTCAR"
    build_potcar(symbols_in_order, potcar_path)

    incar_settings = dict(INCAR_SETTINGS)
    spin_override = EXTENSION_SPIN_OVERRIDES.get(compound_dir.name)
    if spin_override is not None:
        incar_settings["ISPIN"] = 2
        incar_settings["MAGMOM"] = magmom_for(structure, spin_override)

    incar_path = compound_dir / "INCAR"
    Incar(incar_settings).write_file(str(incar_path))

    kpoints = Kpoints.automatic_density_by_vol(structure, kppvol=100)
    kpoints.write_file(str(compound_dir / "KPOINTS"))

    potcar_variants = [
        __import__("prepare_vasp_lobster")._resolve_potcar_variant(el) for el in symbols_in_order
    ]
    basis_lines = Lobsterin.get_basis(structure=structure, potcar_symbols=potcar_variants)
    dict_for_basis = dict(line.split(" ", 1) for line in basis_lines)

    lobsterin = Lobsterin.standard_calculations_from_vasp_files(
        POSCAR_input=str(poscar_path),
        INCAR_input=str(incar_path),
        dict_for_basis=dict_for_basis,
        option="standard",
    )
    lobsterin.write_lobsterin(str(compound_dir / "lobsterin"))
    lobsterin.write_INCAR(
        incar_input=str(incar_path),
        incar_output=str(incar_path),
        poscar_input=str(poscar_path),
        isym=-1,
    )

    job_name = compound_dir.name
    (compound_dir / "submit.sh").write_text(SLURM_TEMPLATE.format(job_name=job_name))

    final_incar = Incar.from_file(str(incar_path))
    spin_note = f" ISPIN=2 MAGMOM={final_incar.get('MAGMOM')}" if spin_override else ""
    print(f"prepared {compound_dir.name}  (elements: {symbols_in_order}, NBANDS={final_incar['NBANDS']}{spin_note})")


def main():
    for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
        if compound_dir.name.startswith("extension_") and (compound_dir / "POSCAR").exists() and not (compound_dir / "lobsterin").exists():
            prepare_one(compound_dir)


if __name__ == "__main__":
    main()
