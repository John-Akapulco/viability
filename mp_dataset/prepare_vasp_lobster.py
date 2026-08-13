"""Prepare VASP (static, LOBSTER-compatible) + LOBSTER input files for each
downloaded compound in mp_dataset/structures/*, and a matching SLURM submit
script following the lab's existing convention (see
a2/Si3N4_SG227_0GPa_lob/vasp6_lobster.sh on Yargla).

Does not submit anything -- see submit_all.sh.
"""

import gzip
import shutil
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.lobster.inputs import Lobsterin
from pymatgen.io.vasp import Incar, Kpoints, Potcar

POTCAR_DIR = Path("/home/gilles/pymatgen_potcars/POT_GGA_PAW_PBE")
STRUCTURES_ROOT = Path(__file__).parent / "structures"

# VASP recommended PAW-PBE potentials for the elements in this dataset.
RECOMMENDED_POTCAR = {
    "Na": "Na_pv",
    "Cl": "Cl",
    "Si": "Si",
    "Al": "Al",
    "Ni": "Ni",
    "Li": "Li_sv",
    "Br": "Br",
    "C": "C",
    "Be": "Be",
    "Cu": "Cu",
}

INCAR_SETTINGS = {
    "PREC": "Accurate",
    "LASPH": True,
    "EDIFF": 1e-6,
    "ENCUT": 600,
    "IBRION": -1,
    "NSW": 0,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "LWAVE": True,
    "LCHARG": True,
    "ISYM": -1,
    "NPAR": 4,
    "KPAR": 2,
    "NBANDS": 100,
}

SLURM_TEMPLATE = """#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=vasp.log
#SBATCH --time=24:00:00

module purge
module load intel
module load impi/2021.13
module load vasp/6.5.0
module load lobster/5.1.1
export OMP_NUM_THREADS=1

time srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std

export OMP_NUM_THREADS=16
lobster-5.1.1
"""


def build_potcar(symbols: list[str], out_path: Path) -> None:
    with open(out_path, "wb") as out:
        for el in symbols:
            variant = RECOMMENDED_POTCAR[el]
            gz_path = POTCAR_DIR / f"POTCAR.{variant}.gz"
            with gzip.open(gz_path, "rb") as f:
                shutil.copyfileobj(f, out)


def prepare_one(compound_dir: Path) -> None:
    poscar_path = compound_dir / "POSCAR"
    structure = Structure.from_file(poscar_path)
    # element symbols in POSCAR block order (pymatgen preserves POSCAR order)
    symbols_in_order = list(dict.fromkeys(str(s) for s in structure.species))

    potcar_path = compound_dir / "POTCAR"
    build_potcar(symbols_in_order, potcar_path)

    incar = Incar(INCAR_SETTINGS)
    incar.write_file(str(compound_dir / "INCAR"))

    kpoints = Kpoints.automatic_density_by_vol(structure, kppvol=100)
    kpoints.write_file(str(compound_dir / "KPOINTS"))

    lobsterin = Lobsterin.standard_calculations_from_vasp_files(
        POSCAR_input=str(poscar_path),
        INCAR_input=str(compound_dir / "INCAR"),
        POTCAR_input=str(potcar_path),
        option="standard",
    )
    lobsterin.write_lobsterin(str(compound_dir / "lobsterin"))

    job_name = compound_dir.name
    (compound_dir / "submit.sh").write_text(SLURM_TEMPLATE.format(job_name=job_name))

    print(f"prepared {compound_dir.name}  (elements: {symbols_in_order})")


def main():
    for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
        if (compound_dir / "POSCAR").exists():
            prepare_one(compound_dir)


if __name__ == "__main__":
    main()
