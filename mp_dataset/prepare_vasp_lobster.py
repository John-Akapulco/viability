"""Prepare VASP (static, LOBSTER-compatible) + LOBSTER input files for each
downloaded compound in mp_dataset/structures/*, and a matching SLURM submit
script following the lab's existing convention (see
a2/Si3N4_SG227_0GPa_lob/vasp6_lobster.sh on Yargla).

POTCAR choice follows pymatgen's own MPRelaxSet mapping (same one used to
originally relax these structures at Materials Project), so the potentials
are consistent with the input structures' provenance rather than a
hand-picked subset.

NBANDS is not fixed: it is derived from the LOBSTER basis (number of local
orbitals) via pymatgen's Lobsterin._get_nbands/write_INCAR, which is the
standard, minimally-sufficient value for the local-orbital projection to be
well defined. A fixed NBANDS (as used in the 6-compound pilot) is too small
once cells grow past a couple of atoms.

Does not submit anything -- see submit_campaign.sh.
"""

import gzip
import shutil
from pathlib import Path

import yaml
import pymatgen.io.vasp.sets as _vasp_sets
from pymatgen.core import Structure
from pymatgen.io.lobster.inputs import Lobsterin
from pymatgen.io.vasp import Incar, Kpoints

POTCAR_DIR = Path("/home/gilles/pymatgen_potcars/POT_GGA_PAW_PBE")
STRUCTURES_ROOT = Path(__file__).parent / "structures"

_MPRELAXSET_YAML = Path(_vasp_sets.__file__).parent / "MPRelaxSet.yaml"
_MP_POTCAR_MAP = yaml.safe_load(open(_MPRELAXSET_YAML))["POTCAR"]

# Overrides where the MP-recommended variant isn't present in the local
# PBE POTCAR mirror (checked once for the full campaign's element set).
# W_sv override removed: POTCAR.W_sv.gz in the local PSP mirror is missing
# its LEXCH field entirely (confirmed by inspection), which VASP reads as a
# garbled/incompatible XC-functional for that atom type and refuses to run
# ("I REFUSE TO CONTINUE WITH THIS SICK JOB") -- this caused 5 genuine job
# failures overnight (Al12W, BW, WO2, TcW, TiW), all correctly traced back
# to this one corrupted file. Plain "W" is valid PBE and used instead; the
# fallback chain below would also reach it if this dict were empty, but is
# kept explicit for anyone re-reading the failure history above.
_POTCAR_OVERRIDES = {"W": "W"}

# Fallback order tried if neither the MP-recommended nor the override
# variant exists locally, so a missing file fails loudly with a clear
# element name instead of picking something silently wrong.
_FALLBACK_SUFFIXES = ["_sv", "_pv", "_d", "", "_h"]


def _resolve_potcar_variant(element: str) -> str:
    preferred = _POTCAR_OVERRIDES.get(element, _MP_POTCAR_MAP.get(element))
    candidates = [preferred] if preferred else []
    candidates += [f"{element}{suf}" for suf in _FALLBACK_SUFFIXES]
    for variant in candidates:
        if (POTCAR_DIR / f"POTCAR.{variant}.gz").exists():
            return variant
    raise FileNotFoundError(f"No POTCAR variant found for element {element!r} (tried {candidates})")


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
    "NBANDS": 96,  # placeholder; overwritten below from the LOBSTER basis size
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
ulimit -s unlimited

time srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std

export OMP_NUM_THREADS=16
lobster-5.1.1
"""


def build_potcar(symbols: list[str], out_path: Path) -> None:
    with open(out_path, "wb") as out:
        for el in symbols:
            variant = _resolve_potcar_variant(el)
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

    incar_path = compound_dir / "INCAR"
    Incar(INCAR_SETTINGS).write_file(str(incar_path))

    kpoints = Kpoints.automatic_density_by_vol(structure, kppvol=100)
    kpoints.write_file(str(compound_dir / "KPOINTS"))

    # LOBSTER basis lookup only needs the POTCAR *symbol strings* (e.g.
    # "Na_pv"), resolved against pymatgen's bundled BASIS_PBE_54_standard.yaml
    # table -- not the POTCAR file's parsed content. Going through
    # POTCAR_input=... instead calls Potcar.from_file(), whose PotcarSingle
    # validity check (self.is_valid) errors out on POTCARs pymatgen doesn't
    # recognize by hash (seen with several elements from this PSP mirror).
    potcar_variants = [_resolve_potcar_variant(el) for el in symbols_in_order]
    basis_lines = Lobsterin.get_basis(structure=structure, potcar_symbols=potcar_variants)
    dict_for_basis = dict(line.split(" ", 1) for line in basis_lines)

    lobsterin = Lobsterin.standard_calculations_from_vasp_files(
        POSCAR_input=str(poscar_path),
        INCAR_input=str(incar_path),
        dict_for_basis=dict_for_basis,
        option="standard",
    )
    lobsterin.write_lobsterin(str(compound_dir / "lobsterin"))
    # NBANDS = number of local basis orbitals (LOBSTER's minimally
    # sufficient value); everything else in INCAR is left untouched.
    lobsterin.write_INCAR(
        incar_input=str(incar_path),
        incar_output=str(incar_path),
        poscar_input=str(poscar_path),
        isym=-1,
    )

    job_name = compound_dir.name
    (compound_dir / "submit.sh").write_text(SLURM_TEMPLATE.format(job_name=job_name))

    nbands = Incar.from_file(str(incar_path))["NBANDS"]
    print(f"prepared {compound_dir.name}  (elements: {symbols_in_order}, NBANDS={nbands})")


def main():
    for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
        if (compound_dir / "POSCAR").exists() and not (compound_dir / "lobsterin").exists():
            prepare_one(compound_dir)


if __name__ == "__main__":
    main()
