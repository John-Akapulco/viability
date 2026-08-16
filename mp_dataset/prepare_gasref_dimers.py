"""Build and relax isolated gas-phase dimer references (H2, N2, O2, F2,
Cl2) as a direct test of the user's concern (2026-08-17): the project's
current ELEMENT_REFERENCE entries for these five elements
(analysis/compute_reaction_icohp_case1.py) are periodic MP structures
described as solid molecular crystals of the frozen gas, not isolated
gas-phase molecules -- and DFT-PBE is known to evaluate diatomic gas
reference energies inconsistently (systematic errors severe enough that
Materials Project itself applies compatibility corrections to
O2/N2/F2/H2 references when computing formation_energy_per_atom).

Direct structural check performed before writing this script (see
project memory / session notes): four of the five current references
(H, O, F, Cl) ARE genuine isolated-dimer molecular crystals -- each atom
has exactly one short contact, matching gas-phase bond lengths almost
exactly (H2 0.741 vs lit. 0.741 A; O2 1.222 vs 1.208; F2 1.414 vs 1.412;
Cl2 1.997 vs 1.988). But the N reference (extension_N2_mp-1059834,
theoretical=True, is_metal=True) is NOT diatomic at all: each N atom has
TWO short N-N contacts at 1.296 A (ICOHP ~-15.9 eV each) -- a polymeric/
network nitrogen phase, not molecular N2 (gas-phase N-N is ~1.10 A). Any
case-1 reaction (analysis/compute_reaction_icohp_case1.py) decomposing a
nitrogen-containing compound into elements currently uses this wrong
phase's ICOHP as the N2 product-side reference.

Methodology, per user's explicit specification: model each gas element
as a single dimer molecule in a large (10x10x10 A), UNoptimized cubic
box -- only the interatomic distance is relaxed (IBRION=2, ISIF=2, cell
fixed), not the box itself, since there is no periodic image
interaction to relax against at this size. Starting distances are
literature gas-phase bond lengths (H2 0.74, N2 1.10, O2 1.21, F2 1.41,
Cl2 1.99 A); O2 gets the same ISPIN=2/MAGMOM=1.5 triplet-ground-state
treatment already used for the periodic O2 reference
(prepare_extension_vasp_lobster.EXTENSION_SPIN_OVERRIDES) -- getting
O2's spin state right is itself part of the DFT-reference-accuracy
story the user is raising.

Two-stage pipeline (LOBSTER needs a converged STATIC wavefunction at
the FINAL geometry, same reasoning as every other compound in this
project): this script prepares and submits ONLY the stage-1 relaxation
(vasp_std only, no LOBSTER). Directories are gasref_<EL>2_dimerbox/.
prepare_gasref_stage2.py (run after stage 1 converges) copies CONTCAR
back to POSCAR and hands the directory to
prepare_extension_vasp_lobster.prepare_one() unchanged -- the exact same
static+LOBSTER convention used for every MP-provided pre-relaxed
structure in this project.

Exploratory/diagnostic only: NOT wired into ELEMENT_REFERENCE or any
analysis pipeline. Results are for comparison against the current
periodic references; swapping ELEMENT_REFERENCE itself is a separate,
deliberate follow-up once the magnitude of the discrepancy is known.
"""

import sys
from pathlib import Path

from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Incar, Kpoints

sys.path.insert(0, str(Path(__file__).parent))
from prepare_vasp_lobster import INCAR_SETTINGS, build_potcar  # noqa: E402  (reused unmodified)

STRUCTURES_ROOT = Path(__file__).parent / "structures"

BOX_LENGTH = 10.0  # Angstrom, cubic, fixed (not relaxed)

DIMERS = {
    "H2": {"element": "H", "d0": 0.74},
    "N2": {"element": "N", "d0": 1.10},
    "O2": {"element": "O", "d0": 1.21, "magmom": 1.5},
    "F2": {"element": "F", "d0": 1.41},
    "Cl2": {"element": "Cl", "d0": 1.99},
}

STAGE1_SLURM_TEMPLATE = """#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=vasp.log
#SBATCH --time=2:00:00

module purge
module load intel
module load impi/2021.13
module load vasp/6.5.0
export OMP_NUM_THREADS=1
ulimit -s unlimited

time srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std
"""


def build_dimer_structure(element: str, d0: float) -> Structure:
    lattice = Lattice.cubic(BOX_LENGTH)
    center = BOX_LENGTH / 2
    z1 = center - d0 / 2
    z2 = center + d0 / 2
    return Structure(
        lattice, [element, element],
        [[center, center, z1], [center, center, z2]],
        coords_are_cartesian=True,
    )


def prepare_one(name: str, spec: dict) -> Path:
    element = spec["element"]
    compound_dir = STRUCTURES_ROOT / f"gasref_{name}_dimerbox"
    compound_dir.mkdir(exist_ok=True)

    structure = build_dimer_structure(element, spec["d0"])
    structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")

    build_potcar([element], compound_dir / "POTCAR")

    incar_settings = dict(INCAR_SETTINGS)
    incar_settings.update({
        "IBRION": 2, "NSW": 50, "ISIF": 2, "EDIFFG": -0.01,
        "LWAVE": False, "LCHARG": False, "ISYM": 0,
        "NBANDS": 20,  # generous for a 2-atom cell; not LOBSTER-basis-derived at this stage
        "KPAR": 1,  # single (Gamma) k-point -- KPAR>1 has nothing to split
    })
    if "magmom" in spec:
        incar_settings["ISPIN"] = 2
        incar_settings["MAGMOM"] = [spec["magmom"], spec["magmom"]]
    Incar(incar_settings).write_file(str(compound_dir / "INCAR"))

    Kpoints.gamma_automatic(kpts=(1, 1, 1)).write_file(str(compound_dir / "KPOINTS"))

    (compound_dir / "submit.sh").write_text(STAGE1_SLURM_TEMPLATE.format(job_name=f"gasref_{name}"))

    (compound_dir / "gasref_meta.json").write_text(
        f'{{"element": "{element}", "d0_angstrom": {spec["d0"]}, "box_length_angstrom": {BOX_LENGTH}, '
        f'"stage": 1}}\n'
    )
    print(f"prepared {compound_dir.name} (d0={spec['d0']} A, box={BOX_LENGTH}A)")
    return compound_dir


def main():
    STRUCTURES_ROOT.mkdir(exist_ok=True)
    for name, spec in DIMERS.items():
        prepare_one(name, spec)


if __name__ == "__main__":
    main()
