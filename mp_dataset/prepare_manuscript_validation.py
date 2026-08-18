"""Prepare a new, independent VASP+LOBSTER campaign that reproduces
Reitz & Dronskowski's own DFT methodology (ic-2026-04181q) as closely as
possible, instead of this project's own default settings
(prepare_vasp_lobster.py: plain PBE, ENCUT=600, EDIFF=1e-6, no D3).

Manuscript methodology (Computational Details + Table S1 of the SI,
~/ic-2026-04181q_Proof_hi.pdf and ~/Above_Convex_Hull_SI_InorgChem.pdf):
- VASP 6.4.1 (this project uses 6.5.0, close enough -- not independently
  reproducible without their exact binary/build).
- PAW pseudopotentials (Blochl) -- same POTCAR files this project already
  uses (functional choice is an INCAR tag, not a different POTCAR set).
- PBEsol functional (GGA=PS) + D3 dispersion correction for every phase
  EXCEPT Mn2O7/MnO2 -- D3 damping variant not stated explicitly in the
  extracted text; BJ-damping (IVDW=12) used per user instruction (modern
  standard, most common today).
- Mn2O7/MnO2: PBEsol+U (Hubbard U=7.09 eV on Mn 3d, LDAUTYPE=2) relaxation,
  then a SEPARATE single-point total-energy calculation with the meta-GGA
  r2SCAN functional (no D3, no +U per the manuscript's literal wording --
  an assumption, not independently confirmed) on the same relaxed geometry.
  LOBSTER/ICOHP still comes from the PBEsol+U wavefunction, not r2SCAN
  (LOBSTER's own basis fit, pbeVASPfit2015, is a PBE-family fit; r2SCAN is
  used here ONLY for the total-energy (Delta E) reproduction, per the
  manuscript's own framing of ICOHP as "calculated ... using exactly the
  same algorithm for all compounds" -- i.e. NOT tied to whichever
  functional produced the total energy).
- MnO2's AFM ordering (Noda et al., not independently retrieved): a
  standard checkerboard AFM on the two symmetry-distinct Mn sites of the
  rutile cell (MAGMOM = 3.0 -3.0 0 0 0 0), replacing this project's
  existing ferromagnetic default (2*3.0 4*0.0) for this one directory.
- ENCUT=700 eV, EDIFF=1e-8 eV, EDIFFG=-1e-6 eV/Angstrom, tetrahedron
  method with Blochl corrections (ISMEAR=-5) for Brillouin-zone
  integration, Monkhorst-Pack k-point meshes exactly as given per-phase
  in Table S1 (not this project's own automatic_density_by_vol scheme).
- N2/O2: isolated molecules in an 8x8x8 Angstrom box (Table S1), not this
  project's existing 10x10x10 gasref_* boxes -- rebuilt here at the
  manuscript's exact box size, same starting bond length.

Starting geometries are NOT hand-built from Table S1's lattice parameters
alone (no atomic coordinates are given there, only cell dimensions) --
each phase reuses this project's own already-verified starting structure
(mp_id/COD-sourced, see report/appendix_reitz_dronskowski_structures_*.tex)
as the pre-relaxation geometry, then lets VASP relax it under the
manuscript's own methodology. Comparing our own converged a/b/c/volume
against Table S1's reported values is the actual validation test this
campaign is for -- copying their numbers into our POSCAR would defeat
that purpose.

Two-VASP-step-plus-LOBSTER pipeline per directory (three steps for
Mn2O7/MnO2), driven by submit.sh, matching this project's established
convention of one sequential SLURM job per compound:
  1. relax  (INCAR.relax: ISIF=3 for real crystals, ISIF=2/fixed-box for
     the two isolated molecules; IBRION=2, NSW=200)
  2. static (INCAR.static: NSW=0, LWAVE/LCHARG=True, ISYM=-1 -- the
     LOBSTER-compatible wavefunction), then lobsterin/LOBSTER
  3. (Mn2O7/MnO2 only) r2scan (INCAR.r2scan: single point, METAGGA=R2SCAN,
     no D3/+U, on the same relaxed CONTCAR) for Delta E_r2SCAN only.

Writes mp_dataset/structures/manuscript_<name>/ for each of the 16 phases.
Does not submit -- see submit_manuscript_validation.sh.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import yaml
import pymatgen.io.vasp.sets as _vasp_sets
from pymatgen.core import Structure
from pymatgen.io.lobster.inputs import Lobsterin
from pymatgen.io.vasp import Incar, Kpoints

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from prepare_vasp_lobster import build_potcar, _resolve_potcar_variant  # noqa: E402  (reused unmodified)

STRUCTURES_ROOT = Path(__file__).parent / "structures"

# name -> (source_dir, kpoints (a,b,c), workflow, relax_mode, extra)
# workflow: "standard" (PBEsol+D3) or "mn" (PBEsol+U then r2SCAN single point)
# relax_mode: "full" (ISIF=3, cell+ions) or "molecule" (ISIF=2, ions only,
#   fixed 8x8x8 box per Table S1)
PHASES = {
    "PbN3_2": dict(source="extension_PbN6_mp-620058", kpoints=(9, 5, 3), workflow="standard", relax="full"),
    "Pb": dict(source="extension_Pb_mp-20483", kpoints=(13, 13, 13), workflow="standard", relax="full"),
    "N2": dict(source="gasref_N2_dimerbox", kpoints=(7, 7, 7), workflow="standard", relax="molecule", box=8.0),
    "S4N2": dict(source="extension_S4N2_cod4031496", kpoints=(5, 5, 15), workflow="standard", relax="full"),
    "S8": dict(source="extension_S_mp-77", kpoints=(5, 5, 3), workflow="standard", relax="full"),
    "S4N4": dict(source="extension_S4N4_cod7017102", kpoints=(7, 9, 7), workflow="standard", relax="full"),
    "ZnSn": dict(source="gasref_ZnSn_NiAs", kpoints=(17, 17, 11), workflow="standard", relax="full"),
    "Zn": dict(source="extension_Zn_mp-79", kpoints=(23, 23, 13), workflow="standard", relax="full"),
    "Sn": dict(source="extension_Snalpha_mp-117", kpoints=(9, 9, 9), workflow="standard", relax="full"),
    "CaO_sphalerite": dict(source="gasref_CaO_sphalerite", kpoints=(11, 11, 11), workflow="standard", relax="full"),
    "CaO_rocksalt": dict(source="extension_CaO_mp-2605", kpoints=(13, 13, 13), workflow="standard", relax="full"),
    "CaN": dict(source="extension_CaN_mp-1058549", kpoints=(11, 11, 11), workflow="standard", relax="full"),
    "Ca3N2": dict(source="extension_Ca3N2_mp-844", kpoints=(5, 5, 5), workflow="standard", relax="full"),
    "Mn2O7": dict(source="extension_Mn2O7_mp-28338", kpoints=(9, 3, 7), workflow="mn", relax="full"),
    "MnO2": dict(source="extension_MnO2_mp-510408", kpoints=(13, 13, 21), workflow="mn", relax="full", afm=True),
    "O2": dict(source="gasref_O2_dimerbox", kpoints=(7, 7, 7), workflow="standard", relax="molecule", box=8.0, spin=True),
}

BASE_TAGS = {
    "PREC": "Accurate",
    "LASPH": True,
    "ENCUT": 700,
    "EDIFF": 1e-8,
    "ISYM": -1,
    "NPAR": 4,
    "KPAR": 2,
}

RELAX_TAGS = {
    "IBRION": 2,
    "NSW": 200,
    "ISMEAR": -5,
    "EDIFFG": -1e-6,
    "LWAVE": False,
    "LCHARG": False,
}

STATIC_TAGS = {
    "IBRION": -1,
    "NSW": 0,
    "ISMEAR": -5,
    "LWAVE": True,
    "LCHARG": True,
}

R2SCAN_TAGS = {
    "IBRION": -1,
    "NSW": 0,
    "ISMEAR": -5,
    "METAGGA": "R2SCAN",
    "LWAVE": False,
    "LCHARG": False,
}

D3_TAGS = {"GGA": "PS", "IVDW": 12}
MN_U_TAGS = {
    "GGA": "PS",
    "LDAU": True,
    "LDAUTYPE": 2,
    "LDAUL": None,  # set per-structure below (2 for Mn, -1 for O)
    "LDAUU": None,
    "LDAUJ": None,
    "LDAUPRINT": 1,
    "ISPIN": 2,
}

SLURM_TEMPLATE_STANDARD = """#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=vasp.log
#SBATCH --time=48:00:00

module purge
module load intel
module load impi/2021.13
module load vasp/6.5.0
module load lobster/5.1.1
export OMP_NUM_THREADS=1
ulimit -s unlimited

cp INCAR.relax INCAR
{{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; }} 2>> vasp.log
cp CONTCAR POSCAR
cp CONTCAR CONTCAR.relaxed

cp INCAR.static INCAR
{{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; }} 2>> vasp.log

export OMP_NUM_THREADS=16
{{ time -p lobster-5.1.1 ; }} 2> lobster_time.txt
"""

SLURM_TEMPLATE_MN_EXTRA = """
export OMP_NUM_THREADS=1
cp INCAR.r2scan INCAR
{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; } 2>> vasp.log
cp OSZICAR OSZICAR.r2scan
cp OUTCAR OUTCAR.r2scan
"""


def resize_molecule_box(structure: Structure, box: float) -> Structure:
    """Keep atoms' Cartesian positions (bond length) fixed, resize the
    lattice to an box x box x box cube -- matches Table S1's box convention
    without changing the molecule's own starting geometry."""
    cart_coords = structure.cart_coords
    new_lattice = [[box, 0, 0], [0, box, 0], [0, 0, box]]
    return Structure(new_lattice, structure.species, cart_coords, coords_are_cartesian=True)


def write_incar(path: Path, tags: dict) -> None:
    clean = {k: v for k, v in tags.items() if v is not None}
    Incar(clean).write_file(str(path))


def prepare_one(name: str, spec: dict) -> None:
    source_dir = STRUCTURES_ROOT / spec["source"]
    out_dir = STRUCTURES_ROOT / f"manuscript_{name}"
    out_dir.mkdir(exist_ok=True)

    structure = Structure.from_file(source_dir / "POSCAR")
    if spec["relax"] == "molecule":
        structure = resize_molecule_box(structure, spec["box"])
    structure.to(filename=str(out_dir / "POSCAR"), fmt="poscar")

    symbols_in_order = list(dict.fromkeys(str(s) for s in structure.species))
    build_potcar(symbols_in_order, out_dir / "POTCAR")

    kx, ky, kz = spec["kpoints"]
    kpoints = Kpoints.gamma_automatic((kx, ky, kz))
    kpoints.write_file(str(out_dir / "KPOINTS"))

    incar_relax = dict(BASE_TAGS)
    incar_relax.update(RELAX_TAGS)
    incar_static = dict(BASE_TAGS)
    incar_static.update(STATIC_TAGS)

    if spec["relax"] == "molecule":
        incar_relax["ISIF"] = 2  # ions only, box fixed at the manuscript's 8x8x8
    else:
        incar_relax["ISIF"] = 3  # full cell + ion relaxation

    if spec.get("spin"):
        incar_relax["ISPIN"] = 2
        incar_static["ISPIN"] = 2
        incar_relax["MAGMOM"] = f"{len(structure)}*1.0"
        incar_static["MAGMOM"] = f"{len(structure)}*1.0"

    if spec["workflow"] == "standard":
        incar_relax.update(D3_TAGS)
        incar_static.update(D3_TAGS)
    elif spec["workflow"] == "mn":
        mn_tags = dict(MN_U_TAGS)
        n_species = len(structure.symbol_set)
        # LDAUL/LDAUU/LDAUJ ordered per POSCAR species block; Mn gets L=2
        # (d), U=7.09, J=0; every other species (O) gets L=-1 (no +U).
        ldaul, ldauu, ldauj = [], [], []
        for el in symbols_in_order:
            if el == "Mn":
                ldaul.append(2)
                ldauu.append(7.09)
                ldauj.append(0)
            else:
                ldaul.append(-1)
                ldauu.append(0)
                ldauj.append(0)
        mn_tags["LDAUL"] = " ".join(str(x) for x in ldaul)
        mn_tags["LDAUU"] = " ".join(str(x) for x in ldauu)
        mn_tags["LDAUJ"] = " ".join(str(x) for x in ldauj)
        incar_relax.update(mn_tags)
        incar_static.update(mn_tags)
        if spec.get("afm"):
            # rutile MnO2: 2 Mn (checkerboard AFM) + 4 O (nonmagnetic).
            n_o = sum(1 for s in structure.species if str(s) == "O")
            incar_relax["MAGMOM"] = f"3.0 -3.0 {n_o}*0.0"
            incar_static["MAGMOM"] = f"3.0 -3.0 {n_o}*0.0"
        else:
            n_mn = sum(1 for s in structure.species if str(s) == "Mn")
            n_o = sum(1 for s in structure.species if str(s) == "O")
            incar_relax["MAGMOM"] = f"{n_mn}*3.0 {n_o}*0.0" if n_mn else None
            incar_static["MAGMOM"] = incar_relax.get("MAGMOM")

    write_incar(out_dir / "INCAR.relax", incar_relax)
    write_incar(out_dir / "INCAR.static", incar_static)

    if spec["workflow"] == "mn":
        incar_r2scan = dict(BASE_TAGS)
        incar_r2scan.update(R2SCAN_TAGS)
        incar_r2scan["GGA"] = None  # METAGGA supersedes GGA; leave unset
        write_incar(out_dir / "INCAR.r2scan", incar_r2scan)

    # LOBSTER basis: derived from the STATIC incar (the wavefunction it
    # will actually analyze), same pattern as prepare_vasp_lobster.py.
    potcar_variants = [_resolve_potcar_variant(el) for el in symbols_in_order]
    basis_lines = Lobsterin.get_basis(structure=structure, potcar_symbols=potcar_variants)
    dict_for_basis = dict(line.split(" ", 1) for line in basis_lines)
    lobsterin = Lobsterin.standard_calculations_from_vasp_files(
        POSCAR_input=str(out_dir / "POSCAR"),
        INCAR_input=str(out_dir / "INCAR.static"),
        dict_for_basis=dict_for_basis,
        option="standard",
    )
    lobsterin.write_lobsterin(str(out_dir / "lobsterin"))
    lobsterin.write_INCAR(
        incar_input=str(out_dir / "INCAR.static"),
        incar_output=str(out_dir / "INCAR.static"),
        poscar_input=str(out_dir / "POSCAR"),
        isym=-1,
    )

    submit_text = SLURM_TEMPLATE_STANDARD.format(job_name=f"mval_{name}")
    if spec["workflow"] == "mn":
        submit_text += SLURM_TEMPLATE_MN_EXTRA
    (out_dir / "submit.sh").write_text(submit_text)

    # Provenance: which existing project structure this started from, and
    # the manuscript's own optimized values for direct comparison later.
    (out_dir / "SOURCE.txt").write_text(
        f"Starting structure: {spec['source']} (mp_dataset/structures/{spec['source']}/POSCAR)\n"
        f"Workflow: {spec['workflow']}\nk-points (Table S1): {spec['kpoints']}\n"
    )

    print(f"prepared manuscript_{name}  (from {spec['source']}, workflow={spec['workflow']}, "
          f"kpoints={spec['kpoints']})")


def main():
    for name, spec in PHASES.items():
        prepare_one(name, spec)


if __name__ == "__main__":
    main()
