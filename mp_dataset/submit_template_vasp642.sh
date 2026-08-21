#!/bin/sh
#SBATCH --job-name=REPLACE_ME
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=vasp.log
#SBATCH --time=24:00:00

# Template for a VASP 6.4.2 run in this project's usual two-stage
# (INCAR.relax -> INCAR.static -> LOBSTER) layout -- see e.g.
# mp_dataset/structures/manuscript_S4N2/submit.sh for the VASP 6.5.0
# version this mirrors. Adapted from
# /home/gilles/Al_C/structures/Al4C3_Pnma_0GPa/vasp6_relax.sh (a separate,
# unrelated project) at the user's request, 2026-08-21.
#
# There is no `module load vasp/6.4.2` on this cluster (only vasp/5.4.4
# and vasp/6.5.0 are registered) -- 6.4.2 exists only as a raw binary tree
# under /opt/ohpc/pub/software/vasp.6.4.2, loaded via PATH below, same as
# the Al_C script does.
#
# To use: copy into a compound directory that has INCAR.relax/INCAR.static
# already prepared, set --job-name, then `sbatch submit_vasp642.sh`.

module purge
module load intel
module load impi/2021.13
module load lobster/5.1.1
export OMP_NUM_THREADS=1
export PATH=/opt/ohpc/pub/software/vasp.6.4.2:$PATH
ulimit -s unlimited

cp INCAR.relax INCAR
{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; } 2>> vasp.log
if ! grep -q "General timing" OUTCAR 2>/dev/null; then
    echo "VASP relax step did not finish cleanly -- aborting" >> vasp.log
    exit 1
fi
cp CONTCAR POSCAR
cp CONTCAR CONTCAR.relaxed

cp INCAR.static INCAR
{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; } 2>> vasp.log
if ! grep -q "General timing" OUTCAR 2>/dev/null; then
    echo "VASP static step did not finish cleanly -- aborting before LOBSTER" >> vasp.log
    exit 1
fi

export OMP_NUM_THREADS=16
{ time -p lobster-5.1.1 ; } 2> lobster_time.txt
