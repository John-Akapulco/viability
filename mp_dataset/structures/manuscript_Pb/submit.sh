#!/bin/sh
#SBATCH --job-name=mval_Pb
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
{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; } 2>> vasp.log
cp CONTCAR POSCAR
cp CONTCAR CONTCAR.relaxed

cp INCAR.static INCAR
{ time -p srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std ; } 2>> vasp.log

export OMP_NUM_THREADS=16
{ time -p lobster-5.1.1 ; } 2> lobster_time.txt
