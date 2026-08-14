#!/bin/sh
#SBATCH --job-name=viability_campaign
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=mp_dataset/slurm_logs/%A_%a.out
#SBATCH --time=24:00:00

# One compound per array task; %8 in the sbatch --array flag on the command
# line caps concurrency at 8 running tasks (shared cluster courtesy).
module purge
module load intel
module load impi/2021.13
module load vasp/6.5.0
module load lobster/5.1.1
export OMP_NUM_THREADS=1
ulimit -s unlimited

COMPOUND_DIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" mp_dataset/job_list.txt)
if [ -z "$COMPOUND_DIR" ]; then
  echo "No directory for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi
cd "$COMPOUND_DIR" || exit 1

# NOTE: `{ time -p srun ...; } 2> vasp_time.txt` (the original approach)
# reliably left vasp_time.txt empty and put the timing report at the tail of
# vasp.log instead -- some interaction between srun's remote-task I/O
# forwarding and the shell's fd bookkeeping under SLURM that didn't
# reproduce in a plain interactive shell. Wall-clock via $SECONDS sidesteps
# the whole redirect question.
vasp_start=$SECONDS
srun --mpi=pmi2 --ntasks=$SLURM_NTASKS --cpus-per-task=1 --threads-per-core=1 vasp_std > vasp.log 2>&1
echo "real $((SECONDS - vasp_start))" > vasp_time.txt

export OMP_NUM_THREADS=16
lobster_start=$SECONDS
lobster-5.1.1
echo "real $((SECONDS - lobster_start))" > lobster_time.txt
