#!/bin/bash
# Overnight campaign orchestrator: polls SLURM until the 178-job array
# either finishes or the 6h submission budget is exhausted (in which case
# pending-not-yet-started tasks are cancelled, running ones are left to
# finish), then runs percolation_path.py + the analysis join/stats scripts.
# Meant to be launched with Bash run_in_background and waited on via its
# completion notification -- not polled manually.
set -u
cd /home/gilles/viability
source .venv/bin/activate

ARRAY_JOB_ID=57469
BUDGET_HOURS=6
POLL_SECONDS=300
GRACE_SECONDS=2700   # 45 min grace after cutoff for already-running tasks to finish
cutoff_cancelled=false
cutoff_time=""

echo "[orchestrator] starting at $(date -u +%FT%TZ)"

while true; do
  python mp_dataset/check_campaign_status.py > /tmp/campaign_status_snapshot.txt 2>&1
  cat /tmp/campaign_status_snapshot.txt

  elapsed_h=$(python3 -c "import json; print(json.load(open('mp_dataset/run_status.json')).get('elapsed_hours') or 0)")
  n_done=$(python3 -c "import json; print(json.load(open('mp_dataset/run_status.json'))['n_done'])")
  n_running=$(python3 -c "import json; print(json.load(open('mp_dataset/run_status.json'))['n_running'])")
  n_pending=$(python3 -c "import json; print(json.load(open('mp_dataset/run_status.json'))['n_pending'])")

  if (( $(echo "$elapsed_h >= $BUDGET_HOURS" | bc -l) )) && [ "$cutoff_cancelled" = false ]; then
    echo "[orchestrator] 6h budget reached (elapsed=${elapsed_h}h) -- cancelling pending (not-yet-started) tasks"
    scancel --state=PENDING "$ARRAY_JOB_ID" 2>&1
    cutoff_cancelled=true
    cutoff_time=$(date +%s)
    python3 -c "
import json
s = json.load(open('mp_dataset/run_status.json'))
s['status'] = 'budget_exceeded_pending_cancelled'
json.dump(s, open('mp_dataset/run_status.json', 'w'), indent=2)
"
  fi

  if [ "$n_running" -eq 0 ] && [ "$n_pending" -eq 0 ]; then
    echo "[orchestrator] nothing running or pending -- campaign finished at $(date -u +%FT%TZ)"
    break
  fi

  if [ "$cutoff_cancelled" = true ]; then
    now=$(date +%s)
    if (( now - cutoff_time > GRACE_SECONDS )); then
      echo "[orchestrator] grace period after budget cutoff elapsed with $n_running still running -- stopping wait"
      break
    fi
  fi

  sleep "$POLL_SECONDS"
done

python3 -c "
import json
s = json.load(open('mp_dataset/run_status.json'))
s['status'] = 'compute_phase_done'
json.dump(s, open('mp_dataset/run_status.json', 'w'), indent=2)
"

echo "[orchestrator] running percolation_path.py over mp_dataset/structures/"
python percolation_path.py --root mp_dataset/structures --metric icohp,icobi \
  --output mp_dataset/results.csv --also-json mp_dataset/results.json -v \
  > mp_dataset/percolation_path_run.log 2>&1
echo "[orchestrator] percolation_path.py exit=$?"

echo "[orchestrator] running analysis/build_dataset.py"
python analysis/build_dataset.py > analysis/build_dataset.log 2>&1
echo "[orchestrator] build_dataset.py exit=$?"

echo "[orchestrator] running analysis/stats_analysis.py"
python analysis/stats_analysis.py > analysis/stats_analysis.log 2>&1
echo "[orchestrator] stats_analysis.py exit=$?"

echo "ORCHESTRATOR_DONE $(date -u +%FT%TZ)"
