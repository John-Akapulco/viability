"""Watches the 5 stage-1 gas-dimer relaxation jobs (prepare_gasref_dimers.py,
SLURM 58188-58192); as each finishes, copies its relaxed CONTCAR back to
POSCAR and hands the directory to
prepare_extension_vasp_lobster.prepare_one() unchanged (the project's
standard static+LOBSTER convention), then submits that job. Polls every
60s (these are tiny 2-atom jobs, expected to be fast); 2h budget.
"""

import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from prepare_extension_vasp_lobster import prepare_one  # noqa: E402

STRUCTURES_ROOT = Path(__file__).parent / "structures"
JOBIDS_FILE = Path("/tmp/claude-1002/-home-gilles-viability/1b601a19-d0c6-4653-a79a-e63b346a333b/scratchpad/gasref_stage1_jobids.txt")
STAGE2_JOBIDS_FILE = Path("/tmp/claude-1002/-home-gilles-viability/1b601a19-d0c6-4653-a79a-e63b346a333b/scratchpad/gasref_stage2_jobids.txt")

POLL_SECONDS = 60
BUDGET_SECONDS = 7200


def squeue_ids() -> set[str]:
    out = subprocess.run(["squeue", "-u", "gilles", "-h", "-o", "%i"], capture_output=True, text=True).stdout
    return set(out.split())


def relaxation_converged(compound_dir: Path) -> bool:
    outcar = compound_dir / "OUTCAR"
    contcar = compound_dir / "CONTCAR"
    if not outcar.exists() or not contcar.exists():
        return False
    text = outcar.read_text(errors="replace")
    return "reached required accuracy" in text or "General timing" in text


def main():
    jobs = {}
    for line in JOBIDS_FILE.read_text().splitlines():
        jid, name = line.split()
        jobs[name] = {"jid": jid, "stage2_done": False}

    stage2_log = []
    start = time.time()
    while time.time() - start < BUDGET_SECONDS:
        running = squeue_ids()
        remaining = 0
        for name, info in jobs.items():
            if info["stage2_done"]:
                continue
            remaining += 1
            compound_dir = STRUCTURES_ROOT / name
            if info["jid"] in running:
                continue  # stage 1 still queued/running
            if not relaxation_converged(compound_dir):
                print(f"[gasref-stage2] {name}: stage-1 job left queue but no converged OUTCAR/CONTCAR yet -- check manually")
                info["stage2_done"] = True  # don't retry forever unattended; flagged for manual review
                continue

            poscar_relaxed = (compound_dir / "CONTCAR").read_text()
            (compound_dir / "POSCAR").write_text(poscar_relaxed)

            prepare_one(compound_dir)  # regenerates POTCAR/INCAR/KPOINTS/lobsterin/submit.sh (static+LOBSTER convention)

            result = subprocess.run(["sbatch", "submit.sh"], cwd=compound_dir, capture_output=True, text=True)
            jid2 = result.stdout.strip().split()[-1]
            stage2_log.append(f"{jid2} {name}")
            print(f"[gasref-stage2] {name}: stage 1 converged, submitted stage 2 (static+LOBSTER) as job {jid2}")
            info["stage2_done"] = True

        if remaining == 0:
            print("[gasref-stage2] all 5 compounds handled")
            break
        time.sleep(POLL_SECONDS)

    STAGE2_JOBIDS_FILE.write_text("\n".join(stage2_log) + "\n")
    print(f"[gasref-stage2] wrote {STAGE2_JOBIDS_FILE}")


if __name__ == "__main__":
    main()
