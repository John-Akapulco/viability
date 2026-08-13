"""Snapshot the overnight campaign's progress: for every compound in
job_list.txt (+ the 2 already-submitted smoke-test compounds), classify as
done / running / pending / failed, and for failures, look at the SLURM
array log + vasp.log + OSZICAR + lobsterout to name a probable cause
instead of just "failed". Updates mp_dataset/run_status.json in place.

Usage: python mp_dataset/check_campaign_status.py
"""

import calendar
import json
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).parent
JOB_LIST = ROOT / "job_list.txt"
STATUS_PATH = ROOT / "run_status.json"
SLURM_LOGS = ROOT / "slurm_logs"


_BRACKET_RE = re.compile(r"^(\d+)_\[(\d+)-(\d+)(?:%\d+)?\]$")
_SINGLE_RE = re.compile(r"^(\d+)_(\d+)$")


def squeue_snapshot():
    """Return (running_task_ids, pending_task_ids) as sets of '<jobid>_<idx>'
    strings, expanding SLURM's bracket range notation for still-queued array
    tasks (e.g. '57469_[9-177%8]' -> 57469_9 .. 57469_177)."""
    out = subprocess.run(
        ["squeue", "-u", "gilles", "-h", "-o", "%i %T"],
        capture_output=True, text=True, check=True,
    ).stdout
    running: set[str] = set()
    pending: set[str] = set()
    for line in out.splitlines():
        jobid, state = line.rsplit(" ", 1)
        target = running if state == "RUNNING" else pending if state == "PENDING" else None
        if target is None:
            continue
        m = _BRACKET_RE.match(jobid)
        if m:
            base, lo, hi = m.group(1), int(m.group(2)), int(m.group(3))
            target.update(f"{base}_{i}" for i in range(lo, hi + 1))
            continue
        m = _SINGLE_RE.match(jobid)
        if m:
            target.add(jobid)
        else:
            target.add(jobid)  # plain (non-array) job id, e.g. the 57468 smoke test
    return running, pending


def diagnose_failure(compound_dir: Path, array_task_id: int | None) -> str:
    reasons = []

    if array_task_id is not None:
        slurm_out = SLURM_LOGS / f"57469_{array_task_id}.out"
        if slurm_out.exists():
            text = slurm_out.read_text(errors="replace")
            if re.search(r"[Ss]egmentation fault", text):
                reasons.append("VASP segfault (SIGSEGV)")
            if "oom-kill" in text.lower() or "out of memory" in text.lower():
                reasons.append("out of memory")
            if "DUE TO TIME LIMIT" in text:
                reasons.append("SLURM walltime exceeded")
            if "CANCELLED" in text:
                reasons.append("job cancelled")

    vasp_log = compound_dir / "vasp.log"
    if vasp_log.exists():
        text = vasp_log.read_text(errors="replace")
        if "ZBRENT" in text:
            reasons.append("VASP ZBRENT (ionic step) error")
        if "EDDDAV" in text and "ZHEGV" not in text:
            reasons.append("VASP EDDDAV: sub-space matrix not hermitian (basis/NBANDS issue)")
        if "internal error in RSPHER" in text:
            reasons.append("VASP RSPHER internal error")
        if "General timing and accounting" not in text and not reasons:
            reasons.append("VASP did not reach normal termination (log truncated)")

    oszicar = compound_dir / "OSZICAR"
    if oszicar.exists() and not reasons:
        lines = oszicar.read_text(errors="replace").strip().splitlines()
        if not lines:
            reasons.append("OSZICAR empty (VASP did not start SCF)")

    lobsterout = compound_dir / "lobsterout"
    if (vasp_log.exists() and "General timing and accounting" in vasp_log.read_text(errors="replace")
            and not lobsterout.exists()):
        reasons.append("VASP finished but LOBSTER never ran / produced no lobsterout")
    elif lobsterout.exists():
        text = lobsterout.read_text(errors="replace")
        if "ERROR" in text.upper():
            err_lines = [ln for ln in text.splitlines() if "error" in ln.lower()]
            reasons.append(f"LOBSTER error: {err_lines[0][:200]}" if err_lines else "LOBSTER reported an error")
        elif "finished" not in text:
            reasons.append("LOBSTER did not report normal completion")

    return "; ".join(reasons) if reasons else "unknown (no diagnostic signature matched)"


def main():
    compound_dirs = [Path(line.strip()) for line in JOB_LIST.read_text().splitlines() if line.strip()]
    compound_dirs.append(Path("mp_dataset/structures/exp_stable_Ta5Ge3_mp-17593"))

    running_ids, pending_ids = squeue_snapshot()

    n_done = n_running = n_pending = n_failed = 0
    failures = []
    for idx, d in enumerate(compound_dirs):
        full = ROOT.parent / d if not d.is_absolute() else d
        if (full / "ICOHPLIST.lobster").exists():
            n_done += 1
            continue
        # crude: array task ids correspond to job_list.txt line order (0-indexed);
        # the appended Ta5Ge3 entry has no array id (submitted individually as 57468)
        array_task_id = idx if idx < len(compound_dirs) - 1 else None
        job_tag = f"57469_{array_task_id}" if array_task_id is not None else "57468"
        if job_tag in running_ids:
            n_running += 1
        elif job_tag in pending_ids:
            n_pending += 1
        elif not (full / "vasp.log").exists() and not (full / "OSZICAR").exists():
            # no SLURM record and no output yet -> hasn't been picked up (or
            # this snapshot raced the scheduler); don't cry failure yet.
            n_pending += 1
        else:
            n_failed += 1
            reason = diagnose_failure(full, array_task_id)
            failures.append({"compound": d.name, "reason": reason})

    status = json.loads(STATUS_PATH.read_text()) if STATUS_PATH.exists() else {}

    elapsed_h = None
    if "campaign_start_utc" in status:
        start = time.strptime(status["campaign_start_utc"], "%Y-%m-%dT%H:%M:%SZ")
        # calendar.timegm (not time.mktime, which assumes localtime) for a
        # correct UTC-struct_time -> UTC-epoch conversion.
        elapsed_h = (calendar.timegm(time.gmtime()) - calendar.timegm(start)) / 3600

    status.update({
        "last_checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_hours": elapsed_h,
        "n_done": n_done,
        "n_running": n_running,
        "n_pending": n_pending,
        "n_failed_or_unaccounted": n_failed,
        "failures": failures,
    })
    STATUS_PATH.write_text(json.dumps(status, indent=2))

    print(f"done={n_done} running={n_running} pending={n_pending} failed={n_failed}"
          f"  elapsed={elapsed_h:.2f}h" if elapsed_h is not None else "")
    for f in failures:
        print(f"  FAILED {f['compound']}: {f['reason']}")


if __name__ == "__main__":
    main()
