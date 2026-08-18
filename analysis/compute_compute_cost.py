"""VASP/LOBSTER wall-clock timing across every structure directory.

Supersedes build_dataset.py's parse_vasp_time_from_log()/parse_time_file(),
which only understand POSIX `time -p` output ("real 123.45"). The
individually-submitted batches (maxhull, marginal-formation-energy, widen
-- everything not run through the original submit_array.sh) use bash's
builtin `time` instead, which formats the same line as "real 1m54.739s".
The old regex ([\\d.]+) partially matches that string too (grabbing just
the leading "1"), silently producing a bogus ~1-second reading instead of
failing -- 424/597 directories are affected. Handles both formats
explicitly here instead.

LOBSTER timing similarly falls back from lobster_time.txt (present for
173/597, POSIX format, from the original submit_array.sh's `time -p`
wrapper) to lobsterout's own "finished in H h M min S s ... of wall time"
line (present in 596/597 lobsterout files, format-independent of how the
job was submitted) -- giving near-complete coverage either way.

Writes analysis/compute_cost.json.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT_JSON = Path(__file__).parent / "compute_cost.json"

POSIX_RE = re.compile(r"^real\s+([\d.]+)\s*$", re.MULTILINE)
BASH_RE = re.compile(r"^real\s+(?:(\d+)m)?([\d.]+)s\s*$", re.MULTILINE)
LOBSTEROUT_RE = re.compile(r"finished in\s+(\d+)\s*h\s+(\d+)\s*min\s+(\d+)\s*s\s+(\d+)\s*ms\s+of wall time")


def parse_real_time(text: str) -> float | None:
    m = POSIX_RE.search(text)
    if m:
        return float(m.group(1))
    m = BASH_RE.search(text)
    if m:
        minutes = float(m.group(1)) if m.group(1) else 0.0
        seconds = float(m.group(2))
        return minutes * 60 + seconds
    return None


def parse_vasp_time(compound_dir: Path) -> float | None:
    log = compound_dir / "vasp.log"
    if not log.exists():
        return None
    return parse_real_time(log.read_text(errors="replace"))


def parse_lobster_time(compound_dir: Path) -> float | None:
    time_file = compound_dir / "lobster_time.txt"
    if time_file.exists():
        t = parse_real_time(time_file.read_text(errors="replace"))
        if t is not None:
            return t
    lobsterout = compound_dir / "lobsterout"
    if lobsterout.exists():
        m = LOBSTEROUT_RE.search(lobsterout.read_text(errors="replace"))
        if m:
            h, mn, s, ms = (int(x) for x in m.groups())
            return h * 3600 + mn * 60 + s + ms / 1000
    return None


def summarize(times: list[float]) -> dict:
    if not times:
        return {"n": 0}
    s = sorted(times)
    n = len(s)
    median = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    return {
        "n": n,
        "mean_s": round(sum(s) / n, 1),
        "median_s": round(median, 1),
        "max_s": round(max(s), 1),
        "total_s": round(sum(s), 1),
        "total_h": round(sum(s) / 3600, 2),
    }


def main() -> None:
    dirs = sorted(d for d in STRUCTURES_ROOT.iterdir() if d.is_dir())
    vasp_times, lobster_times = [], []
    per_compound = []
    for d in dirs:
        vt = parse_vasp_time(d)
        lt = parse_lobster_time(d)
        if vt is not None:
            vasp_times.append(vt)
        if lt is not None:
            lobster_times.append(lt)
        per_compound.append({"compound_id": d.name, "vasp_wall_time_s": vt, "lobster_wall_time_s": lt})

    result = {
        "n_structure_dirs": len(dirs),
        "vasp": summarize(vasp_times),
        "lobster": summarize(lobster_times),
    }
    OUT_JSON.write_text(json.dumps({**result, "per_compound": per_compound}, indent=2))

    print(f"{len(dirs)} structure directories scanned")
    print(f"VASP:    n={result['vasp']['n']}  mean={result['vasp']['mean_s']}s  median={result['vasp']['median_s']}s  "
          f"max={result['vasp']['max_s']}s ({result['vasp']['max_s']/3600:.1f}h)  total={result['vasp']['total_h']}h")
    print(f"LOBSTER: n={result['lobster']['n']}  mean={result['lobster']['mean_s']}s  median={result['lobster']['median_s']}s  "
          f"max={result['lobster']['max_s']}s ({result['lobster']['max_s']/3600:.1f}h)  total={result['lobster']['total_h']}h")
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
