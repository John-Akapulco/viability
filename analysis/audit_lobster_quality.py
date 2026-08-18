"""Systematic LOBSTER-quality audit across every compound directory, to
catch the Cu4Au failure mode (100% k-point orthonormalization failure,
bandOverlaps maxDeviation up to 1.10, physically implausible ICOHP
values) elsewhere in the dataset, not just where it happened to be
noticed by inspecting an outlier in delta_per_atom_eV.

Two independent signals, cross-checked against each other rather than
either one alone (a single soft signal -- e.g. maxDeviation slightly
above the project's existing 0.1 flag threshold for elemental
references -- is NOT grounds for exclusion; this project has always
treated that as a documented caveat, not a hard cut):

1. LOBSTER's own self-report: k-point orthonormalization failure
   fraction (lobsterout) and maxDeviation across all reported band
   overlaps (bandOverlaps.lobster).
2. Physically-implausible output: this compound's own delta_per_atom_eV
   (reaction_analysis_case1_full.csv) compared against the whole case-1
   population via a robust (MAD-based, outlier-resistant) z-score --
   the same kind of check that caught Cu4Au (-380.5 eV/atom against a
   population where the next-worst case is -14.3).

A compound is flagged SEVERE (candidate for exclusion) only if BOTH:
  - k-point failure fraction == 1.0 (every single k-point failed, not
    a handful) AND maxDeviation > 0.5 (five times the existing soft
    threshold), AND
  - its delta_per_atom_eV robust z-score exceeds 10 (i.e. is a
    genuine, order-of-magnitude outlier, not just "a bit large").

Writes analysis/lobster_quality_audit.json and prints a summary table.
Does not modify any structure directory -- exclusion is a separate,
explicit step per compound, same as Cu4Au's.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"

KPOINT_FAIL_RE = re.compile(r"WARNING:\s*(\d+)\s+of\s+(\d+)\s+k-points could not be orthonormalized")
MAXDEV_RE = re.compile(r"maxDeviation is:\s*([\d.]+)")
SPILLING_RE = re.compile(r"abs\. charge spilling:\s*([\d.]+)%")


def parse_lobster_quality(compound_dir: Path) -> dict:
    result = {"kpoint_fail_frac": None, "max_band_deviation": None, "charge_spilling_pct": None}
    lobsterout = compound_dir / "lobsterout"
    if lobsterout.exists():
        text = lobsterout.read_text(errors="replace")
        m = KPOINT_FAIL_RE.search(text)
        if m:
            failed, total = int(m.group(1)), int(m.group(2))
            result["kpoint_fail_frac"] = failed / total if total else None
        m2 = SPILLING_RE.search(text)
        if m2:
            result["charge_spilling_pct"] = float(m2.group(1))
    band_overlaps = compound_dir / "bandOverlaps.lobster"
    if band_overlaps.exists():
        devs = [float(x) for x in MAXDEV_RE.findall(band_overlaps.read_text(errors="replace"))]
        if devs:
            result["max_band_deviation"] = max(devs)
    return result


def main():
    rows = []
    for d in sorted(STRUCTURES_ROOT.iterdir()):
        if not d.is_dir() or d.name.startswith("manuscript_"):
            continue
        if not (d / "lobsterout").exists():
            continue
        q = parse_lobster_quality(d)
        q["compound_id"] = d.name
        rows.append(q)

    df = pd.DataFrame(rows)
    print(f"{len(df)} compounds with a lobsterout file scanned")
    print(f"  with kpoint_fail_frac reported: {df['kpoint_fail_frac'].notna().sum()}")
    print(f"  with max_band_deviation reported: {df['max_band_deviation'].notna().sum()}")

    # Merge against delta_per_atom_eV for the physically-implausible-output check
    ra_path = REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv"
    ra = pd.read_csv(ra_path)[["compound_id", "delta_per_atom_eV"]]
    df = df.merge(ra, on="compound_id", how="left")

    vals = df["delta_per_atom_eV"].dropna()
    median = vals.median()
    mad = (vals - median).abs().median()
    # robust z-score (MAD-based, 1.4826 = consistency constant for normal dist)
    df["robust_z"] = (df["delta_per_atom_eV"] - median).abs() / (1.4826 * mad) if mad > 0 else np.nan

    severe = df[
        (df["kpoint_fail_frac"] == 1.0)
        & (df["max_band_deviation"] > 0.5)
        & (df["robust_z"] > 10)
    ]

    soft_flag = df[
        ((df["kpoint_fail_frac"] > 0) | (df["max_band_deviation"] > 0.1))
        & ~df["compound_id"].isin(severe["compound_id"])
    ]

    print(f"\nSEVERE (100% k-point failure AND maxDeviation>0.5 AND |z|>10): {len(severe)}")
    if len(severe):
        print(severe[["compound_id", "kpoint_fail_frac", "max_band_deviation", "delta_per_atom_eV", "robust_z"]].to_string(index=False))

    print(f"\nSoft-flag only (some k-point failures or maxDeviation>0.1, but not severe): {len(soft_flag)}")
    print(f"  (existing project convention: documented caveat, not exclusion)")

    out = {
        "n_scanned": len(df),
        "severe_candidates": severe["compound_id"].tolist(),
        "soft_flag_count": len(soft_flag),
    }
    (REPO_ROOT / "analysis" / "lobster_quality_audit.json").write_text(json.dumps(out, indent=2, default=str))
    df.to_csv(REPO_ROOT / "analysis" / "lobster_quality_audit.csv", index=False)
    print(f"\nWrote analysis/lobster_quality_audit.json and .csv")


if __name__ == "__main__":
    main()
