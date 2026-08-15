"""Correlation analysis for the antibonding-population-near-frontier metric
on the extension_* batch (analysis/antibonding_extension.csv, written by
compute_antibonding_extension.py) -- the targeted diagnostic for known/
expected-unstable and deliberately off-equilibrium structures. Same
statistical convention as stats_analysis_antibonding.py: Spearman rank
correlation, no assumed linearity, n<15 groups flagged explicitly.

Target is energy_above_hull_eV_per_atom, NOT formation_energy_per_atom --
the extension metadata never fetched formation_energy_per_atom (it wasn't
needed for the original diagnostic purpose, and the 2 COD-sourced compounds
have no MP entry to fetch it from). Per project convention a result on one
target must not be assumed to transfer to the other -- this is a different
test from stats_analysis_antibonding.py, not a re-run of it.

Caveat baked into every summary: this is a hand-picked, non-random set of
compounds spanning very different chemistries (elemental metals, molecular
solids, oxides, nitrides, carbon allotropes...), not a systematic sample
like the 186-compound campaign. Treat any correlation here as weaker
evidence than the equivalent number on the main campaign, and read it
alongside the per-compound table, not instead of it.

Writes analysis/stats_summary_antibonding_extension.json.
"""

import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
CSV_PATH = HERE / "antibonding_extension.csv"

TARGET = "energy_above_hull_eV_per_atom"
SMALL_N_THRESHOLD = 15
ZERO_TOL = 1e-6

PRIMARY_METRICS = {
    "antibond_w_raw": "antibonding population near frontier, dE=1.0 (raw)",
    "antibond_w_normalized": "antibonding population near frontier, dE=1.0 (normalized)",
}

SENSITIVITY_METRICS = {
    "antibond_w_raw_dE0.5": "antibonding population, raw, dE=0.5",
    "antibond_w_raw_dE1.0": "antibonding population, raw, dE=1.0",
    "antibond_w_raw_dE2.0": "antibonding population, raw, dE=2.0",
}

# Known LOBSTER projection-quality issues, documented at compute time --
# flagged explicitly rather than silently trusted or silently dropped.
FLAGGED_COMPOUNDS = {
    "extension_CaN_mp-1058549": (
        "LOBSTER band-overlap maxDeviation up to ~9 (vs <=0.03 typical for "
        "the rest of the batch), 153/729 poorly-orthonormalized k-points -- "
        "projection-quality flag, not a crash. CaN was deliberately chosen "
        "as a metallic/unstable structure, which is a plausible root cause, "
        "but that doesn't itself guarantee the near-E_F ICOHP/ICOBI data "
        "this metric integrates is trustworthy here."
    ),
    "extension_CaO_mp-2605": (
        "Same caveat as CaN: LOBSTER band-overlap maxDeviation up to ~17, "
        "142/729 poorly-orthonormalized k-points. This CaO polymorph was "
        "deliberately chosen as an artificial, non-ground-state structure."
    ),
}


def spearman_row(df: pd.DataFrame, metric: str, group_label: str) -> dict:
    sub = df[[metric, TARGET]].dropna()
    n = len(sub)
    if n < 4:
        return {"group": group_label, "metric": metric, "n": n, "rho": None, "p_value": None,
                "note": "too few points for a correlation"}
    rho, p = spearmanr(sub[metric], sub[TARGET])
    return {
        "group": group_label,
        "metric": metric,
        "n": n,
        "rho": round(float(rho), 4),
        "p_value": round(float(p), 4),
        "note": "small sample (n<15) -- do not over-interpret" if n < SMALL_N_THRESHOLD else None,
    }


def correlation_table(df: pd.DataFrame, metrics: dict) -> list[dict]:
    rows = []
    for metric in metrics:
        rows.append(spearman_row(df, metric, "all"))
        for bond_type, sub in df.groupby("bond_type"):
            rows.append(spearman_row(sub, metric, f"bond_type={bond_type}"))
        for is_metal, sub in df.groupby("is_metal"):
            rows.append(spearman_row(sub, metric, f"is_metal={is_metal}"))
    return rows


def zero_floor_stats(df: pd.DataFrame) -> dict:
    raw = df["antibond_w_raw"]
    near_zero = raw.abs() < ZERO_TOL
    return {
        "n_total": len(df),
        "n_zero_raw": int(near_zero.sum()),
        "frac_zero_raw": round(float(near_zero.mean()), 4) if len(df) else None,
    }


def flagged_compound_check(df: pd.DataFrame) -> dict:
    out = {}
    for cid, reason in FLAGGED_COMPOUNDS.items():
        row = df[df["compound_id"] == cid]
        if row.empty:
            out[cid] = {"present_in_csv": False, "reason": reason}
        else:
            r = row.iloc[0]
            out[cid] = {
                "present_in_csv": True,
                "reason": reason,
                "antibond_w_normalized": r["antibond_w_normalized"],
                "energy_above_hull_eV_per_atom": r.get("energy_above_hull_eV_per_atom"),
            }
    return out


def main():
    df = pd.read_csv(CSV_PATH)
    n_with_target = int(df[TARGET].notna().sum())

    summary = {
        "n_total_rows": len(df),
        "n_with_target": n_with_target,
        "target": TARGET,
        "caveat_sample": (
            "Hand-picked, non-random compound set spanning very different "
            "chemistries -- not a systematic sample like the 186-compound "
            "campaign. Treat correlations here as weaker evidence than the "
            "equivalent test in stats_analysis_antibonding.py, and read "
            "alongside the per-compound table (antibonding_extension.csv), "
            "not instead of it."
        ),
        "caveat_target": (
            "Target is energy_above_hull_eV_per_atom, not "
            "formation_energy_per_atom (unavailable for this batch) -- do "
            "not compare rho/p directly against the main campaign's "
            "formation_energy_per_atom results."
        ),
        "n_by_bond_type": df["bond_type"].value_counts(dropna=False).to_dict(),
        "n_by_is_metal": {str(k): int(v) for k, v in df["is_metal"].value_counts(dropna=False).items()},
        "zero_floor_stats": zero_floor_stats(df),
        "correlations_primary": correlation_table(df, PRIMARY_METRICS),
        "correlations_sensitivity_dE_all_group_only": [
            spearman_row(df, metric, "all") for metric in SENSITIVITY_METRICS
        ],
        "flagged_compounds": flagged_compound_check(df),
    }

    (HERE / "stats_summary_antibonding_extension.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
