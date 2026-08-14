"""Correlation analysis for the antibonding-population-near-frontier metric
(analysis/percolation_vs_antibonding.csv) against formation_energy_per_atom.
Same statistical convention as stats_analysis.py: Spearman rank correlation,
no assumed linearity, n<15 groups flagged explicitly, no SISSO.

Explicitly authorized follow-up to analysis/METRIC_DEFINITION_antibonding.md
step 2 -- tests whether the metric predicts anything at all, which that
document left open.

Writes analysis/stats_summary_antibonding.json and PNG figures under
analysis/figures_antibonding/.
"""

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
CSV_PATH = HERE / "percolation_vs_antibonding.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

TARGET = "formation_energy_per_atom"
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

COMPARISON_METRICS = {
    "icohp_percolation_weight_min": "percolation weight (raw)",
    "icohp_percolation_weight_min_normalized": "percolation weight (normalized)",
    "icohp_sum": "ICOHP sum",
    "icohp_mean": "ICOHP mean",
    "icohp_min": "ICOHP min (strongest bond)",
    "icohp_max": "ICOHP max (weakest bond)",
    "mincut_icohp_min_normalized": "periodic min-cut (normalized, mission #3 headline)",
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
    out = {
        "n_total": len(df),
        "n_zero_raw": int(near_zero.sum()),
        "frac_zero_raw": round(float(near_zero.mean()), 4),
    }
    for is_metal, sub in df.groupby("is_metal"):
        sub_zero = sub["antibond_w_raw"].abs() < ZERO_TOL
        out[f"is_metal={is_metal}"] = {
            "n": len(sub),
            "n_zero_raw": int(sub_zero.sum()),
            "frac_zero_raw": round(float(sub_zero.mean()), 4),
        }
    for bond_type, sub in df.groupby("bond_type"):
        sub_zero = sub["antibond_w_raw"].abs() < ZERO_TOL
        out[f"bond_type={bond_type}"] = {
            "n": len(sub),
            "n_zero_raw": int(sub_zero.sum()),
            "frac_zero_raw": round(float(sub_zero.mean()), 4),
        }
    return out


def make_bondtype_figure(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a"}
    for bond_type, sub in df.groupby("bond_type", dropna=False):
        sub = sub.dropna(subset=["antibond_w_normalized", TARGET])
        if sub.empty:
            continue
        label = bond_type if isinstance(bond_type, str) else "unclassified"
        ax.scatter(
            sub[TARGET],
            sub["antibond_w_normalized"],
            label=f"{label} (n={len(sub)})",
            color=colors.get(bond_type, "gray"),
            alpha=0.75,
            edgecolors="none",
        )
    ax.set_xlabel(r"Formation energy (eV/atom)")
    ax.set_ylabel("Antibonding population near frontier (normalized, dE=1.0)")
    ax.set_title("Antibonding population vs. formation energy, by bond type")
    ax.legend()
    fig.tight_layout()
    out = FIG_DIR / "antibonding_vs_formation_energy_by_bondtype.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def make_ismetal_figure(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {True: "#c0764a", False: "#3b6fa0"}
    labels = {True: "metal", False: "gapped"}
    for is_metal, sub in df.groupby("is_metal"):
        sub = sub.dropna(subset=["antibond_w_normalized", TARGET])
        if sub.empty:
            continue
        ax.scatter(
            sub[TARGET],
            sub["antibond_w_normalized"],
            label=f"{labels.get(is_metal, is_metal)} (n={len(sub)})",
            color=colors.get(is_metal, "gray"),
            alpha=0.75,
            edgecolors="none",
        )
    ax.set_xlabel(r"Formation energy (eV/atom)")
    ax.set_ylabel("Antibonding population near frontier (normalized, dE=1.0)")
    ax.set_title("Antibonding population vs. formation energy, metal vs. gapped")
    ax.legend()
    fig.tight_layout()
    out = FIG_DIR / "antibonding_vs_formation_energy_by_ismetal.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    df = pd.read_csv(CSV_PATH)

    summary = {
        "n_total_compounds": len(df),
        "n_by_bond_type": df["bond_type"].value_counts(dropna=False).to_dict(),
        "n_by_is_metal": {str(k): int(v) for k, v in df["is_metal"].value_counts(dropna=False).items()},
        "zero_floor_stats": zero_floor_stats(df),
        "correlations_primary": correlation_table(df, PRIMARY_METRICS),
        "correlations_sensitivity_dE_all_group_only": [
            spearman_row(df, metric, "all") for metric in SENSITIVITY_METRICS
        ],
        "correlations_comparison_existing_descriptors": correlation_table(df, COMPARISON_METRICS),
    }

    bondtype_fig = make_bondtype_figure(df)
    ismetal_fig = make_ismetal_figure(df)
    summary["bondtype_figure"] = str(bondtype_fig.relative_to(HERE))
    summary["ismetal_figure"] = str(ismetal_fig.relative_to(HERE))

    (HERE / "stats_summary_antibonding.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
