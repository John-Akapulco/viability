"""Correlation analysis for the ICOBI antibonding-population-near-frontier
descriptor (analysis/icobi_antibonding_all.csv), against
formation_energy_per_atom -- same statistical convention as
stats_analysis_antibonding.py (Spearman rank correlation, no assumed
linearity, n<15 groups flagged explicitly, no SISSO), run side by side
with the ICOHP version so the two descriptors are directly comparable on
the same population.

Also analyzes analysis/delta_antibonding_case1.csv (the case-1 reaction
deltas for both descriptors, see compute_delta_antibonding_case1.py for
the extensive/intensive convention note).

Writes analysis/stats_summary_icobi_antibonding.json and PNG figures
under analysis/figures_antibonding/.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
ICOBI_CSV = HERE / "icobi_antibonding_all.csv"
DELTA_CSV = HERE / "delta_antibonding_case1.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

TARGET = "formation_energy_per_atom"
SMALL_N_THRESHOLD = 15

PRIMARY_METRICS = {
    "icobi_antibond_w_raw": "ICOBI antibonding population near frontier, dE=1.0 (raw)",
    "icobi_antibond_w_normalized": "ICOBI antibonding population near frontier, dE=1.0 (normalized)",
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


def make_figure(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a"}
    for bond_type, sub in df.groupby("bond_type", dropna=False):
        sub = sub.dropna(subset=["icobi_antibond_w_normalized", TARGET])
        if sub.empty:
            continue
        label = bond_type if isinstance(bond_type, str) else "unclassified"
        ax.scatter(
            sub[TARGET], sub["icobi_antibond_w_normalized"],
            label=f"{label} (n={len(sub)})", color=colors.get(bond_type, "gray"),
            alpha=0.75, edgecolors="none",
        )
    ax.set_xlabel("Formation energy (eV/atom)")
    ax.set_ylabel("ICOBI antibonding population near frontier (normalized, dE=1.0)")
    ax.set_title("ICOBI antibonding population vs. formation energy, by bond type")
    ax.legend()
    fig.tight_layout()
    out = FIG_DIR / "icobi_antibonding_vs_formation_energy_by_bondtype.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def delta_correlation_table(df: pd.DataFrame) -> list[dict]:
    rows = []
    for target in ("formation_energy_per_atom", "energy_above_hull_eV_per_atom"):
        if target not in df.columns:
            continue
        for metric in ("delta_icohp_antibond", "delta_icobi_antibond"):
            sub = df[[metric, target]].dropna()
            n = len(sub)
            if n < 4:
                rows.append({"metric": metric, "target": target, "n": n, "rho": None, "p_value": None})
                continue
            rho, p = spearmanr(sub[metric], sub[target])
            rows.append({
                "metric": metric, "target": target, "n": n,
                "rho": round(float(rho), 4), "p_value": round(float(p), 4),
                "note": "small sample (n<15)" if n < SMALL_N_THRESHOLD else None,
            })
    return rows


def main():
    df = pd.read_csv(ICOBI_CSV)

    summary = {
        "n_total_compounds": len(df),
        "n_by_bond_type": df["bond_type"].value_counts(dropna=False).to_dict(),
        "n_by_is_metal": {str(k): int(v) for k, v in df["is_metal"].value_counts(dropna=False).items()},
        "correlations_primary": correlation_table(df, PRIMARY_METRICS),
    }

    fig = make_figure(df)
    summary["figure"] = str(fig.relative_to(HERE))

    if DELTA_CSV.exists():
        delta_df = pd.read_csv(DELTA_CSV)
        summary["delta_case1_n"] = len(delta_df)
        summary["delta_case1_correlations"] = delta_correlation_table(delta_df)

    (HERE / "stats_summary_icobi_antibonding.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
