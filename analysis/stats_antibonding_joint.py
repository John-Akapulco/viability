"""Joint statistical analysis of the two antibonding-population-near-
frontier descriptors -- ICOHP (mission #4, cohp_extraction.
antibonding_population_near_frontier) and ICOBI (this session's
addition, cohp_extraction.icobi_antibonding_population_near_frontier)
-- run side by side on the SAME, full 394-compound population, plus
their case-1 reaction deltas (compute_delta_antibonding_case1.py).

Same statistical convention throughout this project: Spearman rank
correlation, no assumed linearity, n<15 groups flagged explicitly, no
SISSO.

Inputs:
  analysis/icohp_antibonding_full.csv   (394: percolation_vs_antibonding.csv
                                          + icohp_antibonding_maxhull.csv)
  analysis/icobi_antibonding_all.csv    (392)
  analysis/delta_antibonding_case1.csv  (320 case-1 reactions, both descriptors)

Writes analysis/stats_summary_antibonding_joint.json and a comparison
figure under analysis/figures_antibonding/.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
ICOHP_CSV = HERE / "icohp_antibonding_full.csv"
ICOBI_CSV = HERE / "icobi_antibonding_all.csv"
DELTA_CSV = HERE / "delta_antibonding_case1.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

TARGET = "formation_energy_per_atom"
SMALL_N_THRESHOLD = 15


def spearman_row(df: pd.DataFrame, metric: str, target: str, group_label: str) -> dict:
    sub = df[[metric, target]].dropna()
    n = len(sub)
    if n < 4:
        return {"group": group_label, "metric": metric, "target": target, "n": n, "rho": None, "p_value": None}
    rho, p = spearmanr(sub[metric], sub[target])
    return {
        "group": group_label, "metric": metric, "target": target, "n": n,
        "rho": round(float(rho), 4), "p_value": round(float(p), 4),
        "note": "small sample (n<15) -- do not over-interpret" if n < SMALL_N_THRESHOLD else None,
    }


def correlation_table(df: pd.DataFrame, metric: str, target: str = TARGET) -> list[dict]:
    rows = [spearman_row(df, metric, target, "all")]
    for bond_type, sub in df.groupby("bond_type"):
        rows.append(spearman_row(sub, metric, target, f"bond_type={bond_type}"))
    for is_metal, sub in df.groupby("is_metal"):
        rows.append(spearman_row(sub, metric, target, f"is_metal={is_metal}"))
    return rows


def make_comparison_figure(icohp_df: pd.DataFrame, icobi_df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    colors = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a"}

    for bond_type, sub in icohp_df.groupby("bond_type", dropna=False):
        sub = sub.dropna(subset=["antibond_w_normalized", TARGET])
        if sub.empty:
            continue
        label = bond_type if isinstance(bond_type, str) else "unclassified"
        axes[0].scatter(sub[TARGET], sub["antibond_w_normalized"], label=f"{label} (n={len(sub)})",
                         color=colors.get(bond_type, "gray"), alpha=0.7, edgecolors="none")
    axes[0].set_xlabel("Formation energy (eV/atom)")
    axes[0].set_ylabel("ICOHP antibonding population, normalized")
    axes[0].set_title("ICOHP antibonding (n=394)")
    axes[0].legend(fontsize=8)

    for bond_type, sub in icobi_df.groupby("bond_type", dropna=False):
        sub = sub.dropna(subset=["icobi_antibond_w_normalized", TARGET])
        if sub.empty:
            continue
        label = bond_type if isinstance(bond_type, str) else "unclassified"
        axes[1].scatter(sub[TARGET], sub["icobi_antibond_w_normalized"], label=f"{label} (n={len(sub)})",
                         color=colors.get(bond_type, "gray"), alpha=0.7, edgecolors="none")
    axes[1].set_xlabel("Formation energy (eV/atom)")
    axes[1].set_ylabel("ICOBI antibonding population, normalized")
    axes[1].set_title("ICOBI antibonding (n=392)")
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    out = FIG_DIR / "icohp_vs_icobi_antibonding_comparison.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    icohp_df = pd.read_csv(ICOHP_CSV)
    icobi_df = pd.read_csv(ICOBI_CSV)
    delta_df = pd.read_csv(DELTA_CSV)

    summary = {
        "icohp_antibonding": {
            "n_total": len(icohp_df),
            "correlations_raw": correlation_table(icohp_df, "antibond_w_raw"),
            "correlations_normalized": correlation_table(icohp_df, "antibond_w_normalized"),
        },
        "icobi_antibonding": {
            "n_total": len(icobi_df),
            "correlations_raw": correlation_table(icobi_df, "icobi_antibond_w_raw"),
            "correlations_normalized": correlation_table(icobi_df, "icobi_antibond_w_normalized"),
        },
        "delta_case1": {
            "n_total": len(delta_df),
            "delta_icohp_antibond_vs_formation_energy_stratified": correlation_table(
                delta_df, "delta_icohp_antibond", "formation_energy_per_atom"),
            "delta_icobi_antibond_vs_formation_energy_stratified": correlation_table(
                delta_df, "delta_icobi_antibond", "formation_energy_per_atom"),
            "delta_icohp_antibond_vs_formation_energy": spearman_row(
                delta_df, "delta_icohp_antibond", "formation_energy_per_atom", "all"),
            "delta_icobi_antibond_vs_formation_energy": spearman_row(
                delta_df, "delta_icobi_antibond", "formation_energy_per_atom", "all"),
            "delta_icohp_antibond_vs_hull": spearman_row(
                delta_df, "delta_icohp_antibond", "energy_above_hull_eV_per_atom", "all"),
            "delta_icobi_antibond_vs_hull": spearman_row(
                delta_df, "delta_icobi_antibond", "energy_above_hull_eV_per_atom", "all"),
            "delta_icohp_vs_delta_icobi_cross_correlation": (
                lambda sub: {
                    "n": len(sub),
                    "rho": round(float(spearmanr(sub["delta_icohp_antibond"], sub["delta_icobi_antibond"])[0]), 4),
                    "p_value": round(float(spearmanr(sub["delta_icohp_antibond"], sub["delta_icobi_antibond"])[1]), 4),
                }
            )(delta_df[["delta_icohp_antibond", "delta_icobi_antibond"]].dropna()),
        },
    }

    fig = make_comparison_figure(icohp_df, icobi_df)
    summary["comparison_figure"] = str(fig.relative_to(HERE))

    (HERE / "stats_summary_antibonding_joint.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
