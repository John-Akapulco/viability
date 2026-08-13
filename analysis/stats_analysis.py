"""Simple, interpretable statistics on analysis/percolation_vs_hull.csv:
Spearman correlation (percolation weight, raw and normalized, vs
energy_above_hull) overall and per bond_type, the same for the classic
ICOHP aggregates for comparison, and a small reference logistic regression
(stable vs metastable) with cross-validated AUC. No SISSO, no symbolic
regression -- see the mission brief for why.

Writes analysis/stats_summary.json (all numbers, for REPORT.md to quote)
and PNG figures under analysis/figures/.
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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).parent
CSV_PATH = HERE / "percolation_vs_hull.csv"
FIG_DIR = HERE / "figures"
FIG_DIR.mkdir(exist_ok=True)

SMALL_N_THRESHOLD = 15

CANDIDATE_METRICS = {
    "icohp_percolation_weight_min": "percolation weight (raw)",
    "icohp_percolation_weight_min_normalized": "percolation weight (normalized by |ICOHP_min|)",
    "icohp_sum": "ICOHP sum",
    "icohp_mean": "ICOHP mean",
    "icohp_min": "ICOHP min (strongest bond)",
    "icohp_max": "ICOHP max (weakest bond)",
}


def spearman_row(df: pd.DataFrame, metric: str, group_label: str) -> dict:
    sub = df[[metric, "energy_above_hull_eV_at"]].dropna()
    n = len(sub)
    if n < 4:
        return {"group": group_label, "metric": metric, "n": n, "rho": None, "p_value": None,
                "note": "too few points for a correlation"}
    rho, p = spearmanr(sub[metric], sub["energy_above_hull_eV_at"])
    return {
        "group": group_label,
        "metric": metric,
        "n": n,
        "rho": round(float(rho), 4),
        "p_value": round(float(p), 4),
        "note": "small sample (n<15) -- do not over-interpret" if n < SMALL_N_THRESHOLD else None,
    }


def correlation_table(df: pd.DataFrame) -> list[dict]:
    rows = []
    for metric in CANDIDATE_METRICS:
        rows.append(spearman_row(df, metric, "all"))
        for bond_type, sub in df.groupby("bond_type"):
            rows.append(spearman_row(sub, metric, f"bond_type={bond_type}"))
    return rows


def make_scatter_figure(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a"}
    for bond_type, sub in df.groupby("bond_type"):
        sub = sub.dropna(subset=["icohp_percolation_weight_min"])
        if sub.empty:
            continue
        ax.scatter(
            sub["energy_above_hull_eV_at"],
            sub["icohp_percolation_weight_min"],
            label=f"{bond_type} (n={len(sub)})",
            color=colors.get(bond_type, "gray"),
            alpha=0.75,
            edgecolors="none",
        )
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_{\mathrm{above\ hull}}$ (eV/atom)")
    ax.set_ylabel("Percolation weight (eV, log scale)")
    ax.set_title("Percolation weight vs. distance to the convex hull")
    ax.legend()
    fig.tight_layout()
    out = FIG_DIR / "percolation_vs_ehull.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def logistic_reference_model(df: pd.DataFrame) -> dict:
    sub = df.dropna(subset=["icohp_percolation_weight_min_normalized", "icohp_mean", "bond_type", "family"]).copy()
    sub["y_stable"] = (sub["family"] == "exp_stable").astype(int)

    if sub["y_stable"].nunique() < 2 or len(sub) < 20:
        return {"skipped": True, "reason": f"insufficient data (n={len(sub)}, classes={sub['y_stable'].nunique()})"}

    dummies = pd.get_dummies(sub["bond_type"], prefix="bond", drop_first=True)
    X = pd.concat(
        [sub[["icohp_percolation_weight_min_normalized", "icohp_mean"]].reset_index(drop=True),
         dummies.reset_index(drop=True)],
        axis=1,
    )
    y = sub["y_stable"].to_numpy()

    X_scaled = StandardScaler().fit_transform(X)

    n_splits = min(5, int(sub["y_stable"].value_counts().min()))
    if n_splits < 2:
        return {"skipped": True, "reason": f"smallest class has <2 members, cannot cross-validate (n={len(sub)})"}

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        aucs = cross_val_score(LogisticRegression(max_iter=1000), X_scaled, y, cv=cv, scoring="roc_auc")

    # also fit on the full set once, for the ROC curve figure (illustrative, not held-out)
    from sklearn.metrics import RocCurveDisplay

    model = LogisticRegression(max_iter=1000).fit(X_scaled, y)
    fig, ax = plt.subplots(figsize=(5, 5))
    RocCurveDisplay.from_estimator(model, X_scaled, y, ax=ax, name="LogisticRegression (full-fit, illustrative)")
    ax.plot([0, 1], [0, 1], linestyle="--", color="gray", linewidth=1)
    ax.set_title("Reference ROC: stable (on-hull) vs metastable")
    fig.tight_layout()
    roc_path = FIG_DIR / "roc_stable_vs_metastable.png"
    fig.savefig(roc_path, dpi=150)
    plt.close(fig)

    return {
        "skipped": False,
        "n": len(sub),
        "n_stable": int(sub["y_stable"].sum()),
        "n_metastable": int((1 - sub["y_stable"]).sum()),
        "features": list(X.columns),
        "cv_folds": n_splits,
        "cv_auc_mean": round(float(np.mean(aucs)), 4),
        "cv_auc_std": round(float(np.std(aucs)), 4),
        "cv_auc_per_fold": [round(float(a), 4) for a in aucs],
        "roc_figure": str(roc_path.relative_to(HERE)),
    }


def main():
    df = pd.read_csv(CSV_PATH)

    summary = {
        "n_total_compounds": len(df),
        "n_by_family": df["family"].value_counts(dropna=False).to_dict(),
        "n_by_bond_type": df["bond_type"].value_counts(dropna=False).to_dict(),
        "n_disconnected_or_missing_percolation": int(df["icohp_percolation_weight_min"].isna().sum()),
        "correlations": correlation_table(df),
        "logistic_reference_model": logistic_reference_model(df),
    }

    scatter_path = make_scatter_figure(df)
    summary["scatter_figure"] = str(scatter_path.relative_to(HERE))

    (HERE / "stats_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
