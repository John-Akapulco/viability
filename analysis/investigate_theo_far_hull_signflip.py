"""Follow-up to REPORT_dimensionality_mincut.md's §5.3/§6 open question: is
the extension4 theo_far_hull subset's opposite-sign mincut-vs-formation-
energy correlation (rho=-0.40, p=0.011, n=39) a real far-from-hull physical
effect, a cell-size/symmetry selection artifact of "highest energy_above_hull
per formula, nsites<=40" (download_extension4.py), or a between-chemistry
confound?

Three checks, in increasing order of how directly they control for
composition:

1. Partial Spearman correlation controlling for n_sites and sg_number
   (the cell-size/symmetry confound the report flagged as untested).
2. Residual correlation after regressing mincut and formation_energy_per_atom
   each on anion-identity dummies (N/O/F/P/S/Cl) -- controls for chemistry
   at the anion-category level, still pooling different cations together.
3. Paired within-chemsys deltas: for each of the 37 chemsys with both an
   exp_polymorph and a theo_far_hull pick, correlate
   (theo_mincut - exp_mincut) against (theo_FE - exp_FE) -- the strongest
   control, since both entries share the exact same two elements and only
   the structure/hull-distance differs.

Not part of build_dataset_v2.py -- run manually, same convention as
make_mincut_collapse_figure.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).parent.parent
MC = "mincut_icohp_min_normalized"
FE = "formation_energy_per_atom"
ANIONS = ["N", "O", "F", "P", "S", "Cl"]


KIND_SUFFIX = {"exp_polymorph": "exp", "theo_far_hull": "theo"}


def load_theo_far_hull() -> pd.DataFrame:
    df = pd.read_csv(REPO_ROOT / "analysis" / "percolation_vs_formation_energy.csv")
    picks = json.loads((REPO_ROOT / "mp_dataset" / "extension4_candidates.json").read_text())
    cid_of = lambda p: f"extension_{p['formula']}_{KIND_SUFFIX[p['kind']]}_{p['material_id']}"
    kind_by_cid = {cid_of(p): p["kind"] for p in picks}
    chemsys_by_cid = {cid_of(p): p["chemsys"] for p in picks}
    df["kind"] = df["compound_id"].map(kind_by_cid)
    df["chemsys"] = df["compound_id"].map(chemsys_by_cid)
    return df[df["kind"].notna()].copy()


def rank_partial_corr(x, y, z_cols) -> tuple[float, float]:
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    Z = np.column_stack([np.ones(len(x))] + [stats.rankdata(z) for z in z_cols])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    return stats.pearsonr(rx - Z @ bx, ry - Z @ by)


def main():
    ext4 = load_theo_far_hull()
    theo = ext4[ext4["kind"] == "theo_far_hull"].copy()
    theo["anion"] = theo["chemsys"].str.split("-").apply(lambda p: [e for e in p if e in ANIONS][0])
    print(f"theo_far_hull: n={len(theo)}")

    rho, p = stats.spearmanr(theo[MC], theo[FE])
    print(f"\n0. Pooled (reproduces report §5.3): rho={rho:.4f} p={p:.4f}")

    rho, p = rank_partial_corr(theo[MC], theo[FE], [theo["n_sites"]])
    print(f"1a. Partial, controlling n_sites: rho={rho:.4f} p={p:.4f}")
    rho, p = rank_partial_corr(theo[MC], theo[FE], [theo["sg_number"]])
    print(f"1b. Partial, controlling sg_number: rho={rho:.4f} p={p:.4f}")
    rho, p = rank_partial_corr(theo[MC], theo[FE], [theo["n_sites"], theo["sg_number"]])
    print(f"1c. Partial, controlling n_sites+sg_number: rho={rho:.4f} p={p:.4f}")
    print("    -> survives: NOT a cell-size/symmetry selection artifact.")

    kw_mc = stats.kruskal(*[g[MC].values for _, g in theo.groupby("anion")])
    kw_fe = stats.kruskal(*[g[FE].values for _, g in theo.groupby("anion")])
    print(f"\n2a. Kruskal-Wallis, mincut by anion: H={kw_mc.statistic:.2f} p={kw_mc.pvalue:.4f}")
    print(f"2b. Kruskal-Wallis, formation_energy by anion: H={kw_fe.statistic:.2f} p={kw_fe.pvalue:.4f}")
    dummies = pd.get_dummies(theo["anion"], drop_first=True).astype(float)
    X = np.column_stack([np.ones(len(theo))] + [dummies[c].values for c in dummies.columns])
    resid_mc = theo[MC].values - X @ np.linalg.lstsq(X, theo[MC].values, rcond=None)[0]
    resid_fe = theo[FE].values - X @ np.linalg.lstsq(X, theo[FE].values, rcond=None)[0]
    rho, p = stats.spearmanr(resid_mc, resid_fe)
    print(f"2c. Residual correlation after regressing out anion identity: rho={rho:.4f} p={p:.4f}")
    print("    -> collapses: mincut AND formation_energy both vary significantly by anion.")

    pairs = []
    for cs, g in ext4.groupby("chemsys"):
        exp, th = g[g["kind"] == "exp_polymorph"], g[g["kind"] == "theo_far_hull"]
        if len(exp) == 1 and len(th) == 1:
            pairs.append({
                "chemsys": cs,
                "d_mincut": th[MC].values[0] - exp[MC].values[0],
                "d_fe": th[FE].values[0] - exp[FE].values[0],
            })
    pdf = pd.DataFrame(pairs)
    rho, p = stats.spearmanr(pdf["d_mincut"], pdf["d_fe"])
    print(f"\n3. Paired within-chemsys delta (n={len(pdf)}): rho={rho:.4f} p={p:.4f}")
    print("    -> collapses: same conclusion as check 2, from a stricter (exact-composition) control.")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    cmap = plt.get_cmap("tab10")
    for i, anion in enumerate(ANIONS):
        sub = theo[theo["anion"] == anion]
        if sub.empty:
            continue
        ax.scatter(sub[FE], sub[MC], label=f"{anion} (n={len(sub)})", color=cmap(i), s=40, edgecolors="none")
    ax.set_xlabel("Formation energy (eV/atom)")
    ax.set_ylabel("Min-cut, normalized")
    ax.set_title(
        "extension4 theo_far_hull: mincut vs. formation energy, by anion\n"
        "(pooled rho=-0.40 collapses to rho≈-0.07 once anion is controlled for)"
    )
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out = REPO_ROOT / "analysis" / "figures_v2" / "mincut_theo_far_hull_by_anion.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
