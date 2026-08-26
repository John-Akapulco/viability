"""Scatter figure for this project's headline correlation: reaction
Delta(antibonding) (mission #4b, analysis/delta_antibonding_case1.csv)
against formation_energy_per_atom -- the strongest single correlation in
the project (ICOHP-based: rho=-0.6225, n=438) and the manuscript's
Table~tab:delta-antibonding result, which so far has no accompanying
figure (manuscript/manuscript.tex Sec.~res-delta-antibonding).

formation_energy_per_atom is used as the decomposition-enthalpy proxy on
the y-axis (per atom formation enthalpy of the compound from its
elements is the same quantity, opposite sign convention, as the per-atom
decomposition enthalpy the case-1 reaction targets) -- the descriptor
itself sits on the x-axis, per the requested layout.

Two panels, ICOHP- and ICOBI-based Delta(antibonding) side by side,
colored by bond_type (same palette as the rest of this project's
antibonding figures, analysis/stats_analysis_antibonding.py). Writes
analysis/figures_antibonding/delta_antibonding_vs_formation_energy.png.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
CSV = HERE / "delta_antibonding_case1.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

COLORS = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a", "mixed": "#a05a9b"}

PANELS = [
    ("delta_icohp_antibond", r"$\Delta$(ICOHP$_\mathrm{antibonding}$) (eV)"),
    ("delta_icobi_antibond", r"$\Delta$(ICOBI$_\mathrm{antibonding}$)"),
]


def main() -> None:
    df = pd.read_csv(CSV)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, (col, xlabel) in zip(axes, PANELS):
        sub_all = df.dropna(subset=[col, "formation_energy_per_atom"])
        rho, p = spearmanr(sub_all[col], sub_all["formation_energy_per_atom"])

        for bond_type, sub in sub_all.groupby("bond_type", dropna=False):
            label = bond_type if isinstance(bond_type, str) else "unclassified"
            ax.scatter(
                sub[col],
                sub["formation_energy_per_atom"],
                label=f"{label} (n={len(sub)})",
                color=COLORS.get(bond_type, "gray"),
                alpha=0.75,
                edgecolors="none",
                s=22,
            )
        ax.axvline(0, color="black", linewidth=0.6, linestyle=":")
        ax.set_xlabel(xlabel)
        ax.set_title(rf"$\rho$={rho:.3f}, $p$={p:.1e}, $n$={len(sub_all)}", fontsize=10)

        # Correlation (rank-based) is computed on the full population above;
        # the x-range is clipped to the robust bulk of the data (1st-99th
        # percentile, padded) so a handful of extreme-magnitude compounds
        # do not compress the rest of the population into an unreadable
        # sliver -- clipped points are marked with arrows at the axis edge
        # rather than silently dropped.
        lo, hi = sub_all[col].quantile([0.01, 0.99])
        pad = 0.15 * (hi - lo)
        xlim = (lo - pad, hi + pad)
        n_clipped_lo = int((sub_all[col] < xlim[0]).sum())
        n_clipped_hi = int((sub_all[col] > xlim[1]).sum())
        ax.set_xlim(*xlim)
        if n_clipped_lo:
            ax.annotate(
                f"{n_clipped_lo} point(s)\noff-scale", xy=(0.02, 0.03), xycoords="axes fraction",
                fontsize=7, color="gray", ha="left",
            )
        if n_clipped_hi:
            ax.annotate(
                f"{n_clipped_hi} point(s)\noff-scale", xy=(0.98, 0.03), xycoords="axes fraction",
                fontsize=7, color="gray", ha="right",
            )

    axes[0].set_ylabel("Formation energy (eV/atom)")
    axes[1].legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.suptitle(r"Reaction $\Delta$(antibonding) vs. formation energy (case-1 decomposition)")
    fig.tight_layout()
    out = FIG_DIR / "delta_antibonding_vs_formation_energy.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
