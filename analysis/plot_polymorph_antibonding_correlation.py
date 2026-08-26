"""Figure for analysis/compute_polymorph_antibonding_correlation.py:
raw ICOHP- and ICOBI-based antibonding population vs. enthalpy above the
same-formula ground state (dH_eV_at), one column per icobi_label
bond-type group, one row per descriptor.

Writes analysis/figures_antibonding/polymorph_antibonding_correlation.png.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
CSV = HERE / "polymorph_antibonding_correlation.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

COLORS = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a", "mixed": "#a05a9b"}
BOND_PANELS = [("ionic", "Ionic"), ("covalent", "Covalent"), ("mixed", "Mixed"), ("metallic", "Metallic")]
DESCRIPTOR_ROWS = [
    ("antibond_w_raw", r"ICOHP$_\mathrm{antibonding}$ (eV)"),
    ("icobi_antibond_w_raw", r"ICOBI$_\mathrm{antibonding}$"),
]


def main() -> None:
    df = pd.read_csv(CSV)

    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    for row, (col, xlabel) in enumerate(DESCRIPTOR_ROWS):
        for ax, (bond, title) in zip(axes[row], BOND_PANELS):
            sub = df[df["bond_type_group"] == bond].dropna(subset=[col, "dH_eV_at"])
            ax.scatter(
                sub[col], sub["dH_eV_at"],
                color=COLORS[bond], alpha=0.75, edgecolors="none", s=26,
            )
            rho, p = spearmanr(sub[col], sub["dH_eV_at"]) if len(sub) >= 3 else (float("nan"), float("nan"))
            ax.set_title(f"{title} ($\\rho$={rho:.2f}, $p$={p:.2f}, n={len(sub)})", fontsize=9)
            ax.set_xlabel(xlabel)
            if ax is axes[row][0]:
                ax.set_ylabel(r"$\Delta H$ vs. ground state (eV/atom)")

    fig.suptitle("Polymorph enthalpy vs. antibonding population, by bond type\n"
                 "top: ICOHP-based; bottom: ICOBI-based -- no significant correlation in any group "
                 "(see stats_summary_polymorph_antibonding.json)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = FIG_DIR / "polymorph_antibonding_correlation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
