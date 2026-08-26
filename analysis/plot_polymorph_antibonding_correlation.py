"""Figure for analysis/compute_polymorph_antibonding_correlation.py:
raw ICOHP antibonding population vs. enthalpy above the same-formula
ground state (dH_eV_at), one panel per icobi_label bond-type group.

Writes analysis/figures_antibonding/polymorph_antibonding_correlation.png.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent
CSV = HERE / "polymorph_antibonding_correlation.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

COLORS = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a", "mixed": "#a05a9b"}
PANELS = [("ionic", "Ionic"), ("covalent", "Covalent"), ("mixed", "Mixed"), ("metallic", "Metallic")]


def main() -> None:
    df = pd.read_csv(CSV)

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    for ax, (bond, title) in zip(axes.flat, PANELS):
        sub = df[df["bond_type_group"] == bond]
        ax.scatter(
            sub["antibond_w_raw"], sub["dH_eV_at"],
            color=COLORS[bond], alpha=0.75, edgecolors="none", s=26,
        )
        ax.set_title(f"{title} (n={len(sub)} rows, {sub['formula'].nunique()} formulas)", fontsize=10)
        ax.set_xlabel(r"ICOHP$_\mathrm{antibonding}$ (eV)")
        ax.set_ylabel(r"$\Delta H$ vs. ground state (eV/atom)")

    fig.suptitle("Polymorph enthalpy vs. ICOHP antibonding population, by bond type\n"
                 "(no significant correlation in any group, see stats_summary_polymorph_antibonding.json)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = FIG_DIR / "polymorph_antibonding_correlation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
