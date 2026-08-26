"""Formation energy (decomposition-enthalpy proxy) vs. reaction
Delta(ICOHP_antibonding), faceted by bond_type (covalent / ionic /
metallic / others), with the viable/non-viable frontier from
classify_viability() (reaction_analysis/classify.py) overlaid.

Data: analysis/case1_viability.csv (viability_label, ground truth, from
mission #5's total delta_icohp + delta_energy=-formation_energy_per_atom)
merged with analysis/delta_antibonding_case1.csv (bond_type,
delta_icohp_antibond, mission #4b).

Two reference lines per panel:
  - y=0 (formation_energy_per_atom=0): EXACT half of the real boundary,
    descriptor-independent -- classify_viability() returns STABLE_ON_HULL
    (viable) for every compound at or below this line regardless of any
    ICOHP-derived quantity (delta_energy=-formation_energy_per_atom>=0).
  - x=0 (delta_icohp_antibond=0): the sign boundary tested by
    analysis/compute_viability_antibonding_correlation.py, shown as a
    visual reference for whether the antibonding descriptor's own sign
    tracks the real labels -- NOT the literal classify_viability()
    boundary, which uses mission #5's total delta_icohp
    (delta_icohp_per_atom_eV, plotted nowhere here) rather than the
    antibonding-only variant. The distinction is called out in the
    caption/legend rather than overstated.

Writes analysis/figures_antibonding/antibonding_viability_by_bondtype.png.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
VIABILITY_CSV = HERE / "case1_viability.csv"
ANTIBOND_CSV = HERE / "delta_antibonding_case1.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)

VIABLE_LABELS = {"stable_on_hull", "metastable_viable"}
NONVIABLE_LABEL = "unstable_nonexistent"

VIABLE_COLOR = "#2a9d5c"
NONVIABLE_COLOR = "#c0392b"
INSUFFICIENT_COLOR = "#999999"

PANELS = [
    ("covalent", "Covalent"),
    ("ionic", "Ionic"),
    ("metallic", "Metallic"),
    ("others", "Others (mixed / unclassified)"),
]


def main() -> None:
    via = pd.read_csv(VIABILITY_CSV)
    ab = pd.read_csv(ANTIBOND_CSV)
    df = via.merge(
        ab[["compound_id", "bond_type", "delta_icohp_antibond"]],
        on="compound_id", how="inner",
    )
    df["bond_type"] = df["bond_type"].fillna("others")
    df.loc[~df["bond_type"].isin(["covalent", "ionic", "metallic"]), "bond_type"] = "others"

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=False, sharey=True)
    for ax, (bt, title) in zip(axes.flat, PANELS):
        sub = df[df["bond_type"] == bt].dropna(subset=["delta_icohp_antibond", "formation_energy_per_atom_eV"])

        for label, color, marker_name in [
            (VIABLE_LABELS, VIABLE_COLOR, "viable"),
            ({NONVIABLE_LABEL}, NONVIABLE_COLOR, "non-viable"),
            ({"insufficient_data"}, INSUFFICIENT_COLOR, "insufficient data"),
        ]:
            g = sub[sub["viability_label"].isin(label)]
            if g.empty:
                continue
            ax.scatter(
                g["delta_icohp_antibond"], g["formation_energy_per_atom_eV"],
                label=f"{marker_name} (n={len(g)})", color=color, alpha=0.75,
                edgecolors="none", s=26,
            )

        ax.axhline(0, color="black", linewidth=0.8, linestyle="-")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.text(
            0.03, 0.96, "non-viable only\npossible here",
            transform=ax.transAxes, fontsize=6.5, color=NONVIABLE_COLOR,
            va="top", ha="left", style="italic",
        )

        # Clip a handful of extreme-magnitude compounds (e.g. one ionic
        # outlier at -0.71 eV) out of the visible range so the rest of the
        # panel's population stays readable; Spearman rho/p is still
        # computed on the unclipped data below.
        if len(sub) >= 10:
            lo, hi = sub["delta_icohp_antibond"].quantile([0.02, 0.98])
            pad = 0.2 * (hi - lo) if hi > lo else 0.01
            xlim = (min(lo - pad, -0.001), hi + pad)
            n_clipped = int((sub["delta_icohp_antibond"] < xlim[0]).sum()) + int((sub["delta_icohp_antibond"] > xlim[1]).sum())
            ax.set_xlim(*xlim)
            if n_clipped:
                ax.text(
                    0.97, 0.04, f"{n_clipped} pt(s) off-scale",
                    transform=ax.transAxes, fontsize=6.5, color="gray", ha="right",
                )

        rho_txt = ""
        vc = sub.dropna(subset=["delta_icohp_antibond", "formation_energy_per_atom_eV"])
        if len(vc) >= 5:
            rho, p = spearmanr(vc["delta_icohp_antibond"], vc["formation_energy_per_atom_eV"])
            rho_txt = f"  ($\\rho$={rho:.3f}, $n$={len(vc)})"
        ax.set_title(title + rho_txt, fontsize=10)
        ax.set_xlabel(r"$\Delta$(ICOHP$_\mathrm{antibonding}$) (eV)")

    axes[0, 0].set_ylabel("Formation energy (eV/atom)")
    axes[1, 0].set_ylabel("Formation energy (eV/atom)")
    axes[0, 1].legend(loc="upper right", fontsize=8, framealpha=0.9)

    fig.suptitle(
        r"$\Delta$(ICOHP$_\mathrm{antibonding}$) vs. formation energy, by bond type"
        "\nsolid line = exact thermodynamic boundary (formation energy = 0); "
        "dashed line = descriptor sign reference",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = FIG_DIR / "antibonding_viability_by_bondtype.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
