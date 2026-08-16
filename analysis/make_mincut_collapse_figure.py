"""Diagnostic figure for REPORT_dimensionality_mincut.md §5 (2026-08-16):
scatter of mincut_icohp_min_normalized vs. formation_energy_per_atom,
colored by which population a compound belongs to (original 186 /
elemental references added for reaction_icohp.py / extension4's
experimental half / extension4's deliberately far-above-hull half) --
makes the four different per-population trends (§5.1-5.3) visible in one
plot, instead of the single-color 186-only scatter this section
superseded.

Not part of the automated build_dataset_v2.py pipeline -- run manually,
same convention as report/gen_appendix.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
OUT = REPO_ROOT / "analysis" / "figures_v2" / "mincut_vs_formation_energy_by_population.png"

MAIN_FAMILIES = {"exp_stable", "exp_metastable", "theo_metastable", "hull", "metastable"}


def classify_population(row, ext4_kind: dict) -> str:
    if row["compound_id"] in ext4_kind:
        kind = ext4_kind[row["compound_id"]]
        return "extension4, experimental" if kind == "exp_polymorph" else "extension4, far-above-hull theoretical"
    if row["family"] in MAIN_FAMILIES:
        return "original 186"
    if abs(row["formation_energy_per_atom"]) < 0.02:
        return "elemental reference (FE~0)"
    return "other extension (batches 1-3)"


def main():
    df = pd.read_csv(REPO_ROOT / "analysis" / "percolation_vs_formation_energy.csv")
    picks = json.loads((REPO_ROOT / "mp_dataset" / "extension4_candidates.json").read_text())
    kind_suffix = {"exp_polymorph": "exp", "theo_far_hull": "theo"}
    ext4_kind = {
        f"extension_{p['formula']}_{kind_suffix[p['kind']]}_{p['material_id']}": p["kind"] for p in picks
    }

    df["population"] = df.apply(classify_population, axis=1, ext4_kind=ext4_kind)

    colors = {
        "original 186": "#3b6fa0",
        "elemental reference (FE~0)": "#999999",
        "other extension (batches 1-3)": "#c9b458",
        "extension4, experimental": "#5a9b5a",
        "extension4, far-above-hull theoretical": "#c0764a",
    }
    order = list(colors)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for pop in order:
        sub = df[df["population"] == pop].dropna(subset=["mincut_icohp_min_normalized", "formation_energy_per_atom"])
        if sub.empty:
            continue
        ax.scatter(
            sub["formation_energy_per_atom"], sub["mincut_icohp_min_normalized"],
            label=f"{pop} (n={len(sub)})", color=colors[pop], alpha=0.75, edgecolors="none", s=22,
        )
    ax.set_xlabel("Formation energy (eV/atom)")
    ax.set_ylabel("Min-cut, normalized")
    ax.set_title("Periodic min-cut vs. formation energy, by population\n(REPORT_dimensionality_mincut.md §5)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=150)
    plt.close(fig)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
