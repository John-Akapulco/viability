"""Correlation between the case-1 ViabilityLabel (analysis/case1_viability.csv,
reaction_analysis.classify.classify_viability(), thermodynamics-based) and
the mission-4b reaction-level antibonding descriptor
(delta_icohp_antibond, analysis/delta_antibonding_case1.csv) -- not yet
tested anywhere else in the project. Mission #4's own correlations
(stats_antibonding_joint.py) only go against formation_energy_per_atom/
energy_above_hull, never against the viable/non-viable label itself;
mission #5's ViabilityLabel correlations (compute_case1_viability.py) only
go against the reaction-ICOHP sign (BondingLabel), never against the
antibonding descriptor.

Same statistical convention as the rest of the project: viable =
STABLE_ON_HULL + METASTABLE_VIABLE vs. UNSTABLE_NONEXISTENT
(insufficient_data excluded); Fisher exact test on sign(delta_icohp_antibond)
> 0 (binary) x viable (binary), and Mann-Whitney U on the continuous
delta_icohp_antibond itself, same split -- mirroring exactly the two tests
already run for the mission-5 sign/continuous split in
analysis/test_delta_icohp_viability.py. Mandatory bond_type/is_metal
stratification, n<15 per group always flagged explicitly, no SISSO.

Writes analysis/stats_summary_viability_antibonding.json and a boxplot
figure under analysis/figures_antibonding/.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

HERE = Path(__file__).parent
VIABILITY_CSV = HERE / "case1_viability.csv"
ANTIBOND_CSV = HERE / "delta_antibonding_case1.csv"
FIG_DIR = HERE / "figures_antibonding"
FIG_DIR.mkdir(exist_ok=True)
OUT_JSON = HERE / "stats_summary_viability_antibonding.json"

VIABLE_LABELS = {"stable_on_hull", "metastable_viable"}
NONVIABLE_LABEL = "unstable_nonexistent"


def _test_group(name: str, g: pd.DataFrame) -> dict:
    vv = g.loc[g["viable"], "delta_icohp_antibond"]
    nv = g.loc[~g["viable"], "delta_icohp_antibond"]
    n_viable, n_nonviable = len(vv), len(nv)
    entry = {"group": name, "n_viable": n_viable, "n_nonviable": n_nonviable}
    if min(n_viable, n_nonviable) < 5:
        entry["note"] = "skipped: fewer than 5 in one class"
        return entry
    u, p_mwu = mannwhitneyu(vv, nv, alternative="two-sided")
    entry["mannwhitney_u"] = float(u)
    entry["mannwhitney_p"] = float(p_mwu)
    entry["median_viable"] = float(vv.median())
    entry["median_nonviable"] = float(nv.median())
    table = pd.crosstab(g["delta_icohp_antibond"] > 0, g["viable"])
    if table.shape == (2, 2):
        odds, p_fisher = fisher_exact(table)
        entry["fisher_odds_ratio_sign_positive_vs_viable"] = float(odds)
        entry["fisher_p"] = float(p_fisher)
    entry["note"] = "n<15 in at least one class" if min(n_viable, n_nonviable) < 15 else None
    return entry


def main() -> None:
    via = pd.read_csv(VIABILITY_CSV)
    ab = pd.read_csv(ANTIBOND_CSV)
    m = via.merge(
        ab[["compound_id", "delta_icohp_antibond", "delta_icobi_antibond", "bond_type", "is_metal"]],
        on="compound_id", how="inner",
    )
    m = m[m["viability_label"].isin(VIABLE_LABELS | {NONVIABLE_LABEL})].copy()
    m["viable"] = m["viability_label"].isin(VIABLE_LABELS)
    sub = m.dropna(subset=["delta_icohp_antibond"]).copy()

    results = {"n_total": len(sub), "correlations": []}
    results["correlations"].append(_test_group("all", sub))
    for bt in sorted(sub["bond_type"].dropna().unique()):
        results["correlations"].append(_test_group(f"bond_type={bt}", sub[sub["bond_type"] == bt]))
    for im in (True, False):
        results["correlations"].append(_test_group(f"is_metal={im}", sub[sub["is_metal"] == im]))

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.boxplot(
        [sub.loc[sub["viable"], "delta_icohp_antibond"], sub.loc[~sub["viable"], "delta_icohp_antibond"]],
        tick_labels=["viable\n(stable or metastable)", "not viable\n(unstable/nonexistent)"],
    )
    ax.set_ylabel(r"$\Delta$(ICOHP antibonding), products$-$reactant (eV)")
    ax.set_title("Case-1 viability vs. reaction $\\Delta$(ICOHP antibonding)")
    fig_path = FIG_DIR / "viability_vs_delta_antibonding.png"
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    results["figure"] = str(fig_path.relative_to(HERE.parent))

    OUT_JSON.write_text(json.dumps(results, indent=2))
    for row in results["correlations"]:
        print(row)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
