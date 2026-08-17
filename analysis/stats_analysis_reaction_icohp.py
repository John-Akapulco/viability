"""Correlation analysis for the reaction-ICOHP metric (case 1: decomposition
into elements, analysis/reaction_icohp_case1.csv) against both stability
targets used in this project. Same statistical convention as
stats_analysis.py / stats_analysis_antibonding.py: Spearman rank
correlation, no assumed linearity, n<15 groups flagged explicitly, no SISSO.

Two targets tested separately, never assumed to transfer (project
convention, see METRIC_DEFINITION_reaction_icohp.md and
project-viability-methodology):
- energy_above_hull_eV_per_atom: available for all case-1 rows.
- formation_energy_per_atom: primarily from percolation_vs_formation_energy.csv
  (186-compound main campaign), with mp_dataset/formation_energies.json used
  as a fallback for extension-only compounds not in that file (Ca3N2,
  Mn2O7, TiO2 polymorphs, elemental references, ...) -- fetch_formation_energy.py
  was rerun 2026-08-15 to cover all structures dirs, closing the coverage
  gap this script's docstring used to describe. Two COD-sourced compounds
  (S4N2, S4N4) have no mp_id and remain excluded, not imputed.

bond_type / is_metal are mostly unpopulated in reaction_icohp_case1.csv's
own mp_metadata.json-sourced columns (184/192 and 183/192 NaN respectively
-- most of these compounds were never run through the main campaign's bond
classification step). Where the compound overlaps with the main 186
campaign, the better-populated bond_type/is_metal from
percolation_vs_formation_energy.csv is used instead of the sparse one
already in reaction_icohp_case1.csv.

Comparison against prior descriptors (percolation weight, min-cut,
antibonding population) is restricted to the overlap subset with the main
campaign, so it is an apples-to-apples comparison on the same compounds and
same target.

Explicitly authorized follow-up to METRIC_DEFINITION_reaction_icohp.md
section 6 ("no statistical test ... has been run") -- this is that test.

Writes analysis/stats_summary_reaction_icohp.json and PNG figures under
analysis/figures_reaction_icohp/.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
REACTION_CSV = HERE / "reaction_icohp_case1.csv"
MAIN_CSV = HERE / "percolation_vs_formation_energy.csv"
ANTIBOND_CSV = HERE / "percolation_vs_antibonding.csv"
FIG_DIR = HERE / "figures_reaction_icohp"
FIG_DIR.mkdir(exist_ok=True)

TARGETS = ["energy_above_hull_eV_per_atom", "formation_energy_per_atom"]
SMALL_N_THRESHOLD = 15
PRIMARY_METRICS = {"delta_icohp_per_atom": "reaction ICOHP, decomposition into elements (per atom)"}
SENSITIVITY_METRICS = {"delta_icohp_total": "reaction ICOHP, decomposition into elements (total, unnormalized)"}
COMPARISON_METRICS = {
    "delta_icohp_per_atom": "reaction ICOHP (this mission)",
    "icohp_percolation_weight_min_normalized": "percolation weight (normalized, mission #1)",
    "mincut_icohp_min_normalized": "periodic min-cut (normalized, mission #3 headline)",
    "antibond_w_normalized": "antibonding population near frontier (normalized, mission #4 headline)",
}


FORMATION_ENERGIES_JSON = HERE.parent / "mp_dataset" / "formation_energies.json"


def load_merged() -> pd.DataFrame:
    r = pd.read_csv(REACTION_CSV)
    # r's own bond_type (written by compute_reaction_icohp_case1.py) is
    # now sourced from the is_metal-first ICOBI classifier (commit
    # 60fe81a via icohp_icobi_bondtype.csv), which covers every compound
    # with LOBSTER data -- keep it under its own name so the merge below
    # doesn't discard it in favor of `main`'s bond_type, which is only
    # classify()'s composition-only heuristic (same flaw, restricted to
    # the 186-compound main campaign to boot).
    # Same reasoning for is_metal: r's own column (also written by
    # compute_reaction_icohp_case1.py) is now sourced from
    # icohp_icobi_bondtype.csv's is_metal (build_is_metal_map(), 0 NaN
    # across the dataset) -- keep it out of the merge's way so it isn't
    # discarded in favor of `antibond`'s is_metal, which only covers the
    # 186-compound main campaign.
    r = r.rename(columns={"bond_type": "bond_type_icobi", "is_metal": "is_metal_icobi"})
    main = pd.read_csv(MAIN_CSV)[[
        "mp_id", "bond_type", "formation_energy_per_atom",
        "icohp_percolation_weight_min_normalized", "mincut_icohp_min_normalized",
    ]]
    antibond = pd.read_csv(ANTIBOND_CSV)[["mp_id", "is_metal", "antibond_w_normalized"]]

    # pandas merge() treats NaN as a matching join key, so any compound
    # missing mp_id (COD-sourced S4N2/S4N4, a few main-campaign rows that
    # never got one populated, gasref_ZnSn_NiAs, ...) would otherwise
    # fan out against every other NaN-mp_id row on the other side --
    # silently inflating n_total_case1_rows and every correlation's n.
    # Give each NaN a placeholder unique across all three frames (a plain
    # per-frame counter would let e.g. r's row 0 collide with main's row 0)
    # so it can never match anything.
    _placeholder_counter = 0
    for frame in (r, main, antibond):
        na_mask = frame["mp_id"].isna()
        n_na = int(na_mask.sum())
        frame.loc[na_mask, "mp_id"] = [
            f"__no_mp_id_{_placeholder_counter + i}__" for i in range(n_na)
        ]
        _placeholder_counter += n_na

    df = r.merge(main, on="mp_id", how="left")
    df = df.merge(antibond, on="mp_id", how="left")
    df["in_main_campaign"] = df["mp_id"].isin(main["mp_id"]) & ~df["mp_id"].str.startswith("__no_mp_id_")

    # Prefer the ICOBI-based label (covers ~all compounds with LOBSTER
    # data, not just the 186-compound main campaign); fall back to
    # `main`'s classify()-heuristic bond_type only where the ICOBI
    # classifier had nothing (compound outside mp_dataset/structures/
    # entirely, or missing ICOHPLIST/ICOBILIST -- see
    # compute_icohp_icobi_bondtype.py).
    df["bond_type"] = df["bond_type_icobi"].fillna(df["bond_type"])
    df = df.drop(columns=["bond_type_icobi"])

    # Same preference for is_metal: ICOBI-classifier value first, main-
    # campaign-only `antibond` value only as fallback.
    df["is_metal"] = df["is_metal_icobi"].fillna(df["is_metal"])
    df = df.drop(columns=["is_metal_icobi"])

    # formation_energy_per_atom for compounds outside the 186-compound main
    # campaign (extension_* only, e.g. Ca3N2, Mn2O7, TiO2 polymorphs):
    # mp_dataset/fetch_formation_energy.py was rerun 2026-08-15 to cover all
    # structures dirs, not just the original 186 -- fill the gap left by the
    # main-campaign-only merge above rather than re-deriving it.
    fe_all = json.loads(FORMATION_ENERGIES_JSON.read_text())
    fallback = df["mp_id"].map(fe_all)
    df["formation_energy_per_atom"] = df["formation_energy_per_atom"].fillna(fallback)
    return df


def spearman_row(df: pd.DataFrame, metric: str, target: str, group_label: str) -> dict:
    sub = df[[metric, target]].dropna()
    n = len(sub)
    if n < 4:
        return {"group": group_label, "metric": metric, "target": target, "n": n,
                "rho": None, "p_value": None, "note": "too few points for a correlation"}
    rho, p = spearmanr(sub[metric], sub[target])
    return {
        "group": group_label,
        "metric": metric,
        "target": target,
        "n": n,
        "rho": round(float(rho), 4),
        "p_value": round(float(p), 4),
        "note": "small sample (n<15) -- do not over-interpret" if n < SMALL_N_THRESHOLD else None,
    }


def correlation_table(df: pd.DataFrame, metrics: dict, target: str) -> list[dict]:
    rows = []
    for metric in metrics:
        rows.append(spearman_row(df, metric, target, "all"))
        for bond_type, sub in df.groupby("bond_type"):
            rows.append(spearman_row(sub, metric, target, f"bond_type={bond_type}"))
        for is_metal, sub in df.groupby("is_metal"):
            rows.append(spearman_row(sub, metric, target, f"is_metal={is_metal}"))
    return rows


def make_bondtype_figure(df: pd.DataFrame, target: str) -> Path:
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"ionic": "#3b6fa0", "covalent": "#5a9b5a", "metallic": "#c0764a"}
    for bond_type, sub in df.groupby("bond_type", dropna=False):
        sub = sub.dropna(subset=["delta_icohp_per_atom", target])
        if sub.empty:
            continue
        label = bond_type if isinstance(bond_type, str) else "unclassified"
        ax.scatter(
            sub[target], sub["delta_icohp_per_atom"],
            label=f"{label} (n={len(sub)})",
            color=colors.get(bond_type, "gray"), alpha=0.75, edgecolors="none",
        )
    xlabel = "Energy above hull (eV/atom)" if target == "energy_above_hull_eV_per_atom" else "Formation energy (eV/atom)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\Delta$ICOHP per atom (decomposition into elements)")
    ax.set_title(f"Reaction ICOHP vs. {target}, by bond type")
    ax.legend()
    fig.tight_layout()
    out = FIG_DIR / f"reaction_icohp_vs_{target}_by_bondtype.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    df = load_merged()

    summary = {
        "n_total_case1_rows": len(df),
        "n_in_main_campaign": int(df["in_main_campaign"].sum()),
        "n_by_family": df["family"].value_counts(dropna=False).to_dict(),
        "n_by_bond_type": df["bond_type"].value_counts(dropna=False).to_dict(),
        "n_by_is_metal": {str(k): int(v) for k, v in df["is_metal"].value_counts(dropna=False).items()},
    }

    for target in TARGETS:
        summary[f"correlations_primary__{target}"] = correlation_table(df, PRIMARY_METRICS, target)
        summary[f"correlations_sensitivity__{target}"] = [
            spearman_row(df, metric, target, "all") for metric in SENSITIVITY_METRICS
        ]

    main_subset = df[df["in_main_campaign"]]
    for target in TARGETS:
        summary[f"comparison_vs_prior_descriptors__{target}"] = correlation_table(
            main_subset, COMPARISON_METRICS, target
        )

    figs = {}
    for target in TARGETS:
        fig_path = make_bondtype_figure(df, target)
        figs[target] = str(fig_path.relative_to(HERE))
    summary["figures"] = figs

    (HERE / "stats_summary_reaction_icohp.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
