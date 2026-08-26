"""Central hypothesis test of this project, as reformulated 2026-08-16:
does the reaction-ICOHP descriptor (delta(ICOHP), decomposition into
elements) distinguish thermodynamically stable, metastable, and unstable
compounds?

Viability prediction test, extended to every element and compound
computed across the whole project to date (281 case-1 reactions,
deliberately not split by which historical campaign/extension batch a
compound came from -- family and theoretical are used as the ground
truth, not provenance): does delta_per_atom_eV (or its sign,
BondingLabel) discriminate compounds by real thermodynamic stability
(energy_above_hull, formation_energy_per_atom, experimental-vs-
theoretical-only, and the exp_stable/exp_metastable/theo_metastable
family split)?

Writes analysis/stats_summary_delta_icohp_viability.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, kruskal, mannwhitneyu, spearmanr

REPO_ROOT = Path(__file__).parent.parent
OUT_JSON = Path(__file__).parent / "stats_summary_delta_icohp_viability.json"


def load_full_history() -> pd.DataFrame:
    ri = pd.read_csv(REPO_ROOT / "analysis" / "reaction_icohp_case1.csv")[
        ["compound_id", "mp_id", "formula", "family", "theoretical", "energy_above_hull_eV_per_atom"]
    ]
    ra = pd.read_csv(REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv")[["compound_id", "delta_per_atom_eV"]]
    df = ri.merge(ra, on="compound_id", how="inner")
    fe = json.loads((REPO_ROOT / "mp_dataset" / "formation_energies.json").read_text())
    df["formation_energy_per_atom"] = df["mp_id"].map(fe)
    df["bonding_label"] = df["delta_per_atom_eV"].apply(lambda x: "endobondic" if x >= 0 else "exobondic")
    return df


def viability_test(df: pd.DataFrame) -> dict:
    out: dict = {"n_total": len(df)}

    for target in ("energy_above_hull_eV_per_atom", "formation_energy_per_atom"):
        sub = df[["delta_per_atom_eV", target]].dropna()
        rho, p = spearmanr(sub["delta_per_atom_eV"], sub[target])
        out[f"spearman_vs_{target}"] = {"n": len(sub), "rho": round(float(rho), 4), "p": round(float(p), 6)}

    sub = df.dropna(subset=["theoretical"])
    ct = pd.crosstab(sub["theoretical"], sub["bonding_label"])
    odds, p_fisher = fisher_exact(ct)
    chi2, p_chi2, _, _ = chi2_contingency(ct)
    out["bonding_label_vs_theoretical"] = {
        "contingency_table": ct.to_dict(),
        "fisher_exact_p": round(float(p_fisher), 4),
        "odds_ratio": round(float(odds), 4),
        "chi2": round(float(chi2), 4),
        "chi2_p": round(float(p_chi2), 4),
    }

    exp_ = df[df.theoretical == False]["delta_per_atom_eV"].dropna()  # noqa: E712
    theo_ = df[df.theoretical == True]["delta_per_atom_eV"].dropna()  # noqa: E712
    u, p_mw = mannwhitneyu(exp_, theo_)
    out["delta_per_atom_eV_experimental_vs_theoretical"] = {
        "n_experimental": len(exp_), "median_experimental": round(float(exp_.median()), 4),
        "n_theoretical": len(theo_), "median_theoretical": round(float(theo_.median()), 4),
        "mann_whitney_p": round(float(p_mw), 4),
    }

    fam3 = df[df.family.isin(["exp_stable", "exp_metastable", "theo_metastable"])]
    groups = [g["delta_per_atom_eV"].dropna().values for _, g in fam3.groupby("family")]
    h, p_kw = kruskal(*groups)
    out["kruskal_wallis_by_family"] = {
        "groups": ["exp_stable", "exp_metastable", "theo_metastable"],
        "H": round(float(h), 4), "p": round(float(p_kw), 4),
        "n_per_group": fam3.groupby("family").size().to_dict(),
        "median_delta_per_atom_eV_per_group": fam3.groupby("family")["delta_per_atom_eV"].median().round(4).to_dict(),
        "endobondic_fraction_per_group": (
            fam3.groupby("family")["bonding_label"].apply(lambda s: round(float((s == "endobondic").mean()), 4)).to_dict()
        ),
    }
    return out


def main():
    df = load_full_history()
    viability = viability_test(df)

    result = {"viability_prediction": viability}
    OUT_JSON.write_text(json.dumps(result, indent=2, default=str))

    print(f"n={viability['n_total']} case-1 reactions (full project history, all campaigns pooled)")
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
