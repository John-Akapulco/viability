"""Central hypothesis test of this project, as reformulated 2026-08-16:
does the reaction-ICOHP descriptor (delta(ICOHP), decomposition into
elements) distinguish thermodynamically stable, metastable, and unstable
compounds?

Two parts:

1. Concordance test against Reitz & Dronskowski (ic-2026-04181q): Lin's
   Concordance Correlation Coefficient (CCC) -- not just Pearson/Spearman
   correlation, which would be satisfied by any monotonic relationship
   even with a large systematic offset; CCC specifically tests agreement
   with the identity line, the right test for "do our numbers match the
   manuscript's," not just "do they move together." Computed on the 7
   worked reactions already validated in
   tests/test_reitz_dronskowski_validation.py.

2. Viability prediction test, extended to every element and compound
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

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, kruskal, mannwhitneyu, pearsonr, spearmanr

REPO_ROOT = Path(__file__).parent.parent
OUT_JSON = Path(__file__).parent / "stats_summary_delta_icohp_viability.json"


def lin_ccc(x: np.ndarray, y: np.ndarray) -> float:
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(), y.var()
    sxy = ((x - mx) * (y - my)).mean()
    return float(2 * sxy / (vx + vy + (mx - my) ** 2))


def manuscript_concordance() -> dict:
    # (reaction_id, manuscript kJ/mol, computed kJ/mol) -- computed values
    # from tests/test_reitz_dronskowski_validation.py, reproduced here as
    # literal numbers (not re-derived) since this script's job is the
    # concordance statistic, not re-running the fixture validation.
    cases = [
        ("Pb(N3)2 -> Pb + 3N2", 1345, 1344.3),
        ("S4N2 -> 1/2 S8 + N2", 258, 258.5),
        ("S4N4 -> 1/2 S8 + 2N2", 399, 397.3),
        ("ZnSn -> Zn + Sn", -337, -336.6),
        ("CaO[sphalerite] -> CaO[rocksalt]", -79, -78.9),
        ("CaN -> 1/3 Ca3N2 + 1/6 N2", -205, -204.8),
        ("Mn2O7 -> 2MnO2 + 3/2 O2", -186, -187.7),
    ]
    manuscript = np.array([c[1] for c in cases], dtype=float)
    computed = np.array([c[2] for c in cases], dtype=float)
    ccc = lin_ccc(manuscript, computed)
    pear_r, pear_p = pearsonr(manuscript, computed)
    sign_agree = int(np.sum(np.sign(manuscript) == np.sign(computed)))
    resid = computed - manuscript
    return {
        "n": len(cases),
        "cases": [{"reaction": c[0], "manuscript_kJ_per_mol": c[1], "computed_kJ_per_mol": c[2]} for c in cases],
        "lin_ccc": round(ccc, 6),
        "pearson_r": round(float(pear_r), 6),
        "pearson_p": pear_p,
        "sign_agreement": f"{sign_agree}/{len(cases)}",
        "mean_abs_residual_kJ_per_mol": round(float(np.mean(np.abs(resid))), 3),
        "max_abs_residual_kJ_per_mol": round(float(np.max(np.abs(resid))), 3),
        "rmse_kJ_per_mol": round(float(np.sqrt(np.mean(resid ** 2))), 3),
    }


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
    concordance = manuscript_concordance()
    df = load_full_history()
    viability = viability_test(df)

    result = {"manuscript_concordance": concordance, "viability_prediction": viability}
    OUT_JSON.write_text(json.dumps(result, indent=2, default=str))

    print(f"Manuscript concordance: Lin CCC={concordance['lin_ccc']}, sign agreement={concordance['sign_agreement']}")
    print(f"n={viability['n_total']} case-1 reactions (full project history, all campaigns pooled)")
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
