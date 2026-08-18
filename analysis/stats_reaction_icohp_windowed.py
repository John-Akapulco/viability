"""Correlation analysis for the near-E_F windowed reaction Delta(ICOHP)
(analysis/reaction_icohp_windowed_case1.csv) against formation_energy_per_atom,
stratified by bond_type/is_metal -- same convention as every other
descriptor in this project (Spearman, n<15 flagged).

bond_type is re-sourced via build_bond_type_map() (not meta.get()
directly) -- same fix applied twice already this session for the same
class of bug (186 main-campaign compounds' own mp_metadata.json lacks
bond_type/is_metal).

Also loads formation_energy_per_atom from mp_dataset/formation_energies.json
(not in the windowed-case1 CSV itself) and compares against the two
other reaction-delta constructions already in the project:
  - reaction_icohp_case1.csv (full ICOHP, mission #5)
  - delta_antibonding_case1.csv (antibonding-only, this session)

Writes analysis/stats_summary_reaction_icohp_windowed.json.
"""

import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent
sys.path.insert(0, str(HERE))
WINDOWED_CSV = HERE / "reaction_icohp_windowed_case1.csv"
BONDTYPE_CSV = HERE / "icohp_icobi_bondtype.csv"
FULL_CSV = HERE / "reaction_icohp_case1.csv"
DELTA_ANTIBOND_CSV = HERE / "delta_antibonding_case1.csv"
FORMATION_ENERGIES = REPO_ROOT / "mp_dataset" / "formation_energies.json"

SMALL_N_THRESHOLD = 15


def spearman_row(df: pd.DataFrame, metric: str, target: str, group_label: str) -> dict:
    sub = df[[metric, target]].dropna()
    n = len(sub)
    if n < 4:
        return {"group": group_label, "metric": metric, "target": target, "n": n, "rho": None, "p_value": None}
    rho, p = spearmanr(sub[metric], sub[target])
    return {
        "group": group_label, "metric": metric, "target": target, "n": n,
        "rho": round(float(rho), 4), "p_value": round(float(p), 4),
        "note": "small sample (n<15)" if n < SMALL_N_THRESHOLD else None,
    }


def correlation_table(df: pd.DataFrame, metric: str, target: str) -> list[dict]:
    rows = [spearman_row(df, metric, target, "all")]
    for bond_type, sub in df.groupby("bond_type"):
        rows.append(spearman_row(sub, metric, target, f"bond_type={bond_type}"))
    for is_metal, sub in df.groupby("is_metal"):
        rows.append(spearman_row(sub, metric, target, f"is_metal={is_metal}"))
    return rows


def main() -> None:
    df = pd.read_csv(WINDOWED_CSV)
    formation_energies = json.loads(FORMATION_ENERGIES.read_text())
    df["formation_energy_per_atom"] = df["mp_id"].map(formation_energies)
    # is_metal-first ICOBI classifier (icohp_icobi_bondtype.csv), not
    # build_bond_type_map()'s classify()-heuristic-first mapping -- same
    # sourcing fix already applied elsewhere in the pipeline this session.
    bt = pd.read_csv(BONDTYPE_CSV)[["compound_id", "icobi_label", "is_metal"]]
    df = df.drop(columns=[c for c in ("bond_type", "is_metal") if c in df.columns])
    df = df.merge(bt, on="compound_id", how="left").rename(columns={"icobi_label": "bond_type"})

    summary = {
        "n_total": len(df),
        "vs_formation_energy_stratified": correlation_table(
            df, "delta_icohp_windowed_per_atom", "formation_energy_per_atom"),
        "vs_hull_stratified": correlation_table(
            df, "delta_icohp_windowed_per_atom", "energy_above_hull_eV_per_atom"),
    }

    # Cross-correlation against the two other reaction-delta constructions
    if FULL_CSV.exists():
        full = pd.read_csv(FULL_CSV)[["compound_id", "delta_icohp_per_atom"]]
        merged = df.merge(full, on="compound_id", how="inner").dropna(
            subset=["delta_icohp_windowed_per_atom", "delta_icohp_per_atom"])
        if len(merged) >= 4:
            rho, p = spearmanr(merged["delta_icohp_windowed_per_atom"], merged["delta_icohp_per_atom"])
            summary["vs_full_reaction_icohp"] = {"n": len(merged), "rho": round(float(rho), 4), "p_value": round(float(p), 4)}

    if DELTA_ANTIBOND_CSV.exists():
        antibond = pd.read_csv(DELTA_ANTIBOND_CSV)[["compound_id", "delta_icohp_antibond"]]
        merged = df.merge(antibond, on="compound_id", how="inner").dropna(
            subset=["delta_icohp_windowed_per_atom", "delta_icohp_antibond"])
        if len(merged) >= 4:
            rho, p = spearmanr(merged["delta_icohp_windowed_per_atom"], merged["delta_icohp_antibond"])
            summary["vs_delta_antibonding"] = {"n": len(merged), "rho": round(float(rho), 4), "p_value": round(float(p), 4)}

    (HERE / "stats_summary_reaction_icohp_windowed.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
