"""Does the antibonding-population descriptor (mission #4, raw ICOHP
antibonding near E_F/VBM) or its reaction-level Delta (mission #4b)
distinguish which polymorph of a same-composition group is
thermodynamically most stable?

Motivated by the dataset's own design (mp_dataset/download_extension.py's
track (ii): "build genuine same-composition polymorph groups and to
probe descriptor behavior far from equilibrium", manuscript.tex
Sec.~dataset) -- this test was never actually run against that design
goal before now. It is the antibonding-descriptor analog of Case 2 in
reaction_icohp.py/REPORT_reaction_icohp.md Sec.6 (does the strongest-
bonding polymorph track the most stable one? -- 47.3% agreement there,
indistinguishable from the 42.4% chance baseline).

Method, per icobi_label bond-type group (metallic/ionic/covalent/mixed):
1. Restrict to formulas with >=2 entries under that group.
2. Deduplicate by (formula, mp_id) keeping one row per physical
   structure -- the raw dataset contains real duplicate compound_ids for
   the same mp_id under different naming (e.g. an "_exp"/"_theo" suffix
   applied inconsistently across batches), which would otherwise inflate
   group sizes and bias the correlation.
3. Within each formula group, set the lowest-energy_above_hull member as
   the enthalpy zero (the "ground state" for that composition); every
   other member's dH_eV_at is its own energy_above_hull minus that
   minimum -- valid directly (no re-referencing needed) because
   energy_above_hull is already computed on a common basis and formula
   is held fixed within a group.
4. Spearman rho between (a) the raw antibond_w_raw and (b) the reaction
   delta_icohp_antibond (mission #4b, decomposition-into-elements; not
   defined for pure elements like carbon allotropes, or where a batch's
   elemental references are incomplete) against dH_eV_at, pooled across
   every formula in the group.

Writes analysis/polymorph_antibonding_correlation.csv (per-row data) and
analysis/stats_summary_polymorph_antibonding.json (per-group Spearman
results).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
BONDTYPE_CSV = HERE / "icohp_icobi_bondtype.csv"
ANTIBOND_CSV = HERE / "icohp_antibonding_full.csv"
DELTA_CSV = HERE / "delta_antibonding_case1.csv"
OUT_CSV = HERE / "polymorph_antibonding_correlation.csv"
OUT_JSON = HERE / "stats_summary_polymorph_antibonding.json"

BOND_GROUPS = ["ionic", "covalent", "mixed", "metallic"]


def build_group(merged: pd.DataFrame, bond: str) -> pd.DataFrame:
    sub_all = merged[merged["icobi_label"] == bond].copy()
    # Prefer keeping the row that has a delta_icohp_antibond value when a
    # duplicate (formula, mp_id) pair exists, so downstream reaction-delta
    # correlations lose as few rows as possible to the dedup step.
    sub_all = sub_all.sort_values("delta_icohp_antibond", na_position="last")
    sub_all = sub_all.drop_duplicates(subset=["formula", "mp_id"], keep="first")

    counts = sub_all["formula"].value_counts()
    poly_formulas = sorted(counts[counts > 1].index.tolist())
    sub = sub_all[sub_all["formula"].isin(poly_formulas)].copy()

    rows = []
    for f in poly_formulas:
        g = sub[sub["formula"] == f].sort_values("energy_above_hull_eV_at").copy()
        e0 = g["energy_above_hull_eV_at"].iloc[0]
        icohp0 = g["antibond_w_raw"].iloc[0]
        g["dH_eV_at"] = g["energy_above_hull_eV_at"] - e0
        g["pct_vs_ground_state"] = ((g["antibond_w_raw"] - icohp0) / icohp0 * 100) if icohp0 != 0 else float("nan")
        rows.append(g)
    return pd.concat(rows) if rows else sub


def correlate(full: pd.DataFrame) -> dict:
    out: dict = {"n_rows": len(full), "n_formulas": full["formula"].nunique()}

    rho, p = spearmanr(full["antibond_w_raw"], full["dH_eV_at"])
    out["raw_icohp_antibonding_vs_dH"] = {"n": len(full), "rho": round(float(rho), 4), "p": round(float(p), 4)}

    d = full.dropna(subset=["delta_icohp_antibond"])
    if len(d) >= 3:
        rho2, p2 = spearmanr(d["delta_icohp_antibond"], d["dH_eV_at"])
        out["delta_icohp_antibonding_reaction_vs_dH"] = {"n": len(d), "rho": round(float(rho2), 4), "p": round(float(p2), 4)}
    else:
        out["delta_icohp_antibonding_reaction_vs_dH"] = {"n": len(d), "note": "too few rows to test"}
    return out


def main() -> None:
    bt = pd.read_csv(BONDTYPE_CSV)
    raw = pd.read_csv(ANTIBOND_CSV)[["compound_id", "antibond_w_raw"]]
    delta = pd.read_csv(DELTA_CSV)[["compound_id", "delta_icohp_antibond"]]
    merged = bt.merge(raw, on="compound_id", how="left").merge(delta, on="compound_id", how="left")

    all_rows = []
    summary = {}
    for bond in BOND_GROUPS:
        full = build_group(merged, bond)
        full = full.assign(bond_type_group=bond)
        all_rows.append(full)
        summary[bond] = correlate(full)
        print(f"{bond}: n_formulas={full['formula'].nunique()} n_rows={len(full)} -- "
              f"raw rho={summary[bond]['raw_icohp_antibonding_vs_dH']['rho']} "
              f"p={summary[bond]['raw_icohp_antibonding_vs_dH']['p']}")

    combined = pd.concat(all_rows)
    keep_cols = [
        "bond_type_group", "formula", "compound_id", "mp_id", "spacegroup_symbol",
        "theoretical", "antibond_w_raw", "pct_vs_ground_state", "dH_eV_at", "delta_icohp_antibond",
    ]
    combined[keep_cols].to_csv(OUT_CSV, index=False)
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_CSV}\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
