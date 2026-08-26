"""Does the antibonding-population descriptor -- ICOHP-based (mission #4)
or ICOBI-based, raw or as a reaction-level Delta (mission #4b) --
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

ICOBI is tested alongside ICOHP, never in isolation, matching this
project's convention throughout mission #4b
(analysis/compute_delta_antibonding_case1.py) -- the two variants run in
formally opposite raw-population sign conventions (README.md, "Antibonding
population near the frontier") but that does not by itself predict how
each behaves on this polymorph-ranking question, which is why both are
tested rather than assuming ICOBI would simply mirror ICOHP.

A fifth group, "elemental", pools every single-element formula (allotropes)
regardless of icobi_label -- pure elements get assigned a bond_type by the
same is_metal-first logic as compounds (Section~bondtype), which is not
always consistent across allotropes of the same element (e.g. carbon:
graphite alone classifies "metallic" via is_metal=True, every other
allotrope "covalent"), so pooling by element rather than by bond_type is
the physically meaningful grouping here. Only 3 elements in the current
dataset have >=2 allotropes (C: 6, Sn: 2, Zn: 2); none has a case-1
reaction (an element cannot decompose into "other elements"), so only the
raw-descriptor correlation is defined for this group.

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
4. Spearman rho between each of the four antibonding quantities (raw
   ICOHP, raw ICOBI, reaction Delta(ICOHP), reaction Delta(ICOBI) --
   the reaction deltas undefined for pure elements like carbon allotropes,
   or where a batch's elemental references are incomplete) and dH_eV_at,
   pooled across every formula in the group.

Writes analysis/polymorph_antibonding_correlation.csv (per-row data) and
analysis/stats_summary_polymorph_antibonding.json (per-group Spearman
results).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).parent
BONDTYPE_CSV = HERE / "icohp_icobi_bondtype.csv"
ICOHP_ANTIBOND_CSV = HERE / "icohp_antibonding_full.csv"
ICOBI_ANTIBOND_CSV = HERE / "icobi_antibonding_all.csv"
DELTA_CSV = HERE / "delta_antibonding_case1.csv"
OUT_CSV = HERE / "polymorph_antibonding_correlation.csv"
OUT_JSON = HERE / "stats_summary_polymorph_antibonding.json"

BOND_GROUPS = ["ionic", "covalent", "mixed", "metallic", "elemental"]


def is_elemental_formula(formula: str) -> bool:
    """True if `formula` (e.g. "C", "Sn", "Fe2O3") names a single element --
    exactly one distinct element symbol, any subscript."""
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", str(formula))
    elements = {sym for sym, _count in tokens if sym}
    return len(elements) == 1

# (raw column, delta column, pct output column, label used in stats keys)
DESCRIPTORS = [
    ("antibond_w_raw", "delta_icohp_antibond", "pct_vs_ground_state_icohp", "icohp"),
    ("icobi_antibond_w_raw", "delta_icobi_antibond", "pct_vs_ground_state_icobi", "icobi"),
]


def build_group(merged: pd.DataFrame, bond: str) -> pd.DataFrame:
    if bond == "elemental":
        sub_all = merged[merged["formula"].apply(is_elemental_formula)].copy()
    else:
        sub_all = merged[merged["icobi_label"] == bond].copy()
    # Prefer keeping the row that has both delta values when a duplicate
    # (formula, mp_id) pair exists, so downstream reaction-delta
    # correlations lose as few rows as possible to the dedup step.
    sub_all = sub_all.sort_values(["delta_icohp_antibond", "delta_icobi_antibond"], na_position="last")
    sub_all = sub_all.drop_duplicates(subset=["formula", "mp_id"], keep="first")

    counts = sub_all["formula"].value_counts()
    poly_formulas = sorted(counts[counts > 1].index.tolist())
    sub = sub_all[sub_all["formula"].isin(poly_formulas)].copy()

    rows = []
    for f in poly_formulas:
        g = sub[sub["formula"] == f].sort_values("energy_above_hull_eV_at").copy()
        e0 = g["energy_above_hull_eV_at"].iloc[0]
        g["dH_eV_at"] = g["energy_above_hull_eV_at"] - e0
        for raw_col, _, pct_col, _ in DESCRIPTORS:
            v0 = g[raw_col].iloc[0]
            g[pct_col] = ((g[raw_col] - v0) / v0 * 100) if v0 != 0 else float("nan")
        rows.append(g)
    return pd.concat(rows) if rows else sub


def correlate(full: pd.DataFrame) -> dict:
    out: dict = {"n_rows": len(full), "n_formulas": full["formula"].nunique()}

    for raw_col, delta_col, _, label in DESCRIPTORS:
        rho, p = spearmanr(full[raw_col], full["dH_eV_at"])
        out[f"raw_{label}_antibonding_vs_dH"] = {"n": len(full), "rho": round(float(rho), 4), "p": round(float(p), 4)}

        d = full.dropna(subset=[delta_col])
        if len(d) >= 3:
            rho2, p2 = spearmanr(d[delta_col], d["dH_eV_at"])
            out[f"delta_{label}_antibonding_reaction_vs_dH"] = {"n": len(d), "rho": round(float(rho2), 4), "p": round(float(p2), 4)}
        else:
            out[f"delta_{label}_antibonding_reaction_vs_dH"] = {"n": len(d), "note": "too few rows to test"}
    return out


def main() -> None:
    bt = pd.read_csv(BONDTYPE_CSV)
    raw_icohp = pd.read_csv(ICOHP_ANTIBOND_CSV)[["compound_id", "antibond_w_raw"]]
    raw_icobi = pd.read_csv(ICOBI_ANTIBOND_CSV)[["compound_id", "icobi_antibond_w_raw"]]
    delta = pd.read_csv(DELTA_CSV)[["compound_id", "delta_icohp_antibond", "delta_icobi_antibond"]]
    merged = (
        bt.merge(raw_icohp, on="compound_id", how="left")
        .merge(raw_icobi, on="compound_id", how="left")
        .merge(delta, on="compound_id", how="left")
    )

    all_rows = []
    summary = {}
    for bond in BOND_GROUPS:
        full = build_group(merged, bond)
        full = full.assign(bond_type_group=bond)
        all_rows.append(full)
        summary[bond] = correlate(full)
        s = summary[bond]
        print(f"{bond}: n_formulas={full['formula'].nunique()} n_rows={len(full)} -- "
              f"ICOHP raw rho={s['raw_icohp_antibonding_vs_dH']['rho']} p={s['raw_icohp_antibonding_vs_dH']['p']}  "
              f"ICOBI raw rho={s['raw_icobi_antibonding_vs_dH']['rho']} p={s['raw_icobi_antibonding_vs_dH']['p']}")

    combined = pd.concat(all_rows)
    keep_cols = [
        "bond_type_group", "formula", "compound_id", "mp_id", "spacegroup_symbol", "theoretical",
        "antibond_w_raw", "pct_vs_ground_state_icohp", "delta_icohp_antibond",
        "icobi_antibond_w_raw", "pct_vs_ground_state_icobi", "delta_icobi_antibond",
        "dH_eV_at",
    ]
    combined[keep_cols].to_csv(OUT_CSV, index=False)
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_CSV}\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
