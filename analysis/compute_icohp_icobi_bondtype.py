"""Independent bonding-type classification from ICOHP/ICOBI data itself,
cross-checked against classify() (fetch_candidates.py, composition +
is_metal only, reused unmodified throughout this project): does the
electronic-structure data agree with the composition-based heuristic,
and does it cover the cases classify() leaves unclassified?

Reuses analysis/percolation_vs_antibonding.csv's per-compound icobi_mean
(mean ICOBI per bond, dimensionless bond order, already computed by
percolation_path.py for the whole 349-compound dataset -- no LOBSTER
re-parsing needed) and is_metal (already fetched from Materials Project
for every compound, including the 186 whose own mp_metadata.json lacks
it, see compute_antibonding_all.py's docstring on that pitfall).

Scheme (is_metal first, since it is the one DFT-derived, unambiguous
input available for every compound; ICOHP/ICOBI magnitude is not a
reliable metal/non-metal discriminant on its own -- some transition-metal
carbides/nitrides show strong ICOHP despite being metallic):
  1. is_metal=True -> "metallic".
  2. Otherwise, mean ICOBI per bond measures covalency/orbital overlap
     directly: high ICOBI -> shared, covalent bonding; low ICOBI ->
     largely electrostatic, ionic bonding. Threshold = midpoint of the
     median ICOBI among compounds classify() itself already labels
     ionic vs covalent (calibrated from data, not chosen a priori).
  3. No ICOBI data -> "not_classified".

Writes analysis/icohp_icobi_bondtype.csv and prints the calibration and
the full classify_label x icobi_label cross-tabulation.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
IN_CSV = REPO_ROOT / "analysis" / "percolation_vs_antibonding.csv"
OUT_CSV = REPO_ROOT / "analysis" / "icohp_icobi_bondtype.csv"


def main() -> None:
    df = pd.read_csv(IN_CSV)
    df = df.dropna(subset=["is_metal"])

    ionic = df[df["bond_type"] == "ionic"]["icobi_mean"].dropna()
    covalent = df[df["bond_type"] == "covalent"]["icobi_mean"].dropna()
    med_ionic, med_covalent = ionic.median(), covalent.median()
    threshold = (med_ionic + med_covalent) / 2
    print(f"Calibration: classify()=ionic n={len(ionic)} median ICOBI/bond={med_ionic:.4f}")
    print(f"             classify()=covalent n={len(covalent)} median ICOBI/bond={med_covalent:.4f}")
    print(f"Threshold (midpoint of medians): {threshold:.4f}")

    def icobi_label(row):
        if row["is_metal"]:
            return "metallic"
        if pd.isna(row["icobi_mean"]):
            return "not_classified"
        return "covalent" if row["icobi_mean"] >= threshold else "ionic"

    df["icobi_label"] = df.apply(icobi_label, axis=1)

    labeled = df[df["bond_type"].isin(["ionic", "covalent", "metallic"])]
    agree = (labeled["bond_type"] == labeled["icobi_label"]).sum()
    print(f"\nConcordance on classify()-labeled subset: {agree}/{len(labeled)} = {agree/len(labeled):.1%}")

    print("\nCross-tab (classify_label rows x icobi_label cols):")
    print(pd.crosstab(df["bond_type"].fillna("unclassified"), df["icobi_label"]))

    unclassified = df[df["bond_type"].isna()]
    newly = unclassified[unclassified["icobi_label"].isin(["ionic", "covalent"])]
    print(f"\nOf {len(unclassified)} classify()-unclassified compounds, {len(newly)} "
          f"get an ionic/covalent label from ICOBI: {newly['icobi_label'].value_counts().to_dict()}")

    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
