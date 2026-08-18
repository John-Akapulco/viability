"""Independent bonding-type classification from ICOHP/ICOBI data itself,
cross-checked against classify() (fetch_candidates.py, composition +
is_metal only, reused unmodified throughout this project).

REWRITTEN 2026-08-17 to fix a real methodological bug, user-reported
after seeing diamond, boron, SiO2, C3N4, Cl2, ClF3 and HCl listed under
"Ioniques" in the SI master list: the previous version's `icobi_mean`
was reused from analysis/percolation_vs_antibonding.csv, which averages
ICOBI over EVERY symmetry-inequivalent bond label LOBSTER reports within
the wide cohpGenerator cutoff (0.1-6.0 A, see report Sec. 2.4) --
correct for percolation-graph connectivity (percolation_path.py's actual
purpose) but catastrophic as a bonding-type descriptor: diamond alone
has 4012 ICOBILIST entries within that cutoff, of which only 4 are the
real nearest-neighbor C-C bond (ICOBI~0.95); the other 4008 are long-
range near-zero pairs that drag the naive mean down to 0.017, well
below any reasonable covalent threshold. Verified directly on
ICOBILIST.lobster before writing this fix (see project memory).

Fix: classify from the STRONGEST first-coordination-shell bond-type
population instead of a flat mean over every LOBSTER-reported pair.
Uses reaction_analysis.parse_lobster.parse_compound_entry(bond_filter=
"nearest_neighbor") (Reitz & Dronskowski's first-shell convention,
already implemented + tested elsewhere in this project) to get
by_bond_type (mean ICOBI per species pair, first shell only), then:

  icobi_primary_mean = max(by_bond_type means)  -- the strongest
  identified bond in the structure.

Second, independently reported bug (also user-flagged): CsN3, NaN3,
LiP5, NaS and similar compounds are Zintl-type solids -- a covalent
polyatomic anion (N3-, P5-chain, S-S dimer) held together by a much
weaker, genuinely ionic link to the counter-cation. Forcing a single
ionic/covalent label is chemically wrong for these regardless of how
the "primary" bond is computed (whichever population wins the max(),
the other one is real too and gets silently discarded). Detected via a
structural signature, not another distance/magnitude threshold: a
HOMOATOMIC pair (e.g. N-N, P-P, S-S) with mean ICOBI >= COVALENT_
THRESHOLD (a genuine intra-anion covalent bond) coexisting with a
HETEROATOMIC pair (e.g. Na-N, Li-P, connecting a different element)
whose mean is measurably weaker -- exactly the Zintl-Klemm picture
(strong homoatomic polyanion backbone + weaker link to a distinct
cation). Compounds meeting this get a 4th label, "mixed", instead of
being forced into ionic or covalent.

Scheme (is_metal first, since it is the one DFT-derived, unambiguous
input available for every compound):
  1. is_metal=True -> "metallic".
  2. Homoatomic pair >= COVALENT_THRESHOLD AND a weaker heteroatomic
     pair also present -> "mixed" (Zintl-type).
  3. Otherwise, icobi_primary_mean >= COVALENT_THRESHOLD -> "covalent";
     < COVALENT_THRESHOLD -> "ionic". Threshold = midpoint of the
     median icobi_primary_mean among compounds classify() itself
     already labels ionic vs covalent (calibrated from data, same
     methodology as before, just on the corrected metric).
  4. No ICOBI data -> "not_classified".

Scans mp_dataset/structures/* directly (every compound with both
ICOHPLIST.lobster and ICOBILIST.lobster present, any batch) instead of
depending on analysis/percolation_vs_antibonding.csv's batch coverage --
avoids the recurring "which batch is/isn't in the base CSV" bug class
(see project memory: this exact bug recurred three times in three
different scripts). marginal_formation_energy_campaign and
marginal_widen_campaign are now included (both campaigns finished,
194 compounds fully computed, already folded into reaction_icohp_case1.csv
and reaction_analysis_case1_full.csv -- excluding them here left those
two files' bond_type column ~70% NaN for no reason tied to data
availability).

Writes analysis/icohp_icobi_bondtype.csv and prints the calibration and
the full classify_label x icobi_label cross-tabulation.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from reaction_analysis.parse_lobster import parse_compound_entry  # noqa: E402
from analysis.compute_icobi_antibonding_all import build_is_metal_map, build_bond_type_map  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT_CSV = REPO_ROOT / "analysis" / "icohp_icobi_bondtype.csv"
ICOHP_ANTIBOND_FULL_CSV = REPO_ROOT / "analysis" / "icohp_antibonding_full.csv"

EXCLUDED_BATCHES: set[str] = set()


def _element_pair(key: str) -> tuple[str, str]:
    a, b = key.split("-")
    return a, b


def _classify_compound(entry) -> dict:
    """Given a CompoundEntry parsed with bond_filter='nearest_neighbor',
    return icobi_primary_mean/pair, icohp_primary_mean/pair, is_mixed_zintl."""
    result = {
        "icobi_primary_mean": None, "icobi_primary_pair": None,
        "icohp_primary_mean": None, "icohp_primary_pair": None,
        "is_mixed_zintl": False,
    }
    if entry.icobi is not None and entry.icobi.by_bond_type:
        by_pair = entry.icobi.by_bond_type
        primary_pair, primary = max(by_pair.items(), key=lambda kv: kv[1].mean_eV)
        result["icobi_primary_mean"] = primary.mean_eV
        result["icobi_primary_pair"] = primary_pair

        homoatomic_high = [
            (pair, s) for pair, s in by_pair.items()
            if _element_pair(pair)[0] == _element_pair(pair)[1] and s.mean_eV >= COVALENT_THRESHOLD
        ]
        if homoatomic_high:
            strongest_homo_mean = max(s.mean_eV for _, s in homoatomic_high)
            weaker_heteroatomic = [
                (pair, s) for pair, s in by_pair.items()
                if _element_pair(pair)[0] != _element_pair(pair)[1] and s.mean_eV < strongest_homo_mean
            ]
            if weaker_heteroatomic:
                result["is_mixed_zintl"] = True

    if entry.icohp.by_bond_type:
        by_pair = entry.icohp.by_bond_type
        primary_pair, primary = max(by_pair.items(), key=lambda kv: abs(kv[1].mean_eV))
        result["icohp_primary_mean"] = primary.mean_eV
        result["icohp_primary_pair"] = primary_pair

    return result


# Provisional threshold, overwritten by calibration in main() before use --
# module-level so _classify_compound (called during the first, uncalibrated
# pass) has something to compare against; the mixed-detection pass is
# re-run after calibration with the real value.
COVALENT_THRESHOLD = 0.25


def _scan_all_compounds(is_metal_map: dict, bond_type_map: dict) -> pd.DataFrame:
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
            meta_path = compound_dir / "mp_metadata.json"
            if not meta_path.exists():
                continue
            meta = json.loads(meta_path.read_text())
            if meta.get("batch") in EXCLUDED_BATCHES:
                continue
            if meta.get("quality_excluded"):
                continue
            if not (compound_dir / "ICOHPLIST.lobster").exists():
                continue
            is_metal = is_metal_map.get(compound_dir.name)
            if is_metal is None:
                continue
            try:
                entry = parse_compound_entry(
                    compound_dir, role="target", compound_id=compound_dir.name,
                    bond_filter="nearest_neighbor",
                )
            except Exception:  # noqa: BLE001 - one bad compound must not kill the batch
                continue
            classified = _classify_compound(entry)
            rows.append({
                "compound_id": compound_dir.name,
                "mp_id": meta.get("mp_id"),
                "formula": meta.get("formula") or entry.formula,
                "bond_type": bond_type_map.get(compound_dir.name),
                "family": meta.get("family"),
                "theoretical": meta.get("theoretical"),
                "energy_above_hull_eV_at": meta.get("energy_above_hull_eV_per_atom"),
                "spacegroup_symbol": meta.get("spacegroup"),
                "is_metal": is_metal,
                "formation_energy_per_atom": meta.get("formation_energy_per_atom"),
                **classified,
                # Kept under the ORIGINAL column names for downstream
                # consumers (report/gen_appendix_master_list.py): these
                # now hold the corrected first-shell PRIMARY bond value,
                # not the old diluted percolation-graph mean.
                "icobi_mean": classified["icobi_primary_mean"],
                "icohp_mean": classified["icohp_primary_mean"],
            })
    df = pd.DataFrame(rows)

    if ICOHP_ANTIBOND_FULL_CSV.exists():
        antibond = pd.read_csv(ICOHP_ANTIBOND_FULL_CSV)[["compound_id", "antibond_w_normalized"]]
        df = df.merge(antibond, on="compound_id", how="left")

    return df


def main() -> None:
    global COVALENT_THRESHOLD

    is_metal_map = build_is_metal_map()
    bond_type_map = build_bond_type_map()

    # Pass 1: uncalibrated (provisional threshold) just to get every
    # compound's icobi_primary_mean for calibration.
    df = _scan_all_compounds(is_metal_map, bond_type_map)
    print(f"Scanned {len(df)} compounds with usable LOBSTER data (excluded batches: {EXCLUDED_BATCHES})")

    ionic = df[df["bond_type"] == "ionic"]["icobi_primary_mean"].dropna()
    covalent = df[df["bond_type"] == "covalent"]["icobi_primary_mean"].dropna()
    med_ionic, med_covalent = ionic.median(), covalent.median()
    COVALENT_THRESHOLD = (med_ionic + med_covalent) / 2
    print(f"Calibration (corrected first-shell primary-bond ICOBI): classify()=ionic n={len(ionic)} median={med_ionic:.4f}")
    print(f"             classify()=covalent n={len(covalent)} median={med_covalent:.4f}")
    print(f"Threshold (midpoint of medians): {COVALENT_THRESHOLD:.4f}")

    # Pass 2: recompute mixed-Zintl detection + final labels with the
    # calibrated threshold (the provisional one above was only a
    # placeholder for pass 1's own homoatomic-high check).
    df = _scan_all_compounds(is_metal_map, bond_type_map)

    def icobi_label(row):
        if row["is_metal"]:
            return "metallic"
        if pd.isna(row["icobi_primary_mean"]):
            return "not_classified"
        if row["is_mixed_zintl"]:
            return "mixed"
        return "covalent" if row["icobi_primary_mean"] >= COVALENT_THRESHOLD else "ionic"

    df["icobi_label"] = df.apply(icobi_label, axis=1)

    labeled = df[df["bond_type"].isin(["ionic", "covalent", "metallic"])]
    agree = (labeled["bond_type"] == labeled["icobi_label"]).sum()
    print(f"\nConcordance on classify()-labeled subset: {agree}/{len(labeled)} = {agree/len(labeled):.1%}")

    print("\nCross-tab (classify_label rows x icobi_label cols):")
    print(pd.crosstab(df["bond_type"].fillna("unclassified"), df["icobi_label"]))

    n_mixed = (df["icobi_label"] == "mixed").sum()
    print(f"\n{n_mixed} compounds flagged mixed (Zintl-type): "
          f"{sorted(df[df['icobi_label']=='mixed']['formula'].unique().tolist())}")

    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(df)} rows to {OUT_CSV}")

    threshold_path = REPO_ROOT / "analysis" / "icohp_icobi_bondtype_threshold.json"
    threshold_path.write_text(json.dumps({"covalent_threshold_icobi_primary_mean": COVALENT_THRESHOLD}, indent=2))
    print(f"Wrote calibrated threshold ({COVALENT_THRESHOLD:.4f}) to {threshold_path}")


if __name__ == "__main__":
    main()
