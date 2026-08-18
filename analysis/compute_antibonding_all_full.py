"""Compute the ICOHP antibonding-population-near-frontier descriptor
(cohp_extraction.antibonding_population_near_frontier, mission #4's
original metric -- REPORT_antibonding.md) across every compound with a
COHPCAR.lobster file, regardless of which batch it came from.

Supersedes the batch-specific top-up scripts (compute_antibonding_all.py:
merged against percolation_vs_formation_energy.csv, 349-compound cap;
compute_antibonding_extension.py and compute_icohp_antibonding_maxhull.py:
one-off 22- and 45-compound batches) whose outputs were manually
concatenated at some point into analysis/icohp_antibonding_full.csv (no
generator script for that concatenation exists in the repo) -- that file
was last regenerated at the 394-compound (maxhull-inclusive) scale and
had not been extended to the marginal-formation-energy/widen batches
added since. This script directly overwrites
analysis/icohp_antibonding_full.csv with one consistent, full-scale
(597-directory) computation instead of another one-off top-up, following
the same scan-every-compound-directly pattern as
compute_icobi_antibonding_all.py (the ICOBI sibling of this metric).

is_metal sourcing follows the same precedence as
compute_icobi_antibonding_all.py: analysis/percolation_vs_antibonding.csv
where available (Materials-Project-fetched, not derived locally -- see
compute_antibonding_all.py's AlNi/BeCu pitfall), mp_metadata.json
otherwise.

Writes analysis/icohp_antibonding_full.csv (one row per compound) and
analysis/icohp_antibonding_full.json (per-compound, all delta_e values).
"""

import json
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
import cohp_extraction as ce  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
PERCOLATION_CSV = Path(__file__).parent / "percolation_vs_antibonding.csv"
FORMATION_ENERGIES = REPO_ROOT / "mp_dataset" / "formation_energies.json"
OUT_JSON = Path(__file__).parent / "icohp_antibonding_full.json"
OUT_CSV = Path(__file__).parent / "icohp_antibonding_full.csv"

DELTA_ES = (0.5, 1.0, 2.0)
PRIMARY_DELTA_E = 1.0


def build_is_metal_map() -> dict[str, bool]:
    is_metal = {}
    if PERCOLATION_CSV.exists():
        base = pd.read_csv(PERCOLATION_CSV)
        for _, row in base.iterrows():
            if pd.notna(row.get("is_metal")):
                is_metal[row["compound_id"]] = bool(row["is_metal"])
    for d in sorted(STRUCTURES_ROOT.iterdir()):
        if d.name in is_metal:
            continue
        meta_path = d / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("is_metal") is not None:
            is_metal[d.name] = bool(meta["is_metal"])
    return is_metal


def build_bond_type_map() -> dict[str, str]:
    bond_type = {}
    if PERCOLATION_CSV.exists():
        base = pd.read_csv(PERCOLATION_CSV)
        for _, row in base.iterrows():
            if pd.notna(row.get("bond_type")):
                bond_type[row["compound_id"]] = row["bond_type"]
    for d in sorted(STRUCTURES_ROOT.iterdir()):
        if d.name in bond_type:
            continue
        meta_path = d / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("bond_type") is not None:
            bond_type[d.name] = meta["bond_type"]
    return bond_type


def main():
    is_metal_map = build_is_metal_map()
    bond_type_map = build_bond_type_map()
    formation_energies = json.loads(FORMATION_ENERGIES.read_text())

    compound_dirs = [
        d for d in sorted(STRUCTURES_ROOT.iterdir())
        if (d / "COHPCAR.lobster").exists() and (d / "mp_metadata.json").exists()
    ]
    print(f"{len(compound_dirs)} compounds have COHPCAR.lobster")

    results = {}
    rows = []
    n_ok = n_failed = n_no_is_metal = 0
    t_start = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, d in enumerate(compound_dirs):
            compound_id = d.name
            meta = json.loads((d / "mp_metadata.json").read_text())
            is_metal = is_metal_map.get(compound_id)
            if is_metal is None:
                results[compound_id] = {"error": "no is_metal available"}
                n_no_is_metal += 1
                continue

            try:
                per_delta = {}
                for delta_e in DELTA_ES:
                    r = ce.antibonding_population_near_frontier(
                        d, is_metal=is_metal, delta_e=delta_e, vasprun_path=d / "vasprun.xml"
                    )
                    per_delta[str(delta_e)] = r
                results[compound_id] = {"is_metal": is_metal, "by_delta_e": per_delta, "error": None}
                n_ok += 1

                primary = per_delta[str(PRIMARY_DELTA_E)]
                mp_id = meta.get("mp_id")
                rows.append({
                    "compound_id": compound_id,
                    "mp_id": mp_id,
                    "formula": meta.get("formula"),
                    "family": meta.get("family"),
                    "bond_type": bond_type_map.get(compound_id),
                    "is_metal": is_metal,
                    "theoretical": meta.get("theoretical"),
                    "energy_above_hull_eV_per_atom": meta.get("energy_above_hull_eV_per_atom"),
                    "formation_energy_per_atom": formation_energies.get(mp_id) if mp_id else None,
                    "antibond_w_raw": primary["w_antibond_raw"],
                    "antibond_w_normalized": primary["w_antibond_normalized"],
                    "antibond_e_ref": primary["e_ref"],
                    "antibond_w_raw_dE0.5": per_delta["0.5"]["w_antibond_raw"],
                    "antibond_w_raw_dE1.0": per_delta["1.0"]["w_antibond_raw"],
                    "antibond_w_raw_dE2.0": per_delta["2.0"]["w_antibond_raw"],
                })
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                results[compound_id] = {"is_metal": is_metal, "error": f"{type(exc).__name__}: {exc}"}
                n_failed += 1

            if (i + 1) % 50 == 0:
                elapsed = time.time() - t_start
                print(f"  {i+1}/{len(compound_dirs)} done ({elapsed:.0f}s elapsed, ok={n_ok} failed={n_failed})")

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWrote {len(results)} entries to {OUT_JSON} "
          f"(ok={n_ok}, failed={n_failed}, no_is_metal={n_no_is_metal}, {time.time()-t_start:.0f}s total)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
