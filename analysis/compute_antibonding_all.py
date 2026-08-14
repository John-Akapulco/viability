"""Extend the antibonding-population-near-frontier metric (validated on
the 6 pilots, see METRIC_DEFINITION_antibonding.md) to the full 186-compound
dataset. Explicitly authorized follow-up to the pilot validation -- not run
automatically as part of any other pipeline.

is_metal for the reference-energy choice comes from Materials Project
(single batched query, same convention as elsewhere in this project) --
not derived locally, per the documented AlNi/BeCu pitfall.

Writes analysis/antibonding_all.json (per-compound, all delta_e values)
and folds delta_e=1.0 (the primary window) into
analysis/percolation_vs_formation_energy.csv -> analysis/percolation_vs_antibonding.csv.
"""

import json
import os
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
import cohp_extraction as ce  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
PERCOLATION_CSV = Path(__file__).parent / "percolation_vs_formation_energy.csv"
OUT_JSON = Path(__file__).parent / "antibonding_all.json"
OUT_CSV = Path(__file__).parent / "percolation_vs_antibonding.csv"

DELTA_ES = (0.5, 1.0, 2.0)
PRIMARY_DELTA_E = 1.0


def fetch_is_metal(mp_ids: list[str]) -> dict[str, bool]:
    from mp_api.client import MPRester

    api_key = open(os.path.expanduser("~/.mp_api_key")).read().strip()
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(material_ids=mp_ids, fields=["material_id", "is_metal"])
    return {str(d.material_id): d.is_metal for d in docs}


def main():
    base = pd.read_csv(PERCOLATION_CSV)
    mp_ids = sorted(base["mp_id"].dropna().unique())
    print(f"Fetching is_metal for {len(mp_ids)} compounds...")
    is_metal_map = fetch_is_metal(mp_ids)
    n_missing = len(set(mp_ids) - set(is_metal_map))
    if n_missing:
        print(f"WARNING: {n_missing} mp_id missing is_metal")

    results = {}
    n_ok = n_failed = 0
    t_start = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, row in base.iterrows():
            compound_id = row["compound_id"]
            mp_id = row["mp_id"]
            is_metal = is_metal_map.get(mp_id)
            if is_metal is None:
                results[compound_id] = {"error": "no is_metal from MP"}
                n_failed += 1
                continue

            d = STRUCTURES_ROOT / compound_id
            try:
                per_delta = {}
                for delta_e in DELTA_ES:
                    r = ce.antibonding_population_near_frontier(
                        d, is_metal=is_metal, delta_e=delta_e, vasprun_path=d / "vasprun.xml"
                    )
                    per_delta[str(delta_e)] = r
                results[compound_id] = {"is_metal": is_metal, "by_delta_e": per_delta, "error": None}
                n_ok += 1
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                results[compound_id] = {"is_metal": is_metal, "error": f"{type(exc).__name__}: {exc}"}
                n_failed += 1

            if (i + 1) % 20 == 0:
                elapsed = time.time() - t_start
                print(f"  {i+1}/{len(base)} done ({elapsed:.0f}s elapsed, ok={n_ok} failed={n_failed})")

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWrote {len(results)} entries to {OUT_JSON} (ok={n_ok}, failed={n_failed}, "
          f"{time.time()-t_start:.0f}s total)")

    rows = []
    for _, row in base.iterrows():
        compound_id = row["compound_id"]
        r = results.get(compound_id, {})
        extra = {"compound_id": compound_id}
        if r.get("error") is None and "by_delta_e" in r:
            primary = r["by_delta_e"][str(PRIMARY_DELTA_E)]
            extra["is_metal"] = r["is_metal"]
            extra["antibond_w_raw"] = primary["w_antibond_raw"]
            extra["antibond_w_normalized"] = primary["w_antibond_normalized"]
            extra["antibond_e_ref"] = primary["e_ref"]
            for delta_e in DELTA_ES:
                extra[f"antibond_w_raw_dE{delta_e}"] = r["by_delta_e"][str(delta_e)]["w_antibond_raw"]
        else:
            extra["antibond_error"] = r.get("error", "unknown")
        rows.append(extra)

    extra_df = pd.DataFrame(rows)
    merged = base.merge(extra_df, on="compound_id", how="left")
    merged.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(merged)} rows to {OUT_CSV}")
    print(merged["antibond_error"].value_counts(dropna=False) if "antibond_error" in merged else "no errors column")


if __name__ == "__main__":
    main()
