"""Extend the ICOHP antibonding-population-near-frontier descriptor
(cohp_extraction.antibonding_population_near_frontier, mission #4) to the
maxhull_binaries_stress_test batch, which predates compute_antibonding_all.py
and isn't in analysis/percolation_vs_antibonding.csv (that CSV's own base,
percolation_vs_formation_energy.csv, only covers the pre-maxhull 349 --
this batch was never run through percolation_path.py). Needed to compare
the ICOHP and ICOBI antibonding descriptors side by side on the SAME,
full 394-compound population rather than two different-sized ones.

Writes analysis/icohp_antibonding_maxhull.csv (45 rows, same columns as
the icobi_antibond_* ones in icobi_antibonding_all.csv so the two can be
merged directly).
"""

import json
import sys
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
import cohp_extraction as ce  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
FORMATION_ENERGIES = REPO_ROOT / "mp_dataset" / "formation_energies.json"
OUT_CSV = Path(__file__).parent / "icohp_antibonding_maxhull.csv"
BATCH = "maxhull_binaries_stress_test"
DELTA_ES = (0.5, 1.0, 2.0)
PRIMARY_DELTA_E = 1.0


def main():
    formation_energies = json.loads(FORMATION_ENERGIES.read_text())
    rows = []
    n_ok = n_failed = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for d in sorted(STRUCTURES_ROOT.iterdir()):
            meta_path = d / "mp_metadata.json"
            if not meta_path.exists():
                continue
            meta = json.loads(meta_path.read_text())
            if meta.get("batch") != BATCH:
                continue
            is_metal = meta.get("is_metal")
            if is_metal is None or not (d / "COHPCAR.lobster").exists():
                continue
            try:
                per_delta = {}
                for delta_e in DELTA_ES:
                    per_delta[str(delta_e)] = ce.antibonding_population_near_frontier(
                        d, is_metal=is_metal, delta_e=delta_e, vasprun_path=d / "vasprun.xml"
                    )
                primary = per_delta[str(PRIMARY_DELTA_E)]
                mp_id = meta.get("mp_id")
                rows.append({
                    "compound_id": d.name,
                    "mp_id": mp_id,
                    "formula": meta.get("formula"),
                    "family": meta.get("family"),
                    "bond_type": meta.get("bond_type"),
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
                n_ok += 1
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                print(f"FAILED {d.name}: {type(exc).__name__}: {exc}")
                n_failed += 1

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV} (ok={n_ok}, failed={n_failed})")


if __name__ == "__main__":
    main()
