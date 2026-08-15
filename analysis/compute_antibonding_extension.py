"""Run the antibonding-population-near-frontier metric (validated on the 6
pilots, extended to the 186-compound campaign in compute_antibonding_all.py
-- mission #4) on the two extension_* batches instead: 14 hand-picked
compounds (magnetic/COD-sourced/missing-from-MP edge cases) plus 8
deliberately-off-equilibrium structures (3 high-pressure TiO2 polymorphs and
2 theoretical cold-compressed-graphite carbon allotropes, both sub-families
computed here at ambient pressure though their true stable regime is
elsewhere) -- a targeted diagnostic for the open question left in
REPORT_antibonding.md: does a structure known/expected to be
electronically or thermodynamically unstable show elevated antibonding
population near the frontier, as the Peierls/Jahn-Teller-adjacent reading
of this metric would predict?

is_metal comes from mp_metadata.json (already recorded per-compound at
download time, sourced from Materials Project for the MP-sourced entries --
same convention as compute_antibonding_all.py, not derived locally).
The 2 COD-sourced compounds (S4N2, S4N4) have is_metal=null in their
metadata -- no Materials Project entry exists for that composition to draw
on, and this project has an explicitly documented pitfall
(METRIC_DEFINITION_antibonding.md) showing the local LOBSTER-oriented
coarse-k-mesh gap estimate is unreliable (spurious small gaps for known
metals AlNi/BeCu). Rather than repeat that mistake, these two are skipped
here with an explicit reason, not silently guessed.

Writes analysis/antibonding_extension.json (per-compound, all delta_e
values) and analysis/antibonding_extension.csv (primary delta_e=1.0 next to
each compound's own mp_metadata.json fields).
"""

import json
import warnings
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
import sys  # noqa: E402
sys.path.insert(0, str(REPO_ROOT))
import cohp_extraction as ce  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT_JSON = Path(__file__).parent / "antibonding_extension.json"
OUT_CSV = Path(__file__).parent / "antibonding_extension.csv"

DELTA_ES = (0.5, 1.0, 2.0)
PRIMARY_DELTA_E = 1.0


def main():
    compound_dirs = sorted(STRUCTURES_ROOT.glob("extension_*"))
    print(f"Found {len(compound_dirs)} extension_* compound directories")

    results = {}
    rows = []
    n_ok = n_skipped = n_failed = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for d in compound_dirs:
            compound_id = d.name
            meta_path = d / "mp_metadata.json"
            meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
            is_metal = meta.get("is_metal")

            vasprun_path = d / "vasprun.xml"
            if not vasprun_path.exists():
                results[compound_id] = {"error": "not computed yet (no vasprun.xml)"}
                n_skipped += 1
                continue
            if is_metal is None:
                results[compound_id] = {"error": "is_metal unknown (no MP reference -- COD-sourced), skipped rather than derived locally"}
                n_skipped += 1
                continue

            try:
                per_delta = {}
                for delta_e in DELTA_ES:
                    r = ce.antibonding_population_near_frontier(
                        d, is_metal=is_metal, delta_e=delta_e, vasprun_path=vasprun_path
                    )
                    per_delta[str(delta_e)] = r
                results[compound_id] = {"is_metal": is_metal, "by_delta_e": per_delta, "error": None}
                n_ok += 1

                primary = per_delta[str(PRIMARY_DELTA_E)]
                rows.append({
                    "compound_id": compound_id,
                    "label": meta.get("label"),
                    "source": meta.get("source"),
                    "bond_type": meta.get("bond_type"),
                    "is_metal": is_metal,
                    "theoretical": meta.get("theoretical"),
                    "energy_above_hull_eV_per_atom": meta.get("energy_above_hull_eV_per_atom"),
                    "note": meta.get("note"),
                    "antibond_w_raw": primary["w_antibond_raw"],
                    "antibond_w_normalized": primary["w_antibond_normalized"],
                    "antibond_e_ref": primary["e_ref"],
                    **{f"antibond_w_raw_dE{de}": per_delta[str(de)]["w_antibond_raw"] for de in DELTA_ES},
                    **{f"antibond_w_normalized_dE{de}": per_delta[str(de)]["w_antibond_normalized"] for de in DELTA_ES},
                })
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                results[compound_id] = {"is_metal": is_metal, "error": f"{type(exc).__name__}: {exc}"}
                n_failed += 1
                print(f"FAILED {compound_id}: {type(exc).__name__}: {exc}")

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWrote {len(results)} entries to {OUT_JSON} (ok={n_ok}, skipped={n_skipped}, failed={n_failed})")

    df = pd.DataFrame(rows).sort_values("energy_above_hull_eV_per_atom", na_position="last")
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV}")
    print(df[["compound_id", "is_metal", "energy_above_hull_eV_per_atom", "antibond_w_normalized"]].to_string(index=False))


if __name__ == "__main__":
    main()
