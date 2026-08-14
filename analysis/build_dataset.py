"""Join percolation_path.py's output with each compound's Materials Project
metadata into a single flat table for statistical analysis.

Inputs:
  - mp_dataset/results.json           (percolation_path.py output, all compounds)
  - mp_dataset/structures/*/mp_metadata.json  (per-compound MP provenance)
  - mp_dataset/structures/*/vasp_time.txt, lobster_time.txt (compute cost, if present)

bond_type (ionic/covalent/metallic) reuses classify() from
mp_dataset/fetch_candidates.py UNCHANGED -- the 6 pilot compounds already
carry it in mp_metadata.json (selected via that function); for the 180
campaign compounds it's derived here from a single batched MP query
(is_metal + elements), same logic, not reimplemented.

Output: analysis/percolation_vs_hull.csv
"""

import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "mp_dataset"))
from fetch_candidates import classify  # noqa: E402  (reused unmodified)

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
RESULTS_JSON = REPO_ROOT / "mp_dataset" / "results.json"
OUT_CSV = Path(__file__).parent / "percolation_vs_hull.csv"

# The 6-compound pilot (fetch_candidates.py) used a 2-way family label
# (hull/metastable, both experimental by construction); the 180-compound
# campaign (select_campaign.py) uses a 3-way label. Map the pilot's onto
# the campaign's so `family` groups consistently across all 186 compounds.
_PILOT_FAMILY_MAP = {"hull": "exp_stable", "metastable": "exp_metastable"}


def load_metadata(compound_dir: Path) -> dict:
    meta = json.loads((compound_dir / "mp_metadata.json").read_text())
    mp_id = meta.get("mp_id") or meta.get("material_id")
    e_hull = meta.get("energy_above_hull_eV_per_atom")
    nsites = meta.get("nsites") or meta.get("num_sites")
    family = _PILOT_FAMILY_MAP.get(meta.get("family"), meta.get("family"))
    return {
        "mp_id": mp_id,
        "formula": meta.get("formula"),
        "family": family,
        "bond_type": meta.get("bond_type"),  # present for the 6 pilots, None for campaign compounds
        "energy_above_hull_eV_at": e_hull,
        "theoretical": meta.get("theoretical"),
        "n_sites": nsites,
        "spacegroup_symbol": meta.get("spacegroup"),
    }


def parse_time_file(path: Path) -> float | None:
    """Parse `time -p` output ('real 123.45\\nuser ...\\nsys ...') -> real seconds."""
    if not path.exists():
        return None
    m = re.search(r"^real\s+([\d.]+)", path.read_text(), re.MULTILINE)
    return float(m.group(1)) if m else None


def parse_vasp_time_from_log(vasp_log: Path) -> float | None:
    """submit_array.sh wraps `srun vasp_std` as
    `{ time -p srun ... > vasp.log 2>&1; } 2> vasp_time.txt`, intending the
    `time -p` report to land in vasp_time.txt. In practice (observed across
    every compound checked) it lands at the tail of vasp.log instead --
    likely an artifact of how srun's remote-task I/O forwarding interacts
    with the shell's own fd bookkeeping under SLURM, not reproducible in a
    plain interactive shell. vasp_time.txt is consistently empty; parse the
    real, always-present location instead of chasing the redirect bug
    further (cosmetic, doesn't affect any actual compute results)."""
    if not vasp_log.exists():
        return None
    m = re.search(r"^real\s+([\d.]+)", vasp_log.read_text(errors="replace"), re.MULTILINE)
    return float(m.group(1)) if m else None


def fetch_bond_types_and_sg_numbers(mp_ids: list[str]) -> dict[str, dict]:
    """Single batched MP query for compounds missing bond_type (campaign
    set) and for sg_number (not stored anywhere yet, pilots included)."""
    import os as _os
    from mp_api.client import MPRester

    api_key = open(_os.path.expanduser("~/.mp_api_key")).read().strip()
    out = {}
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=mp_ids,
            fields=["material_id", "is_metal", "elements", "symmetry"],
        )
    for d in docs:
        els = {str(e) for e in d.elements}
        out[str(d.material_id)] = {
            "bond_type": classify(els, d.is_metal),
            "sg_number": d.symmetry.number if d.symmetry else None,
        }
    return out


def main():
    results = json.loads(RESULTS_JSON.read_text())
    results_by_id = {r["compound_id"]: r for r in results}

    rows = []
    all_mp_ids = []
    for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
        meta_path = compound_dir / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = load_metadata(compound_dir)
        if meta["mp_id"]:
            all_mp_ids.append(meta["mp_id"])
        rows.append({"compound_id": compound_dir.name, "_dir": compound_dir, **meta})

    lookup = fetch_bond_types_and_sg_numbers(sorted(set(m["mp_id"] for m in rows if m["mp_id"])))

    n_missing_result = 0
    final_rows = []
    for r in rows:
        compound_id = r["compound_id"]
        result = results_by_id.get(compound_id)
        if result is None:
            n_missing_result += 1
            continue  # not yet processed by percolation_path.py (job still running/pending/failed)
        if result.get("error"):
            continue  # percolation_path.py itself flagged this compound as failed

        extra = lookup.get(r["mp_id"], {})
        bond_type = r["bond_type"] or extra.get("bond_type")

        icohp = result.get("metrics", {}).get("icohp", {})
        icobi = result.get("metrics", {}).get("icobi", {})
        icohp_agg = icohp.get("aggregates", {})
        icobi_agg = icobi.get("aggregates", {})

        weight_min = icohp.get("min_weight")
        icohp_min = icohp_agg.get("min")
        weight_min_norm = (
            weight_min / abs(icohp_min) if (weight_min is not None and icohp_min) else None
        )

        final_rows.append({
            "compound_id": compound_id,
            "mp_id": r["mp_id"],
            "formula": r["formula"],
            "bond_type": bond_type,
            "family": r["family"],
            "theoretical": r["theoretical"],
            "energy_above_hull_eV_at": r["energy_above_hull_eV_at"],
            "n_sites": r["n_sites"],
            "sg_number": extra.get("sg_number"),
            "spacegroup_symbol": r["spacegroup_symbol"],
            "icohp_percolation_weight_min": weight_min,
            "icohp_percolation_weight_min_normalized": weight_min_norm,
            "icohp_percolation_direction": icohp.get("min_direction"),
            "icohp_sum": icohp_agg.get("sum"),
            "icohp_mean": icohp_agg.get("mean"),
            "icohp_min": icohp_min,
            "icohp_max": icohp_agg.get("max"),
            "icobi_percolation_weight_min": icobi.get("min_weight"),
            "icobi_sum": icobi_agg.get("sum"),
            "icobi_mean": icobi_agg.get("mean"),
            "icobi_min": icobi_agg.get("min"),
            "icobi_max": icobi_agg.get("max"),
            "vasp_wall_time_s": parse_vasp_time_from_log(r["_dir"] / "vasp.log"),
            "lobster_wall_time_s": parse_time_file(r["_dir"] / "lobster_time.txt"),
        })

    df = pd.DataFrame(final_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV} ({n_missing_result} compounds skipped: no percolation_path.py result yet)")
    print(df["family"].value_counts())
    print(df["bond_type"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
