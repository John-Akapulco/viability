"""Fetch formation_energy_per_atom for the 186 compounds already in
mp_dataset/structures/*/mp_metadata.json (not present there -- that file
only carries energy_above_hull_eV_per_atom). Single batched MP query, same
API key convention as fetch_candidates.py/download_selected.py.

Writes mp_dataset/formation_energies.json ({mp_id: formation_energy_per_atom})
as a separate file -- does not rewrite the existing mp_metadata.json files.
"""

import json
import os
from pathlib import Path

from mp_api.client import MPRester

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

STRUCTURES_ROOT = Path(__file__).parent / "structures"
OUT_PATH = Path(__file__).parent / "formation_energies.json"


def collect_mp_ids() -> list[str]:
    mp_ids = set()
    for compound_dir in STRUCTURES_ROOT.iterdir():
        meta_path = compound_dir / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        mp_id = meta.get("mp_id") or meta.get("material_id")
        if mp_id:
            mp_ids.add(mp_id)
    return sorted(mp_ids)


def main():
    mp_ids = collect_mp_ids()
    print(f"Fetching formation_energy_per_atom for {len(mp_ids)} compounds...")

    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=mp_ids,
            fields=["material_id", "formation_energy_per_atom"],
        )

    result = {str(d.material_id): d.formation_energy_per_atom for d in docs}
    missing = sorted(set(mp_ids) - set(result))
    if missing:
        print(f"WARNING: {len(missing)} mp_id(s) returned no data: {missing}")

    OUT_PATH.write_text(json.dumps(result, indent=2))
    print(f"Wrote {len(result)} entries to {OUT_PATH}")


if __name__ == "__main__":
    main()
