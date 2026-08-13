"""Download the compounds selected by select_campaign.py
(mp_dataset/campaign_selection.json) as POSCAR files, one directory per
compound under mp_dataset/structures/, matching the layout used by the
6-compound pilot (download_selected.py).
"""

import json
import os
from pathlib import Path

from mp_api.client import MPRester

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

SELECTION_PATH = Path(__file__).parent / "campaign_selection.json"
OUT_ROOT = Path(__file__).parent / "structures"


def dirname_for(entry: dict) -> str:
    formula = entry["formula"].replace(" ", "")
    return f"{entry['family']}_{formula}_{entry['mp_id']}"


def main():
    selection = json.loads(SELECTION_PATH.read_text())
    OUT_ROOT.mkdir(exist_ok=True)

    mp_ids = [e["mp_id"] for e in selection]
    with MPRester(API_KEY) as mpr:
        # batched lookup instead of one request per compound
        docs = mpr.materials.summary.search(
            material_ids=mp_ids,
            fields=["material_id", "structure"],
        )
    structures_by_id = {str(d.material_id): d.structure for d in docs}

    n_ok, n_missing = 0, 0
    for entry in selection:
        mp_id = entry["mp_id"]
        structure = structures_by_id.get(mp_id)
        if structure is None:
            print(f"MISSING structure for {mp_id} ({entry['formula']})")
            n_missing += 1
            continue

        compound_dir = OUT_ROOT / dirname_for(entry)
        compound_dir.mkdir(exist_ok=True)
        structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")
        with open(compound_dir / "mp_metadata.json", "w") as f:
            json.dump({**entry, "num_sites_structure": len(structure)}, f, indent=2)
        n_ok += 1

    print(f"\nDownloaded {n_ok}/{len(selection)} structures ({n_missing} missing) into {OUT_ROOT}")


if __name__ == "__main__":
    main()
