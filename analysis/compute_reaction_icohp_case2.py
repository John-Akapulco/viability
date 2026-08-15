"""Case 2 (polymorph comparison) of reaction_icohp.py, run across every
group of same-reduced-formula compounds already in mp_dataset/structures/
-- the batch-scale extension flagged as not-yet-done in
METRIC_DEFINITION_reaction_icohp.md section 6 (only n=1-5 worked examples,
carbon allotropes and TiO2 polymorphs, existed before this).

Unlike case 1, no elemental references or reaction balancing are needed
(reaction_icohp.compare_polymorphs, unmodified) -- every group of >=2
compounds sharing a pymatgen reduced_formula (from mp_metadata.json's
"formula" field, same grouping convention as compute_reaction_icohp_case1.py)
is a directly comparable polymorph set.

Only 8 such groups exist in the current 260-compound structures/ tree
(21 compound-dirs total) -- this project's dataset was built for breadth
across chemical systems (select_campaign.py caps 2 entries per chemsys),
not for polymorph density, so this is inherently a small-n exploratory
pass, not a well-powered statistical test.

For each group, reports whether the most-bonding polymorph (min
icohp_per_atom) is also the most-stable one (min energy_above_hull) --
the same "does more ICOHP bonding track more thermodynamic stability
within a fixed composition" question the worked examples in
METRIC_DEFINITION_reaction_icohp.md already answered "no" for twice
(carbon allotropes, TiO2 polymorphs).

Writes analysis/reaction_icohp_case2.json and analysis/reaction_icohp_case2.csv.
"""

import json
import warnings
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

import sys

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
import reaction_icohp as ri  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT_JSON = Path(__file__).parent / "reaction_icohp_case2.json"
OUT_CSV = Path(__file__).parent / "reaction_icohp_case2.csv"


def is_computed(compound_dir: Path) -> bool:
    return (compound_dir / "ICOHPLIST.lobster").exists() and (compound_dir / "CONTCAR").exists()


def find_polymorph_groups() -> dict:
    groups = defaultdict(list)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
            meta_path = compound_dir / "mp_metadata.json"
            if not meta_path.exists() or not is_computed(compound_dir):
                continue
            meta = json.loads(meta_path.read_text())
            formula = meta.get("formula")
            if not formula:
                continue
            try:
                rf = Composition(formula).reduced_formula
            except Exception:
                continue
            groups[rf].append(compound_dir)
    return {rf: dirs for rf, dirs in groups.items() if len(dirs) >= 2}


def main():
    groups = find_polymorph_groups()
    print(f"Found {len(groups)} polymorph groups ({sum(len(v) for v in groups.values())} compounds)")

    results = {}
    rows = []
    n_agree = n_total_groups = 0
    for rf, dirs in sorted(groups.items()):
        try:
            cmp = ri.compare_polymorphs(dirs)
        except Exception as exc:  # noqa: BLE001 - batch must not die on one bad group
            results[rf] = {"error": f"{type(exc).__name__}: {exc}"}
            print(f"FAILED {rf}: {type(exc).__name__}: {exc}")
            continue

        for poly in cmp["polymorphs"]:
            d = Path(poly["compound_dir"])
            meta = json.loads((d / "mp_metadata.json").read_text())
            poly["compound_id"] = d.name
            poly["mp_id"] = meta.get("mp_id")
            poly["family"] = meta.get("family")
            poly["theoretical"] = meta.get("theoretical")
            poly["energy_above_hull_eV_per_atom"] = meta.get("energy_above_hull_eV_per_atom")

        with_eah = [p for p in cmp["polymorphs"] if p["energy_above_hull_eV_per_atom"] is not None]
        most_bonding = cmp["polymorphs"][0]["compound_id"]  # already sorted by icohp_per_atom
        most_stable = None
        agrees = None
        if len(with_eah) >= 2:
            most_stable = min(with_eah, key=lambda p: p["energy_above_hull_eV_per_atom"])["compound_id"]
            agrees = most_bonding == most_stable
            n_total_groups += 1
            n_agree += int(agrees)

        cmp["most_bonding_compound_id"] = most_bonding
        cmp["most_stable_compound_id"] = most_stable
        cmp["most_bonding_matches_most_stable"] = agrees
        results[rf] = cmp
        rows.extend(cmp["polymorphs"])

        print(f"{rf} (n={len(dirs)}): most-bonding={most_bonding}, most-stable={most_stable}, "
              f"agree={agrees}")

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(results)} groups to {OUT_JSON}, {len(df)} rows to {OUT_CSV}")
    print(f"Most-bonding == most-stable in {n_agree}/{n_total_groups} groups with >=2 EAH-labeled members")


if __name__ == "__main__":
    main()
