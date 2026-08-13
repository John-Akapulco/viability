"""Select a ~180-compound campaign from Materials Project:

  - 60 experimental, thermodynamically stable (on hull)
  - 60 experimental, metastable (0 < E_above_hull <= 100 meV/atom)
  - 60 theoretical, metastable (0 < E_above_hull <= 200 meV/atom)

Constraints (all three families): unary or binary chemical systems,
non-magnetic (Ordering.NM), no f-block elements (lanthanides/actinides),
<=20 sites/cell (kept from the pilot run, for VASP+LOBSTER tractability
at this scale), deprecated=False.

Diversity: at most 2 entries per reduced chemical system (chemsys), so 60
slots are not dominated by many polymorphs of the same 1-2 elements;
random sampling (fixed seed) over the filtered pool for a broad spread
across the periodic table.

Writes mp_dataset/campaign_selection.json (list of {family, mp_id,
formula, chemsys, energy_above_hull, theoretical, nsites}) -- does not
download structures (see download_campaign.py).
"""

import json
import os
import random
from collections import defaultdict
from pathlib import Path

from mp_api.client import MPRester

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

MAX_SITES = 20
N_PER_FAMILY = 60
MAX_PER_CHEMSYS = 2
SEED = 20260813

LANTHANIDES = [
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu",
]
ACTINIDES = [
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr",
]
EXCLUDE_F_ELEMENTS = LANTHANIDES + ACTINIDES

# Already used by the 6-compound pilot -- excluded so the campaign adds
# genuinely new compounds rather than reprocessing them.
ALREADY_USED_MP_IDS = {"mp-22862", "mp-149", "mp-1487", "mp-23259", "mp-169", "mp-2323"}

FIELDS = [
    "material_id",
    "formula_pretty",
    "chemsys",
    "elements",
    "energy_above_hull",
    "theoretical",
    "nsites",
    "nelements",
    "symmetry",
    "total_magnetization_normalized_formula_units",
]

_EXCLUDE_F_SET = frozenset(EXCLUDE_F_ELEMENTS)

# The categorical `magnetic_ordering` field is "Unknown" for the vast
# majority of entries (not curated, not necessarily magnetic) -- filtering
# on it drops ~98% of genuinely non-magnetic compounds. Use the actual
# computed total magnetization instead: near-zero net moment per formula
# unit is the direct physical definition of "non-magnetic".
NON_MAGNETIC_WINDOW = (-0.01, 0.01)  # mu_B / formula unit


def fetch_pool(mpr: MPRester, *, theoretical: bool, e_above_hull_window, is_stable=None):
    # exclude_elements hits a server-side "String should have at most 60
    # characters" validation error once the excluded-element list is this
    # long (all lanthanides+actinides serialize to a ~90-char CSV) -- filter
    # f-block elements client-side instead of relying on the API parameter.
    kwargs = dict(
        num_elements=(1, 2),
        num_sites=(1, MAX_SITES),
        theoretical=theoretical,
        deprecated=False,
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    if is_stable is not None:
        kwargs["is_stable"] = is_stable
    else:
        kwargs["energy_above_hull"] = e_above_hull_window
    docs = mpr.materials.summary.search(**kwargs)
    return [d for d in docs if not (_EXCLUDE_F_SET & {str(e) for e in d.elements})]


def diverse_sample(docs, n, rng, exclude_ids=frozenset()):
    by_chemsys = defaultdict(list)
    for d in docs:
        if d.material_id in exclude_ids or str(d.material_id) in exclude_ids:
            continue
        by_chemsys[d.chemsys].append(d)

    pool = []
    for chemsys, items in by_chemsys.items():
        items = sorted(items, key=lambda d: (d.nsites, d.energy_above_hull))
        pool.extend(items[:MAX_PER_CHEMSYS])

    rng.shuffle(pool)
    return pool[:n]


def to_record(family, d):
    return {
        "family": family,
        "mp_id": str(d.material_id),
        "formula": d.formula_pretty,
        "chemsys": d.chemsys,
        "energy_above_hull_eV_per_atom": d.energy_above_hull,
        "theoretical": d.theoretical,
        "nsites": d.nsites,
        "nelements": d.nelements,
        "spacegroup": d.symmetry.symbol if d.symmetry else None,
        "total_magnetization_per_fu": d.total_magnetization_normalized_formula_units,
    }


def main():
    rng = random.Random(SEED)
    selection = []

    with MPRester(API_KEY) as mpr:
        print("Fetching: experimental, stable (on hull)...")
        pool_stable = fetch_pool(mpr, theoretical=False, e_above_hull_window=None, is_stable=True)
        sample = diverse_sample(pool_stable, N_PER_FAMILY, rng, ALREADY_USED_MP_IDS)
        selection += [to_record("exp_stable", d) for d in sample]
        print(f"  pool={len(pool_stable)}  selected={len(sample)}")

        print("Fetching: experimental, metastable (0-100 meV/atom)...")
        pool_meta_exp = fetch_pool(mpr, theoretical=False, e_above_hull_window=(0.0005, 0.100))
        sample = diverse_sample(pool_meta_exp, N_PER_FAMILY, rng, ALREADY_USED_MP_IDS)
        selection += [to_record("exp_metastable", d) for d in sample]
        print(f"  pool={len(pool_meta_exp)}  selected={len(sample)}")

        print("Fetching: theoretical, metastable (0-200 meV/atom)...")
        pool_meta_theo = fetch_pool(mpr, theoretical=True, e_above_hull_window=(0.0005, 0.200))
        sample = diverse_sample(pool_meta_theo, N_PER_FAMILY, rng, ALREADY_USED_MP_IDS)
        selection += [to_record("theo_metastable", d) for d in sample]
        print(f"  pool={len(pool_meta_theo)}  selected={len(sample)}")

    out_path = Path(__file__).parent / "campaign_selection.json"
    with open(out_path, "w") as f:
        json.dump(selection, f, indent=2)

    print(f"\nWrote {len(selection)} entries to {out_path}")
    for fam in ("exp_stable", "exp_metastable", "theo_metastable"):
        n = sum(1 for s in selection if s["family"] == fam)
        print(f"  {fam}: {n}")


if __name__ == "__main__":
    main()
