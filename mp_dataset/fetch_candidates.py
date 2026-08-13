"""Query the Materials Project API for candidate compounds, split into
2 families (on-hull / metastable-experimental) x 3 bonding types
(ionic / covalent / metallic-alloy). Prints a shortlist per bucket;
does not download full structures yet (see download_selected.py).
"""

from __future__ import annotations

import os
from collections import defaultdict

from mp_api.client import MPRester

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

MAX_SITES = 20
METASTABLE_WINDOW = (0.0005, 0.100)  # eV/atom, strictly above hull up to 100 meV/atom
HULL_WINDOW = (0.0, 0.0005)

ALKALI = {"Li", "Na", "K", "Rb", "Cs"}
ALKALINE_EARTH = {"Be", "Mg", "Ca", "Sr", "Ba"}
HALOGENS = {"F", "Cl", "Br", "I"}
IONIC_ANIONS = HALOGENS | {"O"}
COVALENT_FORMERS = {"B", "C", "Si", "Ge", "N", "P", "As", "Al", "Ga", "In", "Sn", "Pb", "S", "Se", "Te"}
ANION_LIKE_FOR_ALLOY_EXCLUSION = IONIC_ANIONS | {"N", "S", "Se", "Te", "H", "C", "P", "As", "B"}

FIELDS = [
    "material_id",
    "formula_pretty",
    "energy_above_hull",
    "is_metal",
    "band_gap",
    "elements",
    "nelements",
    "nsites",
    "theoretical",
    "symmetry",
]


def classify(elements: set, is_metal: bool) -> str | None:
    if is_metal:
        if len(elements) < 2:
            return None
        if elements & ANION_LIKE_FOR_ALLOY_EXCLUSION:
            return None
        return "metallic"
    if (elements & (ALKALI | ALKALINE_EARTH)) and (elements & IONIC_ANIONS):
        return "ionic"
    if (elements & COVALENT_FORMERS) and not (elements & (ALKALI | ALKALINE_EARTH)):
        return "covalent"
    return None


def fetch_bucket(mpr: MPRester, e_above_hull_window, label: str):
    docs = mpr.materials.summary.search(
        energy_above_hull=e_above_hull_window,
        num_sites=(1, MAX_SITES),
        theoretical=False,
        deprecated=False,
        fields=FIELDS,
    )
    buckets: dict[str, list] = defaultdict(list)
    for d in docs:
        els = {str(e) for e in d.elements}
        kind = classify(els, d.is_metal)
        if kind is None:
            continue
        buckets[kind].append(d)

    for kind in ("ionic", "covalent", "metallic"):
        items = sorted(buckets[kind], key=lambda d: (d.nsites, d.nelements, d.energy_above_hull))
        print(f"\n=== {label} / {kind}  ({len(items)} candidates) ===")
        for d in items[:8]:
            print(
                f"  {d.material_id}  {d.formula_pretty:<15} "
                f"E_hull={d.energy_above_hull:.4f} eV/at  nsites={d.nsites}  "
                f"sg={d.symmetry.symbol if d.symmetry else '?'}  "
                f"is_metal={d.is_metal}  gap={d.band_gap}"
            )
    return buckets


def main():
    with MPRester(API_KEY) as mpr:
        print("################ HULL FAMILY (E_above_hull ~ 0) ################")
        fetch_bucket(mpr, HULL_WINDOW, "hull")
        print("\n\n################ METASTABLE FAMILY (0-100 meV/atom, experimental) ################")
        fetch_bucket(mpr, METASTABLE_WINDOW, "metastable")


if __name__ == "__main__":
    main()
