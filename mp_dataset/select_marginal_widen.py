"""Widen the marginal-formation-energy campaign (marginal_candidates_merged.json,
75 compounds already downloaded/launched, see project memory) with more
candidates from the same marginal-FE regime, per explicit user request
(2026-08-17): "elargis la selection" using the same method as
select_marginal_formation_energy.py/select_marginal_ionic.py (same FE
band, same chemsys-diversity discipline) but with wider per-band/per-
chemsys targets, PLUS two new chemistry-targeted pools the first two
queries structurally cannot reach:

  - "organic": binary/ternary compounds built only from {C, H, N, O} --
    covalent molecular/framework solids (the periodic-wide query is
    binary-only and unconstrained in element choice, so it never
    targets this specific covalent-molecular chemistry on purpose).
    Ternary allowed here specifically because C/H/N/O chemistry is
    overwhelmingly ternary/quaternary in practice (few binary C-H-N-O
    combinations exist at all) -- deliberate, explicit exception to
    this project's usual binary-only rule, per this user request.
  - "molecular_pon": binary/ternary compounds among {P, N, O} only
    (P4O10, P3N5, N2O5, PON, ...) -- classic molecular-crystal
    chemistry (weak inter-molecular packing, strong intra-molecular
    covalent bonds), also ternary-allowed for the same reason.

All four pools share the same marginal-FE bands (FE_LO/FE_HI, both
signs) as the original campaign -- this is a widening of the same
regime, not a new one -- and the same non-magnetic/deprecated=False/
nsites<=40/all-elements-in-ELEMENT_REFERENCE constraints. Every pool
excludes mp_ids already in marginal_candidates_merged.json (the 75
already running) so this is a pure addition, no duplicate downloads/
jobs.

Selection-only: writes marginal_widen_candidates.json, does not
download or compute anything.
"""

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

from mp_api.client import MPRester

sys.path.insert(0, str(Path(__file__).parent.parent / "analysis"))
from compute_reaction_icohp_case1 import ELEMENT_REFERENCE  # noqa: E402

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()
HERE = Path(__file__).parent

MAX_SITES = 40
FE_LO, FE_HI = 0.02, 0.15
NON_MAGNETIC_WINDOW = (-0.01, 0.01)
ALLOWED_ELEMENTS = frozenset(ELEMENT_REFERENCE.keys())

ALREADY_SELECTED = {r["mp_id"] for r in json.loads((HERE / "marginal_candidates_merged.json").read_text())}

FIELDS = [
    "material_id", "formula_pretty", "chemsys", "elements",
    "energy_above_hull", "formation_energy_per_atom", "theoretical",
    "nsites", "nelements", "symmetry", "is_metal", "band_gap",
    "total_magnetization_normalized_formula_units",
]


def to_record(pool, band, d):
    return {
        "pool": pool,
        "band": band,
        "mp_id": str(d.material_id),
        "formula": d.formula_pretty,
        "chemsys": d.chemsys,
        "formation_energy_per_atom_eV": d.formation_energy_per_atom,
        "energy_above_hull_eV_per_atom": d.energy_above_hull,
        "theoretical": d.theoretical,
        "nsites": d.nsites,
        "spacegroup": d.symmetry.symbol if d.symmetry else None,
        "is_metal": d.is_metal,
        "band_gap_eV": d.band_gap,
        "total_magnetization_per_fu": d.total_magnetization_normalized_formula_units,
    }


def diversify(records, n_target, max_per_chemsys):
    by_chemsys = defaultdict(list)
    for r in records:
        by_chemsys[r["chemsys"]].append(r)
    kept = []
    for chemsys, group in by_chemsys.items():
        group.sort(key=lambda r: r["theoretical"])  # experimental first
        kept.extend(group[:max_per_chemsys])
    kept.sort(key=lambda r: (r["chemsys"], r["formula"]))
    if len(kept) > n_target:
        step = len(kept) / n_target
        kept = [kept[int(i * step)] for i in range(n_target)]
    return kept


def fetch_periodic_wide(mpr, n_target_per_band, max_per_chemsys):
    """Same as select_marginal_formation_energy.py but wider targets,
    excluding mp_ids already selected in the first pass."""
    pos = mpr.materials.summary.search(
        num_elements=2, num_sites=(1, MAX_SITES), deprecated=False,
        formation_energy_per_atom=(FE_LO, FE_HI),
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    neg = mpr.materials.summary.search(
        num_elements=2, num_sites=(1, MAX_SITES), deprecated=False,
        formation_energy_per_atom=(-FE_HI, -FE_LO),
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    pos = [d for d in pos if ALLOWED_ELEMENTS.issuperset({str(e) for e in d.elements})
           and str(d.material_id) not in ALREADY_SELECTED]
    neg = [d for d in neg if ALLOWED_ELEMENTS.issuperset({str(e) for e in d.elements})
           and str(d.material_id) not in ALREADY_SELECTED]
    pos_records = diversify([to_record("periodic_wide", "positive", d) for d in pos], n_target_per_band, max_per_chemsys)
    neg_records = diversify([to_record("periodic_wide", "negative", d) for d in neg], n_target_per_band, max_per_chemsys)
    print(f"  periodic_wide: raw pos={len(pos)} neg={len(neg)} -> kept {len(pos_records)}+{len(neg_records)}")
    return pos_records + neg_records


def fetch_ionic(mpr, n_target, max_per_chemsys):
    """Same chemistry-targeted alkali/alkaline-earth x
    {halogen,O,N,P,S} sweep as select_marginal_ionic.py, wider target,
    excluding already-selected mp_ids."""
    ALKALI = ["Li", "Na", "K", "Rb", "Cs"]
    ALKALINE_EARTH = ["Be", "Mg", "Ca", "Sr", "Ba"]
    CATIONS = ALKALI + ALKALINE_EARTH
    ANIONS = ["F", "Cl", "Br", "I", "O", "N", "P", "S"]

    all_records = []
    for cation in CATIONS:
        for anion in ANIONS:
            if cation not in ALLOWED_ELEMENTS or anion not in ALLOWED_ELEMENTS:
                continue
            docs = mpr.materials.summary.search(
                chemsys=f"{cation}-{anion}", num_elements=2, num_sites=(1, MAX_SITES),
                deprecated=False, total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
                fields=FIELDS,
            )
            for d in docs:
                if str(d.material_id) in ALREADY_SELECTED:
                    continue
                fe = d.formation_energy_per_atom
                if fe is None:
                    continue
                if FE_LO <= fe <= FE_HI:
                    band = "positive"
                elif -FE_HI <= fe <= -FE_LO:
                    band = "negative"
                else:
                    continue
                all_records.append(to_record("ionic_targeted", band, d))
    print(f"  ionic_targeted: raw in-band={len(all_records)}")
    return diversify(all_records, n_target, max_per_chemsys)


def fetch_by_element_universe(mpr, pool_name, universe, n_target, max_per_chemsys, num_elements=(2, 3)):
    """Binary/ternary compounds drawn ONLY from `universe`, same
    marginal-FE bands, both signs.

    `elements=` in mp_api is an AND-filter (compound must contain
    EVERY listed element), not a subset/OR filter -- passing the whole
    universe there silently requires all of C+H+N+O (or P+N+O)
    simultaneously, which no binary/ternary compound can satisfy
    (returns 0 candidates, the bug this function was rewritten to
    avoid). Instead enumerate every binary/ternary chemsys drawable
    FROM the universe and OR them together via `chemsys=[...]`, same
    approach as select_marginal_ionic.py's per-pair sweep but batched
    into one query per sign via the list form."""
    from itertools import combinations
    universe = sorted(e for e in universe if e in ALLOWED_ELEMENTS)
    sizes = num_elements if isinstance(num_elements, (list, tuple)) else [num_elements]
    chemsyses = ["-".join(combo) for n in sizes for combo in combinations(universe, n)]

    pos = mpr.materials.summary.search(
        chemsys=chemsyses, num_sites=(1, MAX_SITES), deprecated=False,
        formation_energy_per_atom=(FE_LO, FE_HI),
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    neg = mpr.materials.summary.search(
        chemsys=chemsyses, num_sites=(1, MAX_SITES), deprecated=False,
        formation_energy_per_atom=(-FE_HI, -FE_LO),
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    pos = [d for d in pos if str(d.material_id) not in ALREADY_SELECTED]
    neg = [d for d in neg if str(d.material_id) not in ALREADY_SELECTED]
    print(f"  {pool_name}: chemsyses={chemsyses} raw pos={len(pos)} neg={len(neg)}")
    pos_records = diversify([to_record(pool_name, "positive", d) for d in pos], n_target // 2, max_per_chemsys)
    neg_records = diversify([to_record(pool_name, "negative", d) for d in neg], n_target - n_target // 2, max_per_chemsys)
    return pos_records + neg_records


def main():
    with MPRester(API_KEY) as mpr:
        print("Fetching widened periodic-wide FE pool...")
        periodic = fetch_periodic_wide(mpr, n_target_per_band=50, max_per_chemsys=3)

        print("Fetching widened ionic-targeted pool...")
        ionic = fetch_ionic(mpr, n_target=50, max_per_chemsys=4)

        print("Fetching organic C/H/N/O pool...")
        organic = fetch_by_element_universe(mpr, "organic", ["C", "H", "N", "O"], n_target=20, max_per_chemsys=3)

        print("Fetching molecular P/N/O pool...")
        molecular_pon = fetch_by_element_universe(mpr, "molecular_pon", ["P", "N", "O"], n_target=10, max_per_chemsys=3)

    combined = periodic + ionic + organic + molecular_pon

    # Global dedup + re-apply chemsys diversity across the COMBINED
    # widening batch (same discipline as merge_marginal_candidates.py).
    seen_ids = set()
    deduped = []
    for r in combined:
        if r["mp_id"] in seen_ids:
            continue
        seen_ids.add(r["mp_id"])
        deduped.append(r)

    deduped.sort(key=lambda r: (bool(r["theoretical"]), abs(r["formation_energy_per_atom_eV"])))
    by_chemsys = defaultdict(list)
    kept = []
    for r in deduped:
        cap = 4 if r["pool"] in ("organic", "molecular_pon") else 3
        if len(by_chemsys[r["chemsys"]]) >= cap:
            continue
        by_chemsys[r["chemsys"]].append(r)
        kept.append(r)

    kept.sort(key=lambda r: (r["pool"], r["chemsys"], r["formula"]))

    out_path = HERE / "marginal_widen_candidates.json"
    out_path.write_text(json.dumps(kept, indent=2))

    n_by_pool = defaultdict(int)
    n_metal = n_nonmetal = 0
    for r in kept:
        n_by_pool[r["pool"]] += 1
        if r["is_metal"]:
            n_metal += 1
        else:
            n_nonmetal += 1

    print(f"\n{len(deduped)} unique new entries -> {len(kept)} after chemsys diversity")
    print(f"By pool: {dict(n_by_pool)}")
    print(f"metallic={n_metal}, non-metallic={n_nonmetal}")
    print(f"Wrote {len(kept)} entries to {out_path}\n")

    print(f"{'pool':<16}{'band':<10}{'formula':<14}{'mp_id':<12}{'chemsys':<10}{'FE(eV/at)':>11}  "
          f"sg{'':<10}nsites  is_metal  theo")
    for r in kept:
        sg = r["spacegroup"] or "?"
        print(f"{r['pool']:<16}{r['band']:<10}{r['formula']:<14}{r['mp_id']:<12}{r['chemsys']:<10}"
              f"{r['formation_energy_per_atom_eV']:>11.4f}  {sg:<12}{r['nsites']:<8}"
              f"{str(r['is_metal']):<10}{r['theoretical']}")


if __name__ == "__main__":
    main()
