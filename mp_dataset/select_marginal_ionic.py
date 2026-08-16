"""Targeted follow-up to select_marginal_formation_energy.py: that
script's single [FE_LO, FE_HI] query, run across the whole periodic
table, returned essentially zero ionic compounds (46/50 metallic, 4
non-metallic none of them ionic) -- classic ionic bonding (alkali/
alkaline-earth x halogen/O) is stabilizing enough that even MP's worst
metastable polymorphs of these systems tend to keep a strongly negative
formation_energy_per_atom, well outside any near-zero window. Ionic
candidates near the marginal-FE boundary won't show up from an
unconstrained periodic-table-wide query; they need a chemistry-targeted
one instead.

Cations: alkali (Li-Cs) + alkaline-earth (Be-Ba). Anions: not just
{halogen, O} (fetch_candidates.classify()'s IONIC_ANIONS, the
"traditionally ionic" set) but also {N, P, S} -- the exact set
classify() itself is documented (report Sec. 5, Sec 2.1) to NOT
recognize as ionic despite the cation/anion character being the same;
nitrides/sulfides/phosphides of these cations are less electronegative-
driven and empirically more likely to land near a marginal formation
energy than a halide/oxide of the same cation.

Deliberately NOT pre-constrained to a narrow [FE_LO, FE_HI] window up
front like the periodic-table-wide script -- queries this specific
chemistry's FULL formation-energy distribution first (still
non-magnetic, deprecated=False, binary, nsites<=40, both elements in
ELEMENT_REFERENCE) and reports it, so the actual achievable range is
seen before deciding a cutoff, rather than guessing one and getting an
empty result the way a halide/oxide-only version would.

Selection-only, same convention as every prior script in this family:
writes mp_dataset/marginal_ionic_candidates.json, does not download or
compute anything.
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

MAX_SITES = 40
MAX_PER_CHEMSYS = 3
N_TARGET = 30

ALKALI = ["Li", "Na", "K", "Rb", "Cs"]
ALKALINE_EARTH = ["Be", "Mg", "Ca", "Sr", "Ba"]
CATIONS = ALKALI + ALKALINE_EARTH
IONIC_ANIONS_TRADITIONAL = ["F", "Cl", "Br", "I", "O"]
IONIC_ANIONS_EXTENDED = ["N", "P", "S"]
ANIONS = IONIC_ANIONS_TRADITIONAL + IONIC_ANIONS_EXTENDED

ALLOWED_ELEMENTS = frozenset(ELEMENT_REFERENCE.keys())
NON_MAGNETIC_WINDOW = (-0.01, 0.01)

FIELDS = [
    "material_id", "formula_pretty", "chemsys", "elements",
    "energy_above_hull", "formation_energy_per_atom", "theoretical",
    "nsites", "nelements", "symmetry", "is_metal", "band_gap",
    "total_magnetization_normalized_formula_units",
]


def fetch_system(mpr: MPRester, cation: str, anion: str):
    if cation not in ALLOWED_ELEMENTS or anion not in ALLOWED_ELEMENTS:
        return []
    docs = mpr.materials.summary.search(
        chemsys=f"{cation}-{anion}",
        num_elements=2,
        num_sites=(1, MAX_SITES),
        deprecated=False,
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    return docs


def to_record(d, anion_class):
    return {
        "anion_class": anion_class,
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
    }


def main():
    all_records = []
    with MPRester(API_KEY) as mpr:
        for anion_class, anions in (("traditional", IONIC_ANIONS_TRADITIONAL), ("extended", IONIC_ANIONS_EXTENDED)):
            for cation in CATIONS:
                for anion in anions:
                    docs = fetch_system(mpr, cation, anion)
                    all_records.extend(to_record(d, anion_class) for d in docs)

    print(f"Total entries across all alkali/alkaline-earth x {{halogen,O,N,P,S}} systems: {len(all_records)}")

    fe_vals = [r["formation_energy_per_atom_eV"] for r in all_records if r["formation_energy_per_atom_eV"] is not None]
    print(f"formation_energy_per_atom range: [{min(fe_vals):.4f}, {max(fe_vals):.4f}] eV/atom")

    for label, anions in (("traditional (halogen/O)", IONIC_ANIONS_TRADITIONAL), ("extended (N/P/S)", IONIC_ANIONS_EXTENDED)):
        sub = [r for r in all_records if r["anion_class"] == ("traditional" if "traditional" in label else "extended")]
        sub_fe = [r["formation_energy_per_atom_eV"] for r in sub if r["formation_energy_per_atom_eV"] is not None]
        if sub_fe:
            print(f"  {label}: n={len(sub)}, FE range [{min(sub_fe):.4f}, {max(sub_fe):.4f}], "
                  f"median {sorted(sub_fe)[len(sub_fe)//2]:.4f}")

    # Rank by |FE| ascending (closest to marginal), diversify by chemsys.
    ranked = sorted(all_records, key=lambda r: abs(r["formation_energy_per_atom_eV"]) if r["formation_energy_per_atom_eV"] is not None else 1e9)
    by_chemsys = defaultdict(list)
    kept = []
    for r in ranked:
        if len(by_chemsys[r["chemsys"]]) >= MAX_PER_CHEMSYS:
            continue
        by_chemsys[r["chemsys"]].append(r)
        kept.append(r)
        if len(kept) >= N_TARGET:
            break

    out_path = Path(__file__).parent / "marginal_ionic_candidates.json"
    with open(out_path, "w") as f:
        json.dump(kept, f, indent=2)
    print(f"\nWrote {len(kept)} entries to {out_path}\n")

    print(f"{'class':<12}{'formula':<12}{'mp_id':<12}{'chemsys':<10}{'FE(eV/at)':>11}  "
          f"{'EAH':>8}  sg{'':<10}nsites  is_metal  theo")
    for r in kept:
        sg = r["spacegroup"] or "?"
        print(f"{r['anion_class']:<12}{r['formula']:<12}{r['mp_id']:<12}{r['chemsys']:<10}"
              f"{r['formation_energy_per_atom_eV']:>11.4f}  {r['energy_above_hull_eV_per_atom']:>8.4f}  "
              f"{sg:<12}{r['nsites']:<8}{str(r['is_metal']):<10}{r['theoretical']}")


if __name__ == "__main__":
    main()
