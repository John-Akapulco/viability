"""Select candidates for a new campaign targeting the gap the project's
own Conclusion/Lecture (analysis/REPORT_delta_icohp_viability.md and
report/rapport_campagne2_fr.tex \\S\\ref{sec:reframing}) flags as the
concrete next step: almost every compound computed so far is
"comfortably stable or comfortably not" in `formation_energy_per_atom`
(decomposition into elements), which trivially routes
`classify.classify_viability()` to STABLE_ON_HULL and never exercises
the endobondic/exobondic sign discrimination the whole $\\Delta$(ICOHP)
framework is built around. That machinery only becomes non-trivial when
`delta_energy` (= -formation_energy_per_atom here) is close to the sign
flip -- i.e. compounds whose formation energy from the elements is
itself marginal, near zero, not their distance to the convex hull
(`energy_above_hull`), which is a different axis entirely and already
well represented in the dataset.

This is a *selection-only* script, same convention as
select_campaign.py/select_maxhull_binaries.py: writes
mp_dataset/marginal_formation_energy_candidates.json and prints a
summary table, does not download structures or launch any calculation.

Two bands away from (not touching) zero, reported separately (both are
"borderline" in different senses):
  - positive [FE_LO, FE_HI]: decomposition into elements is
    thermodynamically favorable but only marginally so -- exactly where
    a bonding-derived kinetic barrier (endobondic) could be the
    difference between "exists" and "doesn't", the case this project's
    own 7 reference reactions (Pb(N3)2, S4N2/S4N4, Mn2O7, ...) were all
    chosen to represent.
  - negative [-FE_HI, -FE_LO]: barely thermodynamically stable relative
    to elements -- still routes to STABLE_ON_HULL under
    classify_viability(), but useful as boundary controls on the other
    side of the same range.

A first attempt (2026-08-16) simply took the entries CLOSEST to FE=0
(0 to 0.10 eV/atom, sorted by |FE| ascending) and turned out to select
almost exclusively near-ideal metallic alloys (Mg-Zn, Mg-Ca, Cu-In,
...), whose FE~0 comes from weak/near-random mixing thermodynamics --
a different phenomenon from the kinetically-stabilized, non-metallic
reference compounds (azides, sulfur nitrides, oxides) this campaign is
actually meant to sample. Per user decision, fixed by moving the whole
window AWAY from zero (FE_LO/FE_HI below) rather than by excluding
is_metal -- metals stay eligible, but the near-zero alloy cluster is
now out of range entirely.

Constraints: non-magnetic (total magnetization window, same convention
as select_campaign.py/select_maxhull_binaries.py), deprecated=False,
nsites<=40, strictly binary (num_elements=2, same convention as every
prior campaign in this project -- select_campaign.py, extension4,
select_maxhull_binaries.py; an earlier version of this script allowed
ternaries too, corrected on user feedback), and -- the constraint
specific to this script --
every element must already be in ELEMENT_REFERENCE
(analysis/compute_reaction_icohp_case1.py, 62 elements) so every
selected compound's case-1 (decomposition-to-elements) reaction is
immediately computable with no new elemental-reference calculation
needed (avoids repeating the Bi2O3 gap from the maxhull batch, where 4
entries got no case-1 reaction for exactly this missing-reference
reason).

Diversity: at most 2 entries per chemical system (chemsys), preferring
the entry closest to FE=0 within each band, so the list isn't dominated
by many near-duplicate compounds of the same 1-2 elements.
"""

import json
import os
from collections import defaultdict
from pathlib import Path

from mp_api.client import MPRester

import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "analysis"))
from compute_reaction_icohp_case1 import ELEMENT_REFERENCE  # noqa: E402

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

MAX_SITES = 40
FE_LO = 0.02  # eV/atom, band starts this far from zero (excludes the near-ideal-alloy cluster)
FE_HI = 0.15  # eV/atom, band ends here (comparable order of magnitude to the 7 reference reactions)
MAX_PER_CHEMSYS = 2
N_TARGET_PER_BAND = 25

ALLOWED_ELEMENTS = frozenset(ELEMENT_REFERENCE.keys())

NON_MAGNETIC_WINDOW = (-0.01, 0.01)  # mu_B / formula unit

FIELDS = [
    "material_id",
    "formula_pretty",
    "chemsys",
    "elements",
    "energy_above_hull",
    "formation_energy_per_atom",
    "theoretical",
    "nsites",
    "nelements",
    "symmetry",
    "is_metal",
    "band_gap",
    "total_magnetization_normalized_formula_units",
]


def fetch(mpr: MPRester, *, fe_window):
    docs = mpr.materials.summary.search(
        num_elements=2,
        num_sites=(1, MAX_SITES),
        deprecated=False,
        formation_energy_per_atom=fe_window,
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    return [d for d in docs if ALLOWED_ELEMENTS.issuperset({str(e) for e in d.elements})]


def to_record(band, d):
    return {
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


def diversify(records, n_target):
    """At most MAX_PER_CHEMSYS entries per chemsys (experimental
    preferred within a chemsys, since theoretical-vs-experimental is
    itself a variable this project stratifies by), then evenly sampled
    across the full [FE_LO, FE_HI] range rather than biased toward
    either edge -- unlike the superseded closest-to-zero sort, every
    part of this band is equally "borderline" by construction, so no
    ordering bias is applied beyond chemsys diversity."""
    by_chemsys = defaultdict(list)
    for r in records:
        by_chemsys[r["chemsys"]].append(r)
    kept = []
    for chemsys, group in by_chemsys.items():
        group.sort(key=lambda r: r["theoretical"])  # False (experimental) first
        kept.extend(group[:MAX_PER_CHEMSYS])
    kept.sort(key=lambda r: (r["chemsys"], r["formula"]))
    if len(kept) > n_target:
        # Evenly spaced subsample across the sorted-by-chemsys list, not
        # a truncation, so systems from across the alphabet survive.
        step = len(kept) / n_target
        kept = [kept[int(i * step)] for i in range(n_target)]
    return kept


def main():
    with MPRester(API_KEY) as mpr:
        print("Fetching positive-FE (0, +threshold] candidates...")
        pos_docs = fetch(mpr, fe_window=(FE_LO, FE_HI))
        print(f"  pool={len(pos_docs)}")

        print("Fetching negative-FE [-threshold, 0) candidates...")
        neg_docs = fetch(mpr, fe_window=(-FE_HI, -FE_LO))
        print(f"  pool={len(neg_docs)}")

    pos_records = diversify([to_record("positive", d) for d in pos_docs], N_TARGET_PER_BAND)
    neg_records = diversify([to_record("negative", d) for d in neg_docs], N_TARGET_PER_BAND)
    selection = pos_records + neg_records

    out_path = Path(__file__).parent / "marginal_formation_energy_candidates.json"
    with open(out_path, "w") as f:
        json.dump(selection, f, indent=2)

    print(f"\nWrote {len(selection)} entries ({len(pos_records)} positive-FE, "
          f"{len(neg_records)} negative-FE) to {out_path}\n")

    print(f"{'band':<9}{'formula':<12}{'mp_id':<12}{'chemsys':<10}{'FE(eV/at)':>11}  "
          f"{'EAH':>8}  sg{'':<10}nsites  is_metal  theo")
    for r in selection:
        sg = r["spacegroup"] or "?"
        print(f"{r['band']:<9}{r['formula']:<12}{r['mp_id']:<12}{r['chemsys']:<10}"
              f"{r['formation_energy_per_atom_eV']:>11.4f}  {r['energy_above_hull_eV_per_atom']:>8.4f}  "
              f"{sg:<12}{r['nsites']:<8}{str(r['is_metal']):<10}{r['theoretical']}")


if __name__ == "__main__":
    main()
