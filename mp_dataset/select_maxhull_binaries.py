"""Select candidates for a new campaign: the experimental binary compounds
sitting *farthest* above the thermodynamic convex hull (the most extreme
"real but metastable" edge cases Materials Project has), each paired with
2-3 hypothetical (theoretical, never-synthesized) polymorphs of the same
formula that sit even farther above the hull -- a deliberate stress-test
population for the `ViabilityLabel`/endobondic-exobondic discrimination
question (see analysis/REPORT_delta_icohp_viability.md's "concrete next
step": extend deliberately to borderline-stability compounds, not just
whatever the next unrelated mission adds).

This is a *selection-only* script, same convention as select_campaign.py
and download_extension4.py's candidate-freeze step: it writes
mp_dataset/maxhull_binaries_candidates.json and prints a summary table,
but does not download structures or launch any calculation -- that is an
explicitly separate, not-yet-authorized next step.

Constraints (both the experimental anchor and its theoretical polymorphs):
strictly binary (num_elements=2), non-magnetic (total magnetization
window, same convention as select_campaign.py -- the categorical
`magnetic_ordering` field is unreliable), no f-block elements,
deprecated=False, nsites<=40 (extension4's ceiling, looser than the main
campaign's 20 since these are one-off additions, not a 60-wide family).

Diversity: at most one experimental anchor per chemical system (chemsys),
so the top-N list is not dominated by many polymorphs of the same 1-2
elements -- ranked by energy_above_hull descending, take the single
highest-E_hull experimental entry per chemsys, then the top N systems
overall.
"""

import json
import os
from collections import defaultdict
from pathlib import Path

from mp_api.client import MPRester

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

MAX_SITES = 40
N_ANCHORS = 15
N_POLYMORPHS_PER_ANCHOR = 3

LANTHANIDES = [
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu",
]
ACTINIDES = [
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr",
]
EXCLUDE_F_ELEMENTS = frozenset(LANTHANIDES + ACTINIDES)

NON_MAGNETIC_WINDOW = (-0.01, 0.01)  # mu_B / formula unit

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
    "is_metal",
    "band_gap",
    "total_magnetization_normalized_formula_units",
]


def fetch(mpr: MPRester, *, theoretical: bool):
    docs = mpr.materials.summary.search(
        num_elements=2,
        num_sites=(1, MAX_SITES),
        theoretical=theoretical,
        deprecated=False,
        energy_above_hull=(0.0005, None),  # strictly metastable, off the hull
        total_magnetization_normalized_formula_units=NON_MAGNETIC_WINDOW,
        fields=FIELDS,
    )
    return [d for d in docs if not (EXCLUDE_F_ELEMENTS & {str(e) for e in d.elements})]


def to_record(kind, d, anchor_formula=None, anchor_hull=None):
    rec = {
        "kind": kind,
        "mp_id": str(d.material_id),
        "formula": d.formula_pretty,
        "chemsys": d.chemsys,
        "energy_above_hull_eV_per_atom": d.energy_above_hull,
        "theoretical": d.theoretical,
        "nsites": d.nsites,
        "spacegroup": d.symmetry.symbol if d.symmetry else None,
        "is_metal": d.is_metal,
        "band_gap_eV": d.band_gap,
        "total_magnetization_per_fu": d.total_magnetization_normalized_formula_units,
    }
    if anchor_formula is not None:
        rec["anchor_formula"] = anchor_formula
        rec["anchor_energy_above_hull_eV_per_atom"] = anchor_hull
    return rec


def main():
    with MPRester(API_KEY) as mpr:
        print("Fetching experimental metastable binaries...")
        exp_docs = fetch(mpr, theoretical=False)
        print(f"  pool={len(exp_docs)}")

        print("Fetching theoretical metastable binaries...")
        theo_docs = fetch(mpr, theoretical=True)
        print(f"  pool={len(theo_docs)}")

    theo_by_formula = defaultdict(list)
    for d in theo_docs:
        theo_by_formula[d.formula_pretty].append(d)

    # An experimental entry only qualifies as an anchor if AT LEAST ONE
    # theoretical polymorph of the SAME formula sits farther above the
    # hull than it does (per user instruction 2026-08-16: an anchor with
    # only 1 or 2 -- not the full 3 -- qualifying polymorphs above it is
    # still kept; the earlier all-or-nothing >=2 cutoff was too strict).
    # An anchor with ZERO polymorphs above it is dropped, not because
    # it's uninteresting on its own, but because this campaign's whole
    # point is the paired exp-anchor + theoretical-polymorphs-above
    # construction -- an anchor with nothing above it can't be that.
    # Picking the literal global-maximum experimental entry per chemsys
    # (tried first, unfiltered) is self-defeating for this purpose: those
    # entries are so extreme (several eV/atom for several MP entries --
    # itself a data-quality flag, see the report) that essentially
    # nothing in MP's own theoretical search sits above them.
    MIN_POLYMORPHS_REQUIRED = 1

    qualifying_by_chemsys = defaultdict(list)
    for d in exp_docs:
        n_above = sum(
            1 for t in theo_by_formula.get(d.formula_pretty, [])
            if t.energy_above_hull > d.energy_above_hull
        )
        if n_above >= MIN_POLYMORPHS_REQUIRED:
            qualifying_by_chemsys[d.chemsys].append(d)

    # One anchor per chemsys: among qualifying entries, the one with the
    # highest energy_above_hull (still "maximal distance", now subject to
    # the room-above constraint).
    best_by_chemsys = {
        chemsys: max(docs, key=lambda d: d.energy_above_hull)
        for chemsys, docs in qualifying_by_chemsys.items()
    }

    anchors = sorted(best_by_chemsys.values(), key=lambda d: -d.energy_above_hull)[:N_ANCHORS]

    selection = []
    summary_rows = []
    for anchor in anchors:
        selection.append(to_record("exp_anchor", anchor))

        candidates = [
            d for d in theo_by_formula.get(anchor.formula_pretty, [])
            if d.energy_above_hull > anchor.energy_above_hull
        ]
        candidates.sort(key=lambda d: -d.energy_above_hull)
        picks = candidates[:N_POLYMORPHS_PER_ANCHOR]
        for p in picks:
            selection.append(to_record(
                "theo_polymorph_above", p,
                anchor_formula=anchor.formula_pretty,
                anchor_hull=anchor.energy_above_hull,
            ))

        summary_rows.append((anchor, picks))

    out_path = Path(__file__).parent / "maxhull_binaries_candidates.json"
    with open(out_path, "w") as f:
        json.dump(selection, f, indent=2)

    print(f"\nWrote {len(selection)} entries ({len(anchors)} experimental anchors + "
          f"{len(selection) - len(anchors)} theoretical polymorphs) to {out_path}\n")

    print(f"{'formula':<14}{'mp_id':<12}{'chemsys':<10}{'EAH(eV/at)':>12}  "
          f"sg{'':<10}nsites  is_metal  #polymorphs_above")
    for anchor, picks in summary_rows:
        sg = anchor.symmetry.symbol if anchor.symmetry else "?"
        print(f"{anchor.formula_pretty:<14}{str(anchor.material_id):<12}{anchor.chemsys:<10}"
              f"{anchor.energy_above_hull:>12.4f}  {sg:<12}{anchor.nsites:<8}"
              f"{str(anchor.is_metal):<10}{len(picks)}")


if __name__ == "__main__":
    main()
