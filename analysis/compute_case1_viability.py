"""Combine each extension-campaign compound's case-1 (decomposition-to-
elements) delta_icohp -- already computed in reaction_analysis_case1_full.csv
-- with a matching delta_energy, then run both through
reaction_analysis.classify.classify_viability() for a real ViabilityLabel,
instead of the bonding-only endobondic/exobondic split the appendix tables
show today (per user request 2026-08-16: those tables were being read as
"exobondic == not viable", which classify.py explicitly warns against --
exobondic only means UNSTABLE_NONEXISTENT once delta_energy < 0 too; a
compound stable against decomposing into pure elements is STABLE_ON_HULL
regardless of its ICOHP sign).

delta_energy sourced from mp_dataset/formation_energies.json
(formation_energy_per_atom, MP's own DFT-referenced-to-elements
quantity) via delta_energy_per_atom = -formation_energy_per_atom. This
equivalence holds because formation_energy_per_atom is already defined
per atom as E_compound_per_atom - sum_i(x_i * E_element_i_per_atom), the
"reactants minus products" version of exactly the case-1 reaction (1
compound -> constituent elements); negating it gives the
"products - reactants" convention reaction_analysis uses throughout
(same convention as delta_icohp here), so no new DFT energies need
computing.

COD-sourced compounds (S4N2, S4N4: no mp_id, hence no MP formation
energy) are reported with viability_label="insufficient_data" rather
than fed a NaN into classify_viability, which has no NaN guard of its
own and would silently misclassify (NaN comparisons are always False in
Python, i.e. NaN >= 0 is False -- would look like EXOBONDIC/decomposition-
favorable, not "unknown").

ambiguous_ratio_threshold/exobondic_reference_magnitude are intentionally
left unset (classify_viability never defaults these) -- AMBIGUOUS_CHECK_
KINETICS is not produced by this script; that refinement is a separate,
not-yet-requested follow-up.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from reaction_analysis.classify import classify_viability  # noqa: E402

IN_CSV = REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv"
FORMATION_ENERGIES = REPO_ROOT / "mp_dataset" / "formation_energies.json"
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT_CSV = REPO_ROOT / "analysis" / "case1_viability.csv"

# The report's Table 8 (Fisher exact test, Sec. "Etendu a tout compose...")
# is scoped to the pooled 281-reaction population BEFORE the maxhull batch
# was folded into that particular analysis (only Sec. 5's dataset tables
# and the appendix compound list got the maxhull update so far -- see
# project memory). Excluded here too, so this script's counts line up
# EXACTLY with what Table 8 currently shows (281, not the full 322
# reaction_analysis_case1_full.csv now contains with maxhull included).
EXCLUDED_BATCHES = {"maxhull_binaries_stress_test"}


def _excluded_compound_ids() -> set[str]:
    excluded = set()
    for d in STRUCTURES_ROOT.iterdir():
        meta_path = d / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("batch") in EXCLUDED_BATCHES:
            excluded.add(d.name)
    return excluded


def main() -> None:
    formation_energies = json.loads(FORMATION_ENERGIES.read_text())
    excluded_ids = _excluded_compound_ids()

    with IN_CSV.open() as f:
        all_rows = list(csv.DictReader(f))
    rows = [r for r in all_rows if r["compound_id"] not in excluded_ids]

    print(f"{len(rows)} pooled case-1 reactions ({len(all_rows) - len(rows)} maxhull-batch rows excluded "
          f"to match the report's Table 8 population exactly)")

    out_rows = []
    n_by_label: dict[str, int] = {}
    for r in rows:
        mp_id = r["mp_id"]
        delta_icohp = float(r["delta_per_atom_eV"])
        fe = formation_energies.get(mp_id) if mp_id else None

        if fe is None or math.isnan(delta_icohp):
            label = "insufficient_data"
            bonding = None
            delta_energy = None
            warnings_str = "no MP formation_energy_per_atom (COD-sourced compound)" if fe is None else "delta_icohp is NaN"
        else:
            delta_energy = -float(fe)
            result = classify_viability(delta_energy=delta_energy, delta_icohp=delta_icohp)
            label = result.label.value
            bonding = result.bonding_label.value
            warnings_str = "; ".join(result.warnings)

        n_by_label[label] = n_by_label.get(label, 0) + 1
        out_rows.append({
            "compound_id": r["compound_id"],
            "mp_id": mp_id,
            "formula": r["formula"],
            "theoretical": r["theoretical"],
            "energy_above_hull_eV_per_atom": r["energy_above_hull_eV_per_atom"],
            "formation_energy_per_atom_eV": fe,
            "delta_energy_per_atom_eV": delta_energy,
            "delta_icohp_per_atom_eV": delta_icohp,
            "bonding_label": bonding,
            "viability_label": label,
            "warnings": warnings_str,
        })

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"\nWrote {len(out_rows)} rows to {OUT_CSV}")
    print("\nViabilityLabel counts:")
    for label, n in sorted(n_by_label.items(), key=lambda kv: -kv[1]):
        print(f"  {label:28s} {n}")

    print("\nViabilityLabel x theoretical (matches Table 8's experimental/theoretical split):")
    cross: dict[tuple[str, str], int] = {}
    for r in out_rows:
        key = (r["viability_label"], r["theoretical"])
        cross[key] = cross.get(key, 0) + 1
    for label in sorted(n_by_label):
        exp_n = cross.get((label, "False"), 0)
        theo_n = cross.get((label, "True"), 0)
        print(f"  {label:28s} experimental={exp_n:4d}  theoretical={theo_n:4d}")

    print("\nBondingLabel x ViabilityLabel (shows how many 'exobondic' are actually STABLE_ON_HULL):")
    bonding_cross: dict[tuple[str, str], int] = {}
    for r in out_rows:
        key = (r["bonding_label"] or "n/a", r["viability_label"])
        bonding_cross[key] = bonding_cross.get(key, 0) + 1
    for bl in ("endobondic", "exobondic", "n/a"):
        row_labels = {lbl: bonding_cross.get((bl, lbl), 0) for lbl in n_by_label}
        print(f"  {bl:12s} {row_labels}")


if __name__ == "__main__":
    main()
