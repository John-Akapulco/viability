"""NOTE: bond_type/is_metal are sourced via build_bond_type_map()/
build_is_metal_map() (imported from compute_icobi_antibonding_all.py),
NOT read directly from mp_metadata.json -- the 186 main-campaign
compounds' own metadata carries neither field (only computed later, in
build_dataset.py, and stored in percolation_vs_antibonding.csv); reading
meta.get() directly here silently left "is_metal" and "bond_type" mostly
NaN for those compounds in an earlier version of this script (caught and
fixed the same day, see project memory).

Delta(antibonding-near-frontier), for the case-1 (decomposition-to-
elements) reaction, on BOTH the ICOHP and ICOBI antibonding-population
descriptors (analysis/compute_antibonding_all.py /
analysis/compute_icobi_antibonding_all.py + the maxhull top-ups).

Why this is NOT just reaction_analysis.delta.compute_delta() reused:
that module works on EXTENSIVE, already-formula-unit-summed ICOHP/ICOBI
values (CompoundEntry.icohp.sum_per_formula_unit_eV) and combines them
with integer stoichiometric coefficients -- a well-defined mass balance.
The antibonding-near-frontier descriptor is built differently: it
integrates LOBSTER's own "average" trace (already an INTENSIVE,
bond-averaged quantity for the whole structure, not a per-formula-unit
sum) over a fixed energy window. It does not have a natural "coefficient
x value" extensive decomposition the way a summed ICOHP does -- which is
exactly why mission #4 always correlated it directly against
formation_energy_per_atom (itself intensive), never against a raw
formation energy.

Chosen convention (documented here since it is a genuine judgment call,
not given by any prior module): for a compound decomposing into its
constituent elements, define

    Delta = [atomic-fraction-weighted average of the elements' OWN
             antibonding descriptor, each on its own reference structure]
            - [the compound's own antibonding descriptor]

i.e. "products - reactants" (same sign convention as reaction_analysis
throughout this project), but weighting products by atomic fraction
(intensive combination) rather than by integer stoichiometric
coefficient summed extensively (which would not be dimensionally
consistent with an intensive per-structure average).

Uses the same ELEMENT_REFERENCE mapping as
analysis/compute_reaction_icohp_case1.py (reused, not re-derived) so the
population is exactly the case-1 reaction set already established
there.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).parent))

import compute_reaction_icohp_case1 as ri_case1  # noqa: E402
from compute_icobi_antibonding_all import build_bond_type_map, build_is_metal_map  # noqa: E402

ELEMENT_REFERENCE = ri_case1.ELEMENT_REFERENCE

ICOHP_CSV = Path(__file__).parent / "percolation_vs_antibonding.csv"
ICOHP_MAXHULL_CSV = Path(__file__).parent / "icohp_antibonding_maxhull.csv"
ICOBI_CSV = Path(__file__).parent / "icobi_antibonding_all.csv"
OUT_CSV = Path(__file__).parent / "delta_antibonding_case1.csv"
FORMATION_ENERGIES = REPO_ROOT / "mp_dataset" / "formation_energies.json"


def _load_icohp_antibond() -> dict[str, float]:
    frames = [pd.read_csv(ICOHP_CSV)[["compound_id", "antibond_w_raw"]]]
    if ICOHP_MAXHULL_CSV.exists():
        frames.append(pd.read_csv(ICOHP_MAXHULL_CSV)[["compound_id", "antibond_w_raw"]])
    df = pd.concat(frames, ignore_index=True).dropna(subset=["antibond_w_raw"])
    return dict(zip(df["compound_id"], df["antibond_w_raw"]))


def _load_icobi_antibond() -> dict[str, float]:
    df = pd.read_csv(ICOBI_CSV)[["compound_id", "icobi_antibond_w_raw"]].dropna(subset=["icobi_antibond_w_raw"])
    return dict(zip(df["compound_id"], df["icobi_antibond_w_raw"]))


def main() -> None:
    import json

    bond_type_map = build_bond_type_map()
    is_metal_map = build_is_metal_map()
    icohp_map = _load_icohp_antibond()
    icobi_map = _load_icobi_antibond()
    formation_energies = json.loads(FORMATION_ENERGIES.read_text())
    print(f"ICOHP antibonding available for {len(icohp_map)} compounds")
    print(f"ICOBI antibonding available for {len(icobi_map)} compounds")

    rows = []
    n_ok = n_skip_single = n_skip_ref = n_skip_data = 0
    for compound_dir in sorted(ri_case1.STRUCTURES_ROOT.iterdir()):
        meta_path = compound_dir / "mp_metadata.json"
        if not meta_path.exists() or not ri_case1.is_computed(compound_dir):
            continue
        compound_id = compound_dir.name
        meta = json.loads(meta_path.read_text())
        formula = meta.get("formula")
        if not formula:
            continue
        try:
            comp = Composition(formula)
        except Exception:
            continue
        elements = [str(e) for e in comp.elements]
        if len(elements) == 1:
            n_skip_single += 1
            continue

        missing_refs = [e for e in elements if e not in ELEMENT_REFERENCE]
        if missing_refs:
            n_skip_ref += 1
            continue
        ref_ids = {e: ELEMENT_REFERENCE[e] for e in elements}

        total_atoms = sum(comp.get(e) for e in elements)
        frac = {e: comp.get(e) / total_atoms for e in elements}

        def _delta(desc_map: dict[str, float]) -> float | None:
            if compound_id not in desc_map:
                return None
            if any(ref_ids[e] not in desc_map for e in elements):
                return None
            products = sum(frac[e] * desc_map[ref_ids[e]] for e in elements)
            return products - desc_map[compound_id]

        d_icohp = _delta(icohp_map)
        d_icobi = _delta(icobi_map)
        if d_icohp is None and d_icobi is None:
            n_skip_data += 1
            continue

        n_ok += 1
        rows.append({
            "compound_id": compound_id,
            "mp_id": meta.get("mp_id"),
            "formula": formula,
            "family": meta.get("family"),
            "bond_type": bond_type_map.get(compound_id),
            "is_metal": is_metal_map.get(compound_id),
            "theoretical": meta.get("theoretical"),
            "energy_above_hull_eV_per_atom": meta.get("energy_above_hull_eV_per_atom"),
            "formation_energy_per_atom": formation_energies.get(meta.get("mp_id")) if meta.get("mp_id") else None,
            "delta_icohp_antibond": d_icohp,
            "delta_icobi_antibond": d_icobi,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\n{n_ok} compounds with a case-1 delta(antibonding) "
          f"(skipped: single-element={n_skip_single}, no-elemental-ref={n_skip_ref}, no-descriptor-data={n_skip_data})")
    print(f"  with delta_icohp_antibond: {df['delta_icohp_antibond'].notna().sum()}")
    print(f"  with delta_icobi_antibond: {df['delta_icobi_antibond'].notna().sum()}")
    print(f"Wrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
