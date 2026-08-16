"""Download the 119 structures selected by select_marginal_widen.py
(mp_dataset/marginal_widen_candidates.json) -- a widening of the
marginal-formation-energy campaign
(marginal_candidates_merged.json, already downloaded/computed, see
project memory) using the same FE band and chemsys-diversity discipline
but wider per-band/per-chemsys targets, plus two new chemistry-targeted
pools (organic C/H/N/O, molecular P/N/O) select_marginal_formation_energy.py
and select_marginal_ionic.py structurally could not reach. See that
script's docstring for the full selection rationale.

Same directory-naming and batch-tagging convention as
download_marginal_candidates.py, with its own batch value
("marginal_widen_campaign") so it stays excluded from any analysis
until deliberately folded in (see EXCLUDED_BATCHES in
analysis/compute_icohp_icobi_bondtype.py and
analysis/compute_case1_viability.py).
"""

import json
import os
from pathlib import Path

from mp_api.client import MPRester

import sys
sys.path.insert(0, str(Path(__file__).parent))
from fetch_candidates import classify  # noqa: E402  (reused unmodified)

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()
OUT_ROOT = Path(__file__).parent / "structures"
CANDIDATES_PATH = Path(__file__).parent / "marginal_widen_candidates.json"


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    picks = json.loads(CANDIDATES_PATH.read_text())
    mp_ids = [p["mp_id"] for p in picks]

    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=mp_ids,
            fields=[
                "material_id", "structure", "formula_pretty", "energy_above_hull",
                "formation_energy_per_atom", "theoretical", "nsites", "symmetry",
                "is_metal", "band_gap", "total_magnetization_normalized_formula_units",
            ],
        )
    by_id = {str(d.material_id): d for d in docs}

    n_ok = 0
    for p in picks:
        mp_id = p["mp_id"]
        d = by_id.get(mp_id)
        if d is None:
            print(f"MISSING {p['formula']} ({mp_id})")
            continue

        dirname = f"extension_{p['formula']}_{mp_id}"
        compound_dir = OUT_ROOT / dirname
        if compound_dir.exists():
            print(f"SKIP (dir exists) {dirname}")
            continue
        compound_dir.mkdir()
        d.structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")

        elements = {str(e) for e in d.structure.composition.elements}
        bond_type = classify(elements, d.is_metal)

        meta = {
            "label": p["formula"],
            "mp_id": mp_id,
            "formula": d.formula_pretty,
            "family": "extension",
            "batch": "marginal_widen_campaign",
            "campaign_source": p.get("pool"),
            "chemsys": p["chemsys"],
            "source_query_note": p.get("band"),
            "source": "materials_project",
            "bond_type": bond_type,
            "is_metal": d.is_metal,
            "band_gap_eV": d.band_gap,
            "energy_above_hull_eV_per_atom": d.energy_above_hull,
            "formation_energy_per_atom": d.formation_energy_per_atom,
            "theoretical": d.theoretical,
            "nsites": d.nsites,
            "spacegroup": d.symmetry.symbol if d.symmetry else None,
            "note": (
                f"marginal_widen_campaign (pool={p.get('pool')}): selected for "
                f"formation_energy_per_atom={d.formation_energy_per_atom:.4f} eV/atom "
                f"near the marginal boundary; see marginal_widen_candidates.json / "
                f"select_marginal_widen.py for the exact selection rule."
            ),
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {dirname:35s} sg={meta['spacegroup']} nsites={d.nsites} "
              f"FE={d.formation_energy_per_atom:.4f} is_metal={d.is_metal}")
        n_ok += 1

    print(f"\n{n_ok}/{len(picks)} widened marginal-candidate structures downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
