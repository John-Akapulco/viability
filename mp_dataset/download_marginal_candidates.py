"""Download the 75 structures selected across
select_marginal_formation_energy.py (50, periodic-table-wide, mostly
metallic near-marginal-FE alloys) and select_marginal_ionic.py (30,
alkali/alkaline-earth nitrides and polyphosphides targeted specifically
because the wide query returned essentially no ionic-character
compounds), merged and chemsys-diversified down to 75 by
merge_marginal_candidates.py (mp_dataset/marginal_candidates_merged.json,
frozen and committed). See that script and the two selection scripts for
the full selection rationale: this campaign targets compounds with
formation_energy_per_atom NEAR the marginal boundary (not
energy_above_hull, a different axis) -- the axis that actually drives
classify_viability()'s STABLE_ON_HULL vs endobondic/exobondic branch,
which the project's Conclusion (report Sec. "Lecture") flags as
underrepresented in every campaign so far.

Directory naming keeps the flat "extension_" prefix (not
"marginal_" or similar) so prepare_extension_vasp_lobster.py's existing
auto-discovery (`compound_dir.name.startswith("extension_")`) picks
these up with no code change, same as every prior extension batch:
extension_<formula>_<mp_id> (no kind suffix needed -- unlike maxhull's
anchor/poly pairing, there is no natural pairing structure here, mp_id
alone is already unique).

Magnetism: all 75 entries were selected within the same non-magnetic
total_magnetization window as every other campaign in this project, so
none are expected to need EXTENSION_SPIN_OVERRIDES -- not verified
against the override dict here; checked at prepare-time instead
(prepare_extension_vasp_lobster.py reports MAGMOM/ISPIN for anything it
does set an override for).
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
CANDIDATES_PATH = Path(__file__).parent / "marginal_candidates_merged.json"


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
        compound_dir.mkdir(exist_ok=True)
        d.structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")

        elements = {str(e) for e in d.structure.composition.elements}
        bond_type = classify(elements, d.is_metal)

        meta = {
            "label": p["formula"],
            "mp_id": mp_id,
            "formula": d.formula_pretty,
            "family": "extension",
            "batch": "marginal_formation_energy_campaign",
            "campaign_source": p.get("source"),
            "chemsys": p["chemsys"],
            "source_query_note": p.get("anion_class") or p.get("band"),
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
                f"marginal_formation_energy_campaign: selected for "
                f"formation_energy_per_atom={d.formation_energy_per_atom:.4f} eV/atom "
                f"near the marginal boundary; see marginal_candidates_merged.json / "
                f"select_marginal_formation_energy.py / select_marginal_ionic.py for "
                f"the exact selection rule."
            ),
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {dirname:35s} sg={meta['spacegroup']} nsites={d.nsites} "
              f"FE={d.formation_energy_per_atom:.4f} is_metal={d.is_metal}")
        n_ok += 1

    print(f"\n{n_ok}/{len(picks)} marginal-candidate structures downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
