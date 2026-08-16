"""Download the 45 structures selected by select_maxhull_binaries.py
(mp_dataset/maxhull_binaries_candidates.json, frozen and committed) --
15 experimentally known binary compounds sitting farthest above the
convex hull, each paired with up to 3 theoretical (never-synthesized)
polymorphs of the same formula sitting even farther above it. See that
script's docstring and analysis/REPORT_delta_icohp_viability.md's
"concrete next step" for the campaign's purpose.

Directory naming keeps the flat "extension_" prefix (not e.g.
"maxhull_") so prepare_extension_vasp_lobster.py's existing
auto-discovery (`compound_dir.name.startswith("extension_")`) picks
these up with no code change, same as every prior extension batch --
disambiguated by formula + kind + mp_id: extension_<formula>_<kind>_<mp_id>.

Selection was already frozen to maxhull_binaries_candidates.json by
select_maxhull_binaries.py; this script does not re-query MP for
selection, only fetches structures for the exact mp_ids already picked.

Magnetism: all 45 entries were selected within a
total_magnetization_normalized_formula_units window of +/-0.01 mu_B/f.u.
(same non-magnetic filter as select_campaign.py), so none are expected
to need EXTENSION_SPIN_OVERRIDES -- not verified against the override
dict here; check before running VASP if any compound looks off.
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
CANDIDATES_PATH = Path(__file__).parent / "maxhull_binaries_candidates.json"

KIND_SUFFIX = {"exp_anchor": "anchor", "theo_polymorph_above": "poly"}


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    picks = json.loads(CANDIDATES_PATH.read_text())
    mp_ids = [p["mp_id"] for p in picks]

    with MPRester(API_KEY) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=mp_ids,
            fields=[
                "material_id", "structure", "formula_pretty", "energy_above_hull",
                "theoretical", "nsites", "symmetry", "is_metal", "band_gap",
            ],
        )
    by_id = {str(d.material_id): d for d in docs}

    n_ok = 0
    for p in picks:
        mp_id = p["mp_id"]
        d = by_id.get(mp_id)
        if d is None:
            print(f"MISSING {p['formula']} {p['kind']} ({mp_id})")
            continue

        dirname = f"extension_{p['formula']}_{KIND_SUFFIX[p['kind']]}_{mp_id}"
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
            "batch": "maxhull_binaries_stress_test",
            "chemsys": p["chemsys"],
            "kind": p["kind"],
            "source": "materials_project",
            "bond_type": bond_type,
            "is_metal": d.is_metal,
            "band_gap_eV": d.band_gap,
            "energy_above_hull_eV_per_atom": d.energy_above_hull,
            "theoretical": d.theoretical,
            "nsites": d.nsites,
            "spacegroup": d.symmetry.symbol if d.symmetry else None,
            "anchor_formula": p.get("anchor_formula"),
            "anchor_energy_above_hull_eV_per_atom": p.get("anchor_energy_above_hull_eV_per_atom"),
            "note": (
                f"maxhull_binaries_stress_test batch: {p['kind']} for {p['chemsys']}, "
                f"see select_maxhull_binaries.py / maxhull_binaries_candidates.json "
                f"for the exact selection rule."
            ),
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {dirname:45s} sg={meta['spacegroup']} nsites={d.nsites} "
              f"EAH={d.energy_above_hull:.4f} is_metal={d.is_metal}")
        n_ok += 1

    print(f"\n{n_ok}/{len(picks)} maxhull_binaries structures downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
