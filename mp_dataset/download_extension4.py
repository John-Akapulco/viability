"""Fourth extension batch: alkali/alkaline-earth binary compounds against
N/O/F/P/S/Cl, sourced entirely from Materials Project (no COD, no "high
pressure" claim -- an earlier working draft of this request said "HP",
corrected by the user mid-session to "MP" i.e. Materials Project). For
each of the 60 possible (alkali or alkaline-earth) x (N/O/F/P/S/Cl)
binary chemical systems, MP was queried for every 2-element entry
(mp_dataset/extension4_candidates.json, built by a one-off selection
script, not hand-typed -- see below for why), and up to 2 compounds were
kept per system:

  - "exp_polymorph": an experimentally known (theoretical=False) polymorph
    NOT already anywhere in this dataset's mp_metadata.json (main
    campaign + earlier extension batches) -- deliberately targets the
    project's existing polymorph-density gap (only 8 same-composition
    polymorph groups existed before this batch; see
    analysis/REPORT_reaction_icohp.md's case-2 follow-up and
    reaction_analysis's own case-2 support) by preferring, for each
    system, whichever formula ALREADY has an entry in this dataset, so
    the new compound extends a real polymorph group rather than starting
    an unrelated one.
  - "theo_far_hull": a theoretical (DFT-only, no experimental report)
    entry, the single highest energy_above_hull available (nsites <= 40,
    for compute-cost control) for that same formula -- a genuinely
    far-from-equilibrium structure to stress-test the endobondic/exobondic
    and other descriptors, same spirit as the TiO2/carbon high-pressure-
    polymorph batch (download_extension2.py) but chosen by hull distance
    here rather than literature-matched space groups.

Selection was DETERMINISTIC and computed once (not by hand, unlike every
earlier extension batch's MP_SOURCED-style hardcoded tuple list) because
of scale (89 compounds across 44 systems) -- the exact list actually
downloaded is frozen to mp_dataset/extension4_candidates.json (committed
to the repo) specifically so this script stays reproducible without
depending on MP's index being unchanged; re-running the selection logic
against a live MP query is NOT what this script does.

Directory naming keeps the flat "extension_" prefix (not "extension4_")
so prepare_extension_vasp_lobster.py's existing auto-discovery
(`compound_dir.name.startswith("extension_")`) picks these up with no
code change, same as every prior extension batch -- disambiguated by
formula + kind + mp_id: extension_<formula>_<kind>_<mp_id>.

Magnetism: of these 89, only the two CsO2 entries (alkali superoxide,
open-shell O2^- radical, same physical class as the existing O2/MnO2/Mn
ISPIN=2 exceptions) need spin polarization -- added to
prepare_extension_vasp_lobster.py's EXTENSION_SPIN_OVERRIDES in the same
commit as this script. SrO2 was checked and is NOT an exception case:
it's the alkaline-earth PEROXIDE (Sr2+, closed-shell O2^2-), not a
superoxide, despite the superficial formula-pattern similarity to CsO2.
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
CANDIDATES_PATH = Path(__file__).parent / "extension4_candidates.json"

KIND_SUFFIX = {"exp_polymorph": "exp", "theo_far_hull": "theo"}


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    picks = json.loads(CANDIDATES_PATH.read_text())
    mp_ids = [p["material_id"] for p in picks]

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
        mp_id = p["material_id"]
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
            "batch": "extension4_alkali_alkaline_earth_binaries",
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
            "note": (
                f"extension4 batch: {p['kind']} pick for {p['chemsys']}, "
                f"selected by energy_above_hull among MP's {p['formula']} entries "
                f"(nsites<=40); see download_extension4.py / "
                f"extension4_candidates.json for the exact selection rule."
            ),
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {dirname:45s} sg={meta['spacegroup']} nsites={d.nsites} "
              f"EAH={d.energy_above_hull:.4f} is_metal={d.is_metal}")
        n_ok += 1

    print(f"\n{n_ok}/{len(picks)} extension4 compounds downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
