"""Second hand-picked extension batch: high-pressure polymorphs computed here
at ambient pressure (0 GPa), to test whether the antibonding-population and
other descriptors pick up on the metastability of a structure whose stable
regime is a different thermodynamic condition than the one it's evaluated
under -- a different kind of "far from equilibrium" than the existing
extension_* set (those were chosen for element/chemistry reasons: magnetic,
missing from MP, etc; this batch is chosen for structural reasons).

Two sub-families, all sourced from Materials Project (no COD needed this
time -- all 8 entries exist in MP):

  - 3 experimentally-known high-pressure TiO2 polymorphs (rutile's
    ambient-pressure ground state is NOT included here -- only the phases
    only accessible experimentally under pressure, quenchable/metastable at
    ambient): TiO2-II (columbite/alpha-PbO2-type, Pbcn), baddeleyite-type
    (P2_1/c), TiO2-OII (cotunnite-type, Pnma). Identified among MP's ~90
    TiO2 polymorph entries by matching space group + lattice parameters
    against the experimental literature (Dubrovinsky et al. 2001 and
    related high-pressure XRD studies), not by name lookup -- MP carries no
    common-name field. A 4th named high-pressure phase, akaogeite/TiO2-OI,
    was searched for but not confidently found among MP's Pbca entries
    (only brookite, confirmed by matching its literature lattice
    parameters, was present) -- not included, flagged as a gap rather than
    guessed.

  - The carbon allotrope family: 3 real ambient-pressure phases (graphite,
    diamond, lonsdaleite) as the baseline, plus 2 theoretical superhard sp3
    allotropes predicted from cold-compressed graphite (M-carbon, W-carbon)
    -- these are computational predictions (MP theoretical=True, no direct
    single-crystal structure solution exists, only indirect XRD-pattern
    support from "superhard graphite" recovery experiments), identified by
    space group (M-carbon: C2/m, W-carbon: Pnma, per Li et al. 2009 PRL and
    Wang/Chen/Kawazoe 2011 PRL) cross-checked against sp3-diamond-like
    density (~5.8-5.9 A^3/atom) since MP reports these in a primitized
    monoclinic/orthorhombic setting that doesn't match the papers' raw
    conventional-cell lattice constants directly.

None of the 8 are magnetic (TiO2: Ti4+ is d0; carbon: closed-shell) --
no ISPIN=2 exception needed, unlike the O2/MnO2/Mn entries in the first
extension batch.

Writes one directory per compound under mp_dataset/structures/, same
layout as download_extension.py, family="extension" (this is still the
"hand-picked, not from the automated campaign selection" family).
"""

import json
import os
import sys
from pathlib import Path

from mp_api.client import MPRester

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from fetch_candidates import classify  # noqa: E402  (reused unmodified)

API_KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()
OUT_ROOT = Path(__file__).parent / "structures"

# (label, mp_id, expected_spacegroup, note)
MP_SOURCED = [
    ("TiO2II", "mp-1439", "Pbcn",
     "TiO2-II, columbite/alpha-PbO2-type -- experimental high-pressure "
     "polymorph (~4-10 GPa transition from rutile), quenchable to ambient; "
     "computed here at 0 GPa to test the descriptor against a structure "
     "whose stable regime is elsewhere. 0.0045 eV/at above hull."),
    ("TiO2baddeleyite", "mp-430", "P2_1/c",
     "baddeleyite-type TiO2 (ZrO2 structure), experimental high-pressure "
     "polymorph (~12-30 GPa); computed here at 0 GPa. 0.0395 eV/at above hull."),
    ("TiO2OII", "mp-9173", "Pnma",
     "TiO2-OII, cotunnite/PbCl2-type -- experimental ultra-high-pressure "
     "polymorph (~40+ GPa, historically claimed to rival diamond hardness, "
     "later disputed); computed here at 0 GPa. 0.0575 eV/at above hull. "
     "Akaogeite/TiO2-OI (the named phase between baddeleyite-type and this "
     "one in the real pressure sequence) was searched for but not "
     "confidently identified among MP's TiO2 Pbca entries -- gap, not "
     "included."),
    ("Cgraphite", "mp-48", "P6_3/mmc",
     "hexagonal (2H) graphite, real ambient-pressure ground-state-adjacent "
     "reference for the carbon family below. Semimetal (MP is_metal=True). "
     "0.0031 eV/at above hull."),
    ("Cdiamond", "mp-66", "Fd-3m",
     "diamond, real ambient-metastable sp3 reference. 0.1123 eV/at above hull."),
    ("Clonsdaleite", "mp-47", "P6_3/mmc",
     "lonsdaleite (hexagonal diamond), real but rare natural polymorph "
     "(meteorite impact sites). 0.1395 eV/at above hull."),
    ("CMcarbon", "mp-1080826", "C2/m",
     "M-carbon -- theoretical superhard sp3 allotrope predicted from "
     "cold-compressed graphite (Li et al. 2009 PRL); no confirmed direct "
     "crystal-structure solution, MP theoretical=True. Identified by space "
     "group + diamond-like density, not literature lattice-parameter match "
     "(MP's primitivized cell uses different axes than the paper's "
     "conventional C2/m cell). 0.3009 eV/at above hull -- the least stable "
     "compound in the whole extension set so far."),
    ("CWcarbon", "mp-1190171", "Pnma",
     "W-carbon -- theoretical superhard sp3 allotrope, same cold-compressed-"
     "graphite family as M-carbon (Wang, Chen & Kawazoe 2011 PRL, PRL 106, "
     "075501); MP theoretical=True, same density-based identification "
     "caveat as M-carbon. 0.2927 eV/at above hull."),
]


def download_mp_sourced() -> int:
    mp_ids = [mp_id for _, mp_id, _, _ in MP_SOURCED]
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
    for label, mp_id, expected_sg, note in MP_SOURCED:
        d = by_id.get(mp_id)
        if d is None:
            print(f"MISSING {label} ({mp_id})")
            continue
        actual_sg = d.symmetry.symbol if d.symmetry else None
        if expected_sg and actual_sg != expected_sg:
            print(f"WARNING {label} ({mp_id}): expected sg={expected_sg} but MP reports {actual_sg}")

        compound_dir = OUT_ROOT / f"extension_{label}_{mp_id}"
        compound_dir.mkdir(exist_ok=True)
        d.structure.to(filename=str(compound_dir / "POSCAR"), fmt="poscar")

        elements = {str(e) for e in d.structure.composition.elements}
        bond_type = classify(elements, d.is_metal)

        meta = {
            "label": label,
            "mp_id": mp_id,
            "formula": d.formula_pretty,
            "family": "extension",
            "source": "materials_project",
            "bond_type": bond_type,
            "is_metal": d.is_metal,
            "band_gap_eV": d.band_gap,
            "energy_above_hull_eV_per_atom": d.energy_above_hull,
            "theoretical": d.theoretical,
            "nsites": d.nsites,
            "spacegroup": actual_sg,
            "expected_spacegroup": expected_sg,
            "note": note,
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {compound_dir.name:35s} sg={actual_sg} nsites={d.nsites} "
              f"EAH={d.energy_above_hull:.4f} is_metal={d.is_metal}")
        n_ok += 1
    return n_ok


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    n_mp = download_mp_sourced()
    print(f"\n{n_mp}/{len(MP_SOURCED)} extension-batch-2 compounds downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
