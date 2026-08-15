"""Fourth extension batch: elemental decomposition references for the full
186-compound campaign + the 22-compound extension set, needed to compute
reaction_icohp.py case 1 (compound -> constituent elements) across the
whole dataset, not just the 8 already-computable extension reactions.

62 distinct elements appear across mp_dataset/structures/*/mp_metadata.json
(186 main campaign + extension + pilot). 12 already have a reference
compound in the dataset (Ca, C, K, Mn, N, Na, O, Pb, S, Si, Ti, Zn -- Ti
added just before this script, download_extension3.py). This script adds
the remaining 50.

Selection rule (base case): among MP's non-theoretical entries for that
element (num_elements=1), take the lowest energy_above_hull. This is a
BULK, RULE-BASED selection, unlike the earlier hand-picked extension
batches -- NOT individually literature-verified for every element. Where
the chosen entry sits within ~20 meV/at of another candidate, that's
recorded in the per-element note as an explicit "not individually
verified" flag (same self-skeptical spirit as flagging CaN/CaO's LOBSTER
quality issues elsewhere in this project) rather than silently presenting
a single confident answer.

7 elements got a MANUAL override, where the automated rule's pick
disagreed with a real, well-documented standard-state structure by a
DFT-noise-level energy margin (the same kind of PBE-GGA near-degeneracy
already documented for Ti alpha/omega, download_extension3.py): Ag (fcc),
In (tetragonal), Rb/Cs (bcc alkali), Se (trigonal), Sn (beta/white,
'tin pest'), Ta (bcc refractory). These 7 are flagged individually in
their own note text, not the generic rule-based one.

Fe, Co, Ni, Cr are real room-temperature magnetic elements (Fe/Co/Ni
ferromagnetic, Cr antiferromagnetic/SDW) -- same ISPIN=2 exception
mechanism as O2/MnO2/Mn in download_extension.py, added to
EXTENSION_SPIN_OVERRIDES in prepare_extension_vasp_lobster.py alongside
this script (simple collinear FM initial guesses, not an attempt at the
real magnetic structure, same documented caveat as alpha-Mn).

Writes one directory per element under mp_dataset/structures/, same
layout as the other extension_* batches, family="extension".
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

# (label, mp_id, note)
ELEMENTS = [
    ("Ag", "mp-124", "fcc Ag (Fm-3m), the universally-cited real standard state (Ag/Au/Cu/Ni/Pd/Pt/Al pattern); only 2.1 meV/at above DFT's nominal minimum mp-8566 (hcp-like P6_3/mmc), within typical GGA noise."),
    ("Al", "mp-134", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m); near-degenerate with ['mp-1183144'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("As", "mp-158", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Cmce); near-degenerate with ['mp-11'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("Au", "mp-81", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m); near-degenerate with ['mp-1008634'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("B", "mp-160", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=R-3m); near-degenerate with ['mp-161'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("Ba", "mp-122", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m); near-degenerate with ['mp-1096840', 'mp-10679', 'mp-56', 'mp-1058581'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("Be", "mp-87", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc)"),
    ("Br", "mp-23154", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Cmce)"),
    ("Cd", "mp-94", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc)"),
    ("Cl", "mp-22848", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Cmce)"),
    ("Co", "mp-102", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m); near-degenerate with ['mp-1271142', 'mp-1183710'] (both hcp-like) (within 20 meV/at) -- real Co is a famous close fcc/hcp call (hcp stable below ~422C, fcc above); NOT overridden, flagged instead of guessed"),
    ("Cr", "mp-90", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m) -- matches the real bcc standard state; Cr is a genuinely complex SDW antiferromagnet in reality, same caveat as alpha-Mn (simple collinear FM initial guess used, not the real magnetic structure)"),
    ("Cs", "mp-1", "bcc Cs (Im-3m), the real standard state (alkali metal pattern); DFT's nominal minima are a cluster of near-degenerate theoretical/exotic entries, none matching the well-established real bcc structure."),
    ("Cu", "mp-30", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m); near-degenerate with ['mp-989695', 'mp-989782'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("F", "mp-561203", "on-hull among non-theoretical MP entries (EAH=0.0031 eV/at, sg=C2/c); near-degenerate with ['mp-1525632', 'mp-760482', 'mp-561367'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default; periodic solid F2 proxy, same treatment as N2/O2/Cl2/Br2 elsewhere in this project"),
    ("Fe", "mp-13", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m) -- matches the real bcc (alpha-Fe) standard state; ferromagnetic, needs the ISPIN=2 exception (same mechanism as O2/MnO2/Mn)"),
    ("Ga", "mp-142", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Cmce); near-degenerate with ['mp-1067880'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("Ge", "mp-32", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fd-3m); near-degenerate with ['mp-1091415', 'mp-1007760'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("H", "mp-730101", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P2_12_12_1); near-degenerate with 7 other candidates (within 20 meV/at) -- not individually literature-verified; periodic solid H2 proxy, same molecular-solid treatment as N2/O2 elsewhere in this project"),
    ("Hf", "mp-103", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc)"),
    ("Hg", "mp-10861", "on-hull among non-theoretical MP entries (EAH=0.003 eV/at, sg=P6/mmm); near-degenerate with 7 other candidates (within 20 meV/at) -- not individually literature-verified; Hg is liquid at RT in reality, this is the low-T periodic solid proxy (same spirit as the N2/O2/Cl2/F2/Br2 gas-to-periodic-solid treatment)"),
    ("I", "mp-23153", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Cmce)"),
    ("In", "mp-1055994", "body-centered tetragonal In (I4/mmm), the real standard state (In's well-known tetragonally-distorted fcc); only 4.5 meV/at above DFT's nominal minimum mp-85 (ideal cubic Fm-3m), within typical GGA noise."),
    ("Ir", "mp-101", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m)"),
    ("Li", "mp-1018134", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=R-3m); near-degenerate with 6 other candidates (within 20 meV/at) -- not individually literature-verified; note this is plausibly the real low-T 9R Li ground state (bcc Li is a higher-temperature phase), not overridden"),
    ("Mg", "mp-153", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state"),
    ("Mo", "mp-129", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m) -- matches the real bcc standard state"),
    ("Nb", "mp-75", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m) -- matches the real bcc standard state; near-degenerate with ['mp-2739273'] (within 20 meV/at)"),
    ("Ni", "mp-23", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m) -- matches the real fcc standard state; ferromagnetic, needs the ISPIN=2 exception (same mechanism as O2/MnO2/Mn)"),
    ("Os", "mp-49", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state"),
    ("P", "mp-568348", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P2/c); near-degenerate with ['mp-1198724'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default; large cell (84 sites) -- accepted as the real DFT ground state despite compute cost, same precedent as PbN6/Ca3N2 elsewhere in this project"),
    ("Pd", "mp-2", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m) -- matches the real fcc standard state"),
    ("Pt", "mp-126", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m) -- matches the real fcc standard state"),
    ("Rb", "mp-70", "bcc Rb (Im-3m), the real standard state (alkali metals Li-Cs are all bcc at RT); only 1.8 meV/at above DFT's nominal minimum mp-639755 (I4/mmm), essentially degenerate."),
    ("Re", "mp-8", "on-hull among non-theoretical MP entries (EAH=0.0036 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state"),
    ("Rh", "mp-74", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Fm-3m) -- matches the real fcc standard state"),
    ("Ru", "mp-33", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state"),
    ("Sb", "mp-104", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=R-3m) -- matches the real rhombohedral standard state"),
    ("Sc", "mp-67", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state"),
    ("Se", "mp-14", "trigonal (gray, chain-structure) selenium, P3_121 -- the universally-cited real standard state; DFT's nominal minimum mp-570481 (P2_1/c, 64 sites) is only 1.1 meV/at lower, within typical DFT/k-mesh noise, not a credible distinct real ground state."),
    ("Sn", "mp-623511", "beta-Sn (white tin, I4/mmm), the real room-temperature standard state ('tin pest' alpha/beta transition at 13.2C); DFT's lowest-energy entry mp-117 (alpha/gray tin, Fd-3m) is real but only stable below 13.2C."),
    ("Sr", "mp-139", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc); near-degenerate with ['mp-1187073'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("Ta", "mp-50", "bcc Ta (Im-3m), the real standard state for this refractory metal (same structural family as Cr/V/Nb/Mo/W, all bcc); only 9.1 meV/at above DFT's nominal minimum mp-569794 (30-site P4_2/mnm), within typical GGA noise for this kind of near-degeneracy."),
    ("Tc", "mp-113", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc); near-degenerate with ['mp-867351'] (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("Te", "mp-19", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P3_121) -- matches the real trigonal standard state"),
    ("Tl", "mp-82", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc); near-degenerate with 4 other candidates (within 20 meV/at) -- not individually literature-verified, on-hull/lowest-EAH pick used by default"),
    ("V", "mp-146", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m) -- matches the real bcc standard state"),
    ("W", "mp-91", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=Im-3m) -- matches the real bcc standard state"),
    ("Y", "mp-112", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state"),
    ("Zr", "mp-131", "on-hull among non-theoretical MP entries (EAH=0.0 eV/at, sg=P6_3/mmc) -- matches the real hcp standard state; near-degenerate with ['mp-1077723'] (within 20 meV/at)"),
]


def main():
    OUT_ROOT.mkdir(exist_ok=True)
    mp_ids = [mp_id for _, mp_id, _ in ELEMENTS]
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
    for label, mp_id, note in ELEMENTS:
        d = by_id.get(mp_id)
        if d is None:
            print(f"MISSING {label} ({mp_id})")
            continue
        actual_sg = d.symmetry.symbol if d.symmetry else None

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
            "note": note,
        }
        (compound_dir / "mp_metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"OK  {compound_dir.name:35s} sg={actual_sg} nsites={d.nsites} "
              f"EAH={d.energy_above_hull:.4f} is_metal={d.is_metal}")
        n_ok += 1
    print(f"\n{n_ok}/{len(ELEMENTS)} element references downloaded into {OUT_ROOT}")


if __name__ == "__main__":
    main()
