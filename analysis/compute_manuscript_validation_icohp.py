"""Recompute the ZnSn/Zn/Sn rows of report Table S18 (the DFT+LOBSTER ICOHP
comparison against Reitz & Dronskowski, ic-2026-04181q) under a strict
single-shell convention, diagnosing and resolving the Zn 36% residual left
open by the 2026-08-18 manuscript-methodology campaign.

Root cause (confirmed on mp_dataset/structures/manuscript_Zn/ICOHPLIST.lobster):
reaction_analysis.nearest_neighbor.detect_first_shell() cuts at the single
LARGEST gap anywhere in a pair-type's whole distance spectrum -- correct and
already regression-tested for the general 588-compound bond_type
classification pipeline (analysis/compute_icohp_icobi_bondtype.py), where
species pairs typically have one well-separated first shell. Zn's hcp
lattice is the counter-example this general rule mishandles: its
anomalously low c/a ratio puts a second, chemically distinct Zn-Zn shell
(interlayer, 2.729-2.734 A) only ~0.12 A past the true first shell (basal,
2.608 A) -- much closer than the ~1.04 A gap to the real third shell. The
"largest gap" heuristic locks onto that later gap and lumps both Zn-Zn
shells into one, roughly doubling the summed ICOHP. This is a real limit of
applying one general-purpose algorithm dataset-wide, not a genuine
DFT/LOBSTER methodology difference -- see report/appendix_reitz_dronskowski
_structures_{en,fr}.tex for the write-up. Deliberately NOT changed in
nearest_neighbor.py itself: that module's "largest gap" rule is
tests/test_nearest_neighbor.py-pinned behavior serving the whole dataset,
and Zn's near-degenerate double shell is not representative of the general
case -- changing it there would need its own dataset-wide validation pass,
out of scope here (see project memory: real goal is descriptor validation,
not classification refinement).

Fix applied HERE ONLY (a one-off, manuscript-comparison-specific
recomputation): cluster each pair-type's distances within `tol` (0.02 A) of
the minimum observed distance -- large enough to swallow the ~1e-4 A
numerical splitting LOBSTER reports for symmetry-equivalent bonds, small
enough to stay well inside every real single-shell cluster seen in these
three compounds (largest observed intra-shell spread: 0.034 A, ZnSn's
Sn-Sn) and well short of every observed inter-shell gap (smallest: 0.121 A,
Zn-Zn). Not a claim that a flat 0.02 A cutoff is right for other pair types
or other compounds -- see nearest_neighbor.py's own docstring for why no
project-wide module uses an absolute cutoff.

ZnSn also needed a second correction, already footnoted in the table before
this script existed but not previously applied to the headline number: the
manuscript's own accounting for ZnSn is 6 Zn-Sn + 1 Zn-Zn bonds per formula
unit, with NO Sn-Sn term, even though a real, well-separated first-shell
Sn-Sn contact exists in the DFT+LOBSTER data (3.54-3.57 A). This is a
scope-of-accounting difference, not a shell-detection error, so it is
handled here simply by excluding the Sn-Sn pair type from ZnSn's sum.
"""

from __future__ import annotations

import sys
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Structure

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from reaction_analysis.parse_lobster import raw_bond_records, _element_from_atom_label  # noqa: E402
from reaction_analysis.units import ev_to_kj_per_mol  # noqa: E402

STRUCTURES = REPO_ROOT / "mp_dataset" / "structures"

STRICT_SHELL_TOL_ANGSTROM = 0.02

# (compound_dir_name, manuscript_eV_per_FU, pair_types_to_include)
CASES = [
    ("manuscript_Zn", -12.686, ["Zn-Zn"]),
    ("manuscript_Sn", -7.672, ["Sn-Sn"]),
    ("manuscript_ZnSn", -16.869, ["Sn-Zn", "Zn-Zn"]),  # Sn-Sn excluded, see module docstring
]


def _load_by_pair(compound_dir: Path) -> tuple[dict[str, list[dict]], int]:
    records = raw_bond_records(compound_dir, "ICOHPLIST.lobster", are_cobis=False)
    structure = Structure.from_file(compound_dir / "CONTCAR")
    _, z = structure.composition.get_reduced_composition_and_factor()
    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        key = "-".join(sorted((_element_from_atom_label(r["atom1"]), _element_from_atom_label(r["atom2"]))))
        by_pair[key].append(r)
    return by_pair, round(z)


def _strict_first_shell_sum(records: list[dict], tol: float = STRICT_SHELL_TOL_ANGSTROM) -> tuple[int, float]:
    d_min = min(r["length"] for r in records)
    kept = [r for r in records if r["length"] <= d_min + tol]
    return len(kept), sum(r["value"] for r in kept)


def main() -> None:
    results: dict[str, float] = {}
    print(f"{'compound':<18}{'pair types':<16}{'n_bonds':>8}{'sum/FU (eV)':>14}{'manuscript':>13}{'diff %':>9}")
    for name, manuscript_eV, pair_types in CASES:
        by_pair, z = _load_by_pair(STRUCTURES / name)
        n_total = 0
        sum_total = 0.0
        for pt in pair_types:
            n, s = _strict_first_shell_sum(by_pair[pt])
            n_total += n
            sum_total += s
        per_fu = sum_total / z
        results[name] = per_fu
        diff_pct = abs(per_fu - manuscript_eV) / abs(manuscript_eV) * 100
        print(f"{name:<18}{','.join(pair_types):<16}{n_total:>8}{per_fu:>14.5f}{manuscript_eV:>13.3f}{diff_pct:>8.2f}%")

    delta_eV = (results["manuscript_Zn"] + results["manuscript_Sn"]) - results["manuscript_ZnSn"]
    delta_kJ = ev_to_kj_per_mol(delta_eV)
    manuscript_delta_eV, manuscript_delta_kJ = -3.489, -337
    diff_pct = abs(delta_eV - manuscript_delta_eV) / abs(manuscript_delta_eV) * 100
    print(f"\nDelta ICOHP, ZnSn -> Zn + Sn: {delta_eV:.3f} eV ({delta_kJ:.1f} kJ/mol) "
          f"vs manuscript {manuscript_delta_eV} eV ({manuscript_delta_kJ} kJ/mol), diff {diff_pct:.2f}%")


if __name__ == "__main__":
    main()
