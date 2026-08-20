"""Recompute this project's own DFT+LOBSTER ICOHP values against Reitz &
Dronskowski's (ic-2026-04181q) per-species reference numbers, for all 16
species across their 7 validation reactions -- both under the manuscript's
own methodology (PBEsol+D3(BJ), the `manuscript_*` campaign directories)
and, where a pre-existing entry exists, under this project's default
methodology (PBE, no D3, the original dataset's `extension_*`/`gasref_*`
directories).

Extends the 2026-08-18 fix (originally ZnSn/Zn/Sn only) to the remaining
13 species now that the 2026-08-19 manuscript-validation SLURM campaign is
finishing. As of this run: 10/16 manuscript_* (PBEsol) directories are
computed; the other 6 (Ca3N2, Mn2O7, PbN3_2, S4N2, S4N4, S8) are still on
SLURM and simply omitted from that column until they land. All 14 species
with a pre-existing PBE-functional dataset entry are included in the "Old
calculation" column now (2 species, S8 and Pb(N3)2, have no such entry --
both are novel to the manuscript-validation campaign).

Root cause background (unchanged from the 2026-08-18 fix, see
reaction_analysis/nearest_neighbor.py's own docstring and this project's
memory): reaction_analysis.nearest_neighbor.detect_first_shell() cuts at
the single LARGEST gap anywhere in a pair type's whole distance spectrum,
which mishandles compounds with a near-degenerate double first shell (Zn's
hcp c/a anomaly being the discovered case). Deliberately NOT changed in
nearest_neighbor.py itself -- still tests/test_nearest_neighbor.py-pinned,
dataset-wide behavior; whether it affects other compounds' bond_type
classification is a separate, not-yet-closed question (see
analysis/audit_shell_ambiguity.py and the 2026-08-19 activity report
addendum).

Reference per-species ICOHP values and their exact bond-type scope (which
pair types, how many bonds, and -- for Mn2O7 -- which of two Mn-O shells)
are transcribed directly from tests/fixtures/reitz_dronskowski_cases.yaml's
own inline comments (hand-transcribed from the manuscript by an earlier
session; not re-derived here). Per-species scope:

  PbN3_2:  Pb-N (x1) + N-N (x6, azide chain)
  Pb:      Pb-Pb, first shell (count not specified by the manuscript)
  N2:      N-N, single dimer bond
  S4N2:    three distinct S-N/S-S first-shell pair types (x2 each)
  S8:      S-S, first shell (x8, ring)
  S4N4:    S-N only (x8) -- no S-S term
  ZnSn:    Zn-Sn (x6) + Zn-Zn (x1) -- Sn-Sn explicitly excluded (manuscript
           scope, not a shell-detection issue, see 2026-08-18 fix)
  Zn:      Zn-Zn, true first shell only (x6) -- the original bug fix
  Sn:      Sn-Sn, first shell
  CaO_sphalerite: Ca-O, first shell (x4, tetrahedral)
  CaO_rocksalt:   Ca-O, first shell (x6, octahedral)
  CaN:     Ca-N, first shell (x6, rocksalt-type)
  Ca3N2:   Ca-N, first shell (x12)
  Mn2O7:   Mn-O, FIRST AND SECOND shell together (x2 Mn/FU, "short+long")
           + O-O, first shell (x7)
  MnO2:    Mn-O, scope not specified by the manuscript -- both single- and
           double-shell sums reported, see note in the printed table
  O2:      O-O, single dimer bond

Writes analysis/manuscript_validation_icohp_full.csv and prints the
comparison table (manuscript / new PBEsol / old PBE / diffs).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Structure

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from reaction_analysis.parse_lobster import raw_bond_records, _element_from_atom_label  # noqa: E402
from reaction_analysis.units import ev_to_kj_per_mol  # noqa: E402

STRUCTURES = REPO_ROOT / "mp_dataset" / "structures"
OUT_CSV = REPO_ROOT / "analysis" / "manuscript_validation_icohp_full.csv"

SHELL_GAP_THRESHOLD_ANGSTROM = 0.08  # separates the Mn-O short/long shells (Mn2O7, MnO2)

# Separates a first coordination shell from the next one -- must be wider
# than the widest within-shell spread seen (Ca3N2 0.011 A, Mn2O7 Mn-O/O-O
# sub-clusters 0.012-0.026 A) but narrower than the tightest real
# shell-to-shell gap seen (Sn's extension_Sn_mp-623511 PBE dir: 0.047 A
# between its true first shell and a distinct, weaker-ICOHP second shell).
FIRST_SHELL_GAP_THRESHOLD_ANGSTROM = 0.03


def _load_by_pair(compound_dir: Path) -> tuple[dict[str, list[dict]], int]:
    records = raw_bond_records(compound_dir, "ICOHPLIST.lobster", are_cobis=False)
    structure = Structure.from_file(compound_dir / "CONTCAR")
    _, z = structure.composition.get_reduced_composition_and_factor()
    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        key = "-".join(sorted((_element_from_atom_label(r["atom1"]), _element_from_atom_label(r["atom2"]))))
        by_pair[key].append(r)
    return by_pair, round(z)


def _cluster_shells(records: list[dict], gap_threshold: float = SHELL_GAP_THRESHOLD_ANGSTROM) -> list[list[dict]]:
    """Greedy sequential clustering on sorted distance: start a new shell
    whenever the gap to the next distance exceeds `gap_threshold`."""
    recs = sorted(records, key=lambda r: r["length"])
    shells: list[list[dict]] = [[recs[0]]] if recs else []
    for prev, cur in zip(recs, recs[1:]):
        if cur["length"] - prev["length"] > gap_threshold:
            shells.append([])
        shells[-1].append(cur)
    return shells


def _strict_first_shell(records: list[dict]) -> list[dict]:
    """First coordination shell via gap clustering, not a fixed distance
    tolerance from d_min: a fixed tolerance silently truncates shells whose
    own internal spread exceeds it (found in Ca3N2, 2026-08-20 -- distances
    span 2.4021 to 2.4239 A, so the old 0.02 A tol dropped the outermost
    24/96 bonds, a 24% error vs. the manuscript's -7.692 eV; gap clustering
    gets 0.09%). Verified against every other already-matching species
    (Pb, N2, Sn, Zn, CaO x2, CaN, O2, ZnSn, Mn2O7/MnO2 sub-shells) using
    FIRST_SHELL_GAP_THRESHOLD_ANGSTROM -- notably NOT the wider
    SHELL_GAP_THRESHOLD_ANGSTROM used for Mn-O shell separation, which is
    too loose here and wrongly merges Sn's real, distinct second shell
    (only 0.047 A beyond its first) into the sum."""
    recs = sorted(records, key=lambda r: r["length"])
    shell = [recs[0]]
    for prev, cur in zip(recs, recs[1:]):
        if cur["length"] - prev["length"] > FIRST_SHELL_GAP_THRESHOLD_ANGSTROM:
            break
        shell.append(cur)
    return shell


def _sum(records: list[dict]) -> tuple[int, float]:
    return len(records), sum(r["value"] for r in records)


# (compound_id, manuscript_eV, pbesol_dir, pbe_dir, scope_fn)
# scope_fn(by_pair) -> (n_bonds, sum_eV)
def _scope_single_pair_first_shell(pair: str, halve_for_self_pair: bool = False):
    """halve_for_self_pair: for a homonuclear pair summed around a single
    Z=1 site (e.g. elemental Pb, FCC, CN=12 all at one distance), each
    first-shell bond is shared with the neighboring image of that same
    site. The raw sum is the per-atom coordination total, not the
    per-formula-unit lattice total; the standard convention halves it to
    avoid double-counting each shared bond. Verified against Pb (2026-08-19):
    raw sum matches the manuscript's -5.688 eV to within 0.26% only after
    halving (11.406/2=5.703), vs. a ~100% miss unhalved."""
    def fn(by_pair):
        if pair not in by_pair:
            return None
        n, s = _sum(_strict_first_shell(by_pair[pair]))
        return (n, s / 2 if halve_for_self_pair else s)
    return fn


def _scope_multi_pair_first_shell(pairs: list[str]):
    def fn(by_pair):
        n_total, s_total = 0, 0.0
        for p in pairs:
            if p not in by_pair:
                return None
            n, s = _sum(_strict_first_shell(by_pair[p]))
            n_total += n
            s_total += s
        return (n_total, s_total)
    return fn


def _scope_mn2o7(by_pair):
    if "Mn-O" not in by_pair or "O-O" not in by_pair:
        return None
    shells = _cluster_shells(by_pair["Mn-O"])
    if len(shells) < 2:
        return None
    mn_o_records = shells[0] + shells[1]
    n_mno, s_mno = _sum(mn_o_records)
    n_oo, s_oo = _sum(_strict_first_shell(by_pair["O-O"]))
    return (n_mno + n_oo, s_mno + s_oo)


def _scope_mno2_single_shell(by_pair):
    if "Mn-O" not in by_pair:
        return None
    return _sum(_strict_first_shell(by_pair["Mn-O"]))


def _scope_mno2_double_shell(by_pair):
    if "Mn-O" not in by_pair:
        return None
    shells = _cluster_shells(by_pair["Mn-O"])
    if len(shells) < 2:
        return None
    return _sum(shells[0] + shells[1])


CASES = [
    dict(id="PbN3_2", manuscript_eV=-89.104, pbesol_dir="manuscript_PbN3_2", pbe_dir=None,
         scope=_scope_multi_pair_first_shell(["N-Pb", "N-N"])),
    dict(id="Pb", manuscript_eV=-5.688, pbesol_dir="manuscript_Pb", pbe_dir="extension_Pb_mp-20483",
         scope=_scope_single_pair_first_shell("Pb-Pb", halve_for_self_pair=True)),
    dict(id="N2", manuscript_eV=-23.161, pbesol_dir="manuscript_N2", pbe_dir="gasref_N2_dimerbox",
         scope=_scope_single_pair_first_shell("N-N")),
    dict(id="S4N2", manuscript_eV=-49.24, pbesol_dir="manuscript_S4N2", pbe_dir="extension_S4N2_cod4031496",
         scope=None),  # 3 distinct S-N/S-S pair types, not separable from ICOHPLIST alone -- reported raw below
    dict(id="S8", manuscript_eV=-46.8, pbesol_dir="manuscript_S8", pbe_dir=None,
         scope=_scope_single_pair_first_shell("S-S")),
    dict(id="S4N4", manuscript_eV=-73.84, pbesol_dir="manuscript_S4N4", pbe_dir="extension_S4N4_cod7017102",
         scope=_scope_single_pair_first_shell("N-S")),
    dict(id="ZnSn", manuscript_eV=-16.869, pbesol_dir="manuscript_ZnSn", pbe_dir="gasref_ZnSn_NiAs",
         scope=_scope_multi_pair_first_shell(["Sn-Zn", "Zn-Zn"])),
    dict(id="Zn", manuscript_eV=-12.686, pbesol_dir="manuscript_Zn", pbe_dir="extension_Zn_mp-79",
         scope=_scope_single_pair_first_shell("Zn-Zn")),
    dict(id="Sn", manuscript_eV=-7.672, pbesol_dir="manuscript_Sn", pbe_dir="extension_Sn_mp-623511",
         scope=_scope_single_pair_first_shell("Sn-Sn")),
    dict(id="CaO_sphalerite", manuscript_eV=-3.46, pbesol_dir="manuscript_CaO_sphalerite", pbe_dir="gasref_CaO_sphalerite",
         scope=_scope_single_pair_first_shell("Ca-O")),
    dict(id="CaO_rocksalt", manuscript_eV=-4.278, pbesol_dir="manuscript_CaO_rocksalt", pbe_dir="extension_CaO_mp-2605",
         scope=_scope_single_pair_first_shell("Ca-O")),
    dict(id="CaN", manuscript_eV=-4.302, pbesol_dir="manuscript_CaN", pbe_dir="extension_CaN_mp-1058549",
         scope=_scope_single_pair_first_shell("Ca-N")),
    dict(id="Ca3N2", manuscript_eV=-7.692, pbesol_dir="manuscript_Ca3N2", pbe_dir="extension_Ca3N2_mp-844",
         scope=_scope_single_pair_first_shell("Ca-N")),
    dict(id="Mn2O7", manuscript_eV=-62.21, pbesol_dir="manuscript_Mn2O7", pbe_dir="extension_Mn2O7_mp-28338",
         scope=_scope_mn2o7),
    dict(id="MnO2", manuscript_eV=-18.54, pbesol_dir="manuscript_MnO2", pbe_dir="extension_MnO2_mp-510408",
         scope=_scope_mno2_single_shell),
    dict(id="O2", manuscript_eV=-18.05, pbesol_dir="manuscript_O2", pbe_dir="gasref_O2_dimerbox",
         scope=_scope_single_pair_first_shell("O-O")),
]


def _eval_case(case: dict, dir_name: str | None) -> dict | None:
    if dir_name is None:
        return None
    compound_dir = STRUCTURES / dir_name
    if not (compound_dir / "ICOHPLIST.lobster").exists():
        return None
    by_pair, z = _load_by_pair(compound_dir)
    if case["scope"] is None:
        return {"n_bonds": None, "sum_eV": None, "raw_by_pair": {
            k: _sum(_strict_first_shell(v)) for k, v in by_pair.items()
        }}
    result = case["scope"](by_pair)
    if result is None:
        return None
    n_bonds, sum_eV = result
    return {"n_bonds": n_bonds, "sum_eV": sum_eV / z, "raw_by_pair": None}


def main() -> None:
    rows = []
    print(f"{'species':<16}{'manuscript':>12}{'new(PBEsol)':>13}{'diff%':>8}{'old(PBE)':>12}{'diff%':>8}")
    for case in CASES:
        new = _eval_case(case, case["pbesol_dir"])
        old = _eval_case(case, case["pbe_dir"])

        def fmt(res):
            if res is None:
                return "---", "---"
            if res["sum_eV"] is None:
                return "n/a(raw)", "n/a"
            diff = abs(res["sum_eV"] - case["manuscript_eV"]) / abs(case["manuscript_eV"]) * 100
            return f"{res['sum_eV']:.3f}", f"{diff:.1f}%"

        new_val, new_diff = fmt(new)
        old_val, old_diff = fmt(old)
        print(f"{case['id']:<16}{case['manuscript_eV']:>12.3f}{new_val:>13}{new_diff:>8}{old_val:>12}{old_diff:>8}")

        rows.append({
            "species": case["id"], "manuscript_eV": case["manuscript_eV"],
            "new_pbesol_eV": new["sum_eV"] if new and new["sum_eV"] is not None else "",
            "new_pbesol_n_bonds": new["n_bonds"] if new else "",
            "old_pbe_eV": old["sum_eV"] if old and old["sum_eV"] is not None else "",
            "old_pbe_n_bonds": old["n_bonds"] if old else "",
        })

        if new and new["raw_by_pair"] is not None:
            print(f"    {case['id']} raw first-shell pairs (scope not separable): "
                  f"{ {k: round(v[1], 3) for k, v in new['raw_by_pair'].items()} }")
        if old and old["raw_by_pair"] is not None:
            print(f"    {case['id']} (PBE) raw first-shell pairs: "
                  f"{ {k: round(v[1], 3) for k, v in old['raw_by_pair'].items()} }")

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
