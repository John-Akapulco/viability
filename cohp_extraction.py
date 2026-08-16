#!/usr/bin/env python3
"""COHP extraction + antibonding-population-near-frontier metric.

Steps 0+1 (extraction, cross-validation, metal/gap classification) plus
step 2 (window + metric definition), validated on the 6 pilot compounds
only -- extension to the full 186-compound dataset is a separate,
not-yet-authorized mission, per the same pattern used for every other
descriptor in this project (validate small, then decide whether to scale).

Reuses pymatgen's own COHPCAR parsers throughout -- no hand-rolled parsing
of the LOBSTER COHPCAR format, since a hand-rolled parser is exactly the
kind of place a sign-convention or energy-reference bug would go
unnoticed. Provides:
  - loading COHPCAR.lobster via pymatgen (both the raw Cohpcar reader and
    the higher-level CompleteCohp wrapper);
  - a cross-validation check against the already-validated ICOHPLIST.lobster
    values that percolation_path.py's own tests rely on;
  - metal-vs-gap classification, cross-checked against Materials Project's
    own is_metal/band_gap rather than derived from a coarse local k-mesh
    (see the "known pitfall" note on AlNi/BeCu in
    analysis/REPORT_cohp_feasibility.md: a naive local
    eigenvalue_band_properties/DOS-at-Ef check on our LOBSTER-oriented
    k-mesh spuriously suggests small gaps for both, when MP's own
    converged calculation confirms both are metals, band_gap=0.0);
  - antibonding_population_near_frontier(): the step-2 metric. Window and
    reference point are chosen per compound (E_F for metals, VBM for
    gapped compounds -- see analysis/METRIC_DEFINITION_antibonding.md for
    the full rationale), one-sided (only occupied states, since only
    occupied antibonding character is energetically destabilizing in the
    ground state), integrating the positive (antibonding) part of the
    "average" COHP trace.

Uses pymatgen 2026.5.4 (see requirements.txt).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.outputs import Cohpcar, Icohplist


def load_cohpcar(compound_dir: Path) -> Cohpcar:
    """Raw per-line COHPCAR.lobster reader (pymatgen.io.lobster.outputs.Cohpcar).
    Energies as written by LOBSTER are already E - E_F (LOBSTER's own file
    convention, confirmed empirically below -- not a pymatgen-side shift)."""
    return Cohpcar(filename=str(compound_dir / "COHPCAR.lobster"))


def load_complete_cohp(compound_dir: Path) -> CompleteCohp:
    """Higher-level pymatgen.electronic_structure.cohp.CompleteCohp wrapper,
    which additionally needs the structure (CONTCAR) to attach site/bond
    information. Its internal `energies` array is on an absolute-eV scale
    (efermi is a nonzero absolute value, not 0) -- a different convention
    from the raw Cohpcar's E-E_F-shifted array, but both are internally
    self-consistent (see cross_validate_against_icohplist, which checks
    both against the same ICOHPLIST.lobster reference and gets an exact
    match through either)."""
    return CompleteCohp.from_file(
        fmt="LOBSTER",
        filename=str(compound_dir / "COHPCAR.lobster"),
        structure_file=str(compound_dir / "CONTCAR"),
    )


def cross_validate_against_icohplist(compound_dir: Path) -> Dict[str, object]:
    """For every bond label in ICOHPLIST.lobster, compare the ICOHP value
    there against the ICOHP trace in COHPCAR.lobster evaluated at E=E_F
    (raw Cohpcar, energies pre-shifted by LOBSTER so E_F is at 0). Returns
    per-label differences and summary statistics. Both files are read
    independently through separate pymatgen parsers (Icohplist vs Cohpcar)
    so this is a genuine cross-check, not comparing a value against itself.
    """
    icohplist = Icohplist(filename=str(compound_dir / "ICOHPLIST.lobster"))
    cohpcar = load_cohpcar(compound_dir)

    energies = np.array(cohpcar.energies)
    idx_ef = int(np.argmin(np.abs(energies - 0.0)))
    energy_at_idx = float(energies[idx_ef])

    diffs = {}
    for label in icohplist.icohplist:
        reference = icohplist.icohpcollection.get_icohp_by_label(label)
        if label not in cohpcar.cohp_data:
            diffs[label] = {"reference": reference, "from_cohpcar": None, "diff": None}
            continue
        spins = list(cohpcar.cohp_data[label]["ICOHP"].keys())
        from_cohpcar = sum(cohpcar.cohp_data[label]["ICOHP"][s][idx_ef] for s in spins)
        diffs[label] = {
            "reference": reference,
            "from_cohpcar": float(from_cohpcar),
            "diff": float(from_cohpcar - reference),
        }

    finite_diffs = [abs(d["diff"]) for d in diffs.values() if d["diff"] is not None]
    return {
        "n_labels": len(diffs),
        "n_matched": len(finite_diffs),
        "energy_at_ef_index": energy_at_idx,
        "max_abs_diff": max(finite_diffs) if finite_diffs else None,
        "mean_abs_diff": float(np.mean(finite_diffs)) if finite_diffs else None,
        "per_label": diffs,
    }


def sign_convention_check(compound_dir: Path, label: str = "average") -> Dict[str, object]:
    """Empirical (not assumed) sign-convention statement: confirms that
    negative ICOHP/COHP = bonding, positive = antibonding, by checking
    that pymatgen's own Cohp.has_antibnd_states_below_efermi (which
    internally tests `cohp > positive_limit` for "antibonding") operates
    on the same-signed values as ICOHPLIST.lobster's already-established
    convention (most negative = strongest bond, used throughout this
    project since the pilot report)."""
    complete_cohp = load_complete_cohp(compound_dir)
    cohp_obj = complete_cohp.all_cohps[label] if label != "average" else complete_cohp
    antibnd = cohp_obj.has_antibnd_states_below_efermi(spin=Spin.up, limit=0.01)
    return {
        "label": label,
        "has_antibonding_states_below_efermi": antibnd,
        "convention": (
            "Confirmed by inspection of pymatgen's own "
            "Cohp.has_antibnd_states_below_efermi implementation (cohp_vals > "
            "positive limit => antibonding) applied to the same array that "
            "cross_validate_against_icohplist() shows is numerically "
            "identical (same sign) to the ICOHPLIST.lobster convention "
            "already used throughout this project: NEGATIVE ICOHP/COHP = "
            "bonding, POSITIVE = antibonding."
        ),
    }


def metal_or_gap_from_mp(mp_id: str, api_key_path: str = "~/.mp_api_key") -> Optional[Dict[str, object]]:
    """Metal/gap classification from Materials Project's own converged
    calculation, NOT derived from our local LOBSTER-oriented k-mesh -- see
    module docstring for why the local mesh is unreliable for this
    specific question (confirmed empirically for AlNi/BeCu)."""
    import os

    from mp_api.client import MPRester

    api_key = open(os.path.expanduser(api_key_path)).read().strip()
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=[mp_id], fields=["material_id", "is_metal", "band_gap"]
        )
    if not docs:
        return None
    d = docs[0]
    return {"mp_id": mp_id, "is_metal": d.is_metal, "band_gap": d.band_gap}


# ---------------------------------------------------------------------------
# Step 2: antibonding-population-near-frontier window + metric.
# See analysis/METRIC_DEFINITION_antibonding.md for the full rationale.
# ---------------------------------------------------------------------------

DEFAULT_DELTA_E = 1.0  # eV, one-sided window width below the reference energy


def frontier_reference_energy(
    compound_dir: Path, is_metal: bool, vasprun_path: Optional[Path] = None
) -> float:
    """Reference energy E_ref on the SAME axis as Cohpcar.energies (which
    LOBSTER writes already shifted so E_F = 0):
      - metal: E_ref = 0.0 (E_F itself; COHPCAR's own zero).
      - gapped: E_ref = VBM_absolute(VASP) - E_F_absolute(LOBSTER). Both
        absolute values live on the same VASP-internal eigenvalue scale
        for a given compound (confirmed empirically: VASP's and LOBSTER's
        own reported E_F agree to <1e-4 eV across all 6 pilots -- see
        analysis/METRIC_DEFINITION_antibonding.md). Requires vasprun.xml
        (not COHPCAR-only) since VBM is a band-structure quantity LOBSTER
        itself does not report.
      - `is_metal` MUST come from an external, reliable source (Materials
        Project's own is_metal, per metal_or_gap_from_mp) -- NOT derived
        locally, since a naive local VBM/CBM or DOS-at-E_F check on this
        project's LOBSTER-oriented coarse k-mesh is confirmed unreliable
        for exactly this question (AlNi/BeCu spuriously read as gapped;
        see analysis/REPORT_cohp_feasibility.md).
    """
    if is_metal:
        return 0.0

    import warnings

    from pymatgen.io.vasp import Vasprun

    if vasprun_path is None:
        vasprun_path = compound_dir / "vasprun.xml"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vr = Vasprun(str(vasprun_path), parse_dos=True, parse_eigen=True, parse_potcar_file=False)
        _, _, vbm, _ = vr.eigenvalue_band_properties
        vasp_efermi = vr.efermi

    lobster_efermi = load_cohpcar(compound_dir).efermi
    if abs(vasp_efermi - lobster_efermi) > 1e-2:
        raise ValueError(
            f"VASP E_F ({vasp_efermi:.4f}) and LOBSTER E_F ({lobster_efermi:.4f}) disagree "
            f"by more than 0.01 eV for {compound_dir} -- do not trust the VBM alignment "
            "without investigating first (see analysis/METRIC_DEFINITION_antibonding.md)."
        )
    return vbm - lobster_efermi


def integrate_antibonding_in_window(
    energies: np.ndarray, cohp_values: np.ndarray, e_ref: float, delta_e: float
) -> float:
    """Pure numerical core of the step-2 metric, factored out from
    antibonding_population_near_frontier() so it can be validated on
    hand-crafted synthetic arrays (tests/test_cohp_extraction.py) with a
    known analytic answer, independent of any real-file parsing.

    Integrates max(cohp_values, 0) (the antibonding part only) over the
    one-sided window (e_ref - delta_e, e_ref] via the trapezoidal rule on
    whatever energy grid is supplied.
    """
    energies = np.asarray(energies)
    cohp_values = np.asarray(cohp_values)
    mask = (energies > e_ref - delta_e) & (energies <= e_ref)
    if not mask.any():
        raise ValueError(
            f"Window ({e_ref - delta_e:.4f}, {e_ref:.4f}] contains no grid points "
            f"(grid spacing ~{energies[1]-energies[0]:.4f}) -- delta_e too small."
        )
    antibonding_only = np.clip(cohp_values[mask], a_min=0.0, a_max=None)
    return float(np.trapezoid(antibonding_only, energies[mask]))


def antibonding_population_near_frontier(
    compound_dir: Path,
    is_metal: bool,
    delta_e: float = DEFAULT_DELTA_E,
    label: str = "average",
    vasprun_path: Optional[Path] = None,
) -> Dict[str, object]:
    """Integrate the antibonding (positive) part of COHP(E) over the
    one-sided window (E_ref - delta_e, E_ref], E_ref = E_F (metals) or VBM
    (gapped compounds), on the `label` bond trace (default: LOBSTER's own
    "average" over all bonds).

    Only occupied states can be antibonding-destabilizing in the ground
    state, hence the one-sided window (below the reference only), not a
    window straddling it.

    Returns both the raw integral (LOBSTER's ICOHP-style units, i.e. what
    LOBSTER itself calls "eV" by convention though it's a Hamilton-
    population integral, not literally an energy) and a version normalized
    by the total occupied ICOHP magnitude at E_ref on the same trace (the
    same normalized/raw-pair pattern used throughout this project for the
    percolation weight and the min-cut descriptor).
    """
    cohpcar = load_cohpcar(compound_dir)
    e_ref = frontier_reference_energy(compound_dir, is_metal, vasprun_path)

    energies = np.array(cohpcar.energies)
    spins = list(cohpcar.cohp_data[label]["COHP"].keys())
    cohp_total = sum(cohpcar.cohp_data[label]["COHP"][s] for s in spins)
    icohp_total = sum(cohpcar.cohp_data[label]["ICOHP"][s] for s in spins)

    w_antibond = integrate_antibonding_in_window(energies, cohp_total, e_ref, delta_e)
    n_grid_points_in_window = int(((energies > e_ref - delta_e) & (energies <= e_ref)).sum())

    idx_ref = int(np.argmin(np.abs(energies - e_ref)))
    total_occupied_magnitude = abs(float(icohp_total[idx_ref]))
    w_antibond_normalized = (
        w_antibond / total_occupied_magnitude if total_occupied_magnitude > 0 else None
    )

    return {
        "label": label,
        "is_metal": is_metal,
        "e_ref": e_ref,
        "delta_e": delta_e,
        "window": [e_ref - delta_e, e_ref],
        "n_grid_points_in_window": n_grid_points_in_window,
        "w_antibond_raw": w_antibond,
        "total_occupied_icohp_magnitude": total_occupied_magnitude,
        "w_antibond_normalized": w_antibond_normalized,
    }


# ---------------------------------------------------------------------------
# Step 2b: the same near-frontier antibonding-population metric, on ICOBI
# instead of ICOHP.
#
# CRITICAL, EMPIRICALLY VERIFIED SIGN-CONVENTION DIFFERENCE: ICOBI does NOT
# share ICOHP's "negative = bonding" convention -- it uses the OPPOSITE
# sign, "positive = bonding" (like COOP, not like COHP). Verified two ways
# on real project data (not assumed from memory or documentation):
#   1. cross_validate against ICOBILIST.lobster (see
#      tests/test_cohp_extraction.py's ICOBI cross-check): pymatgen's
#      Cohpcar(are_cobis=True) reproduces every ICOBILIST.lobster value
#      exactly (max abs diff 0.0 across 66 bond labels on the SrSi pilot),
#      confirming pymatgen does NOT silently renormalize the sign when
#      switching from COHP to COBI parsing -- the raw LOBSTER file
#      convention passes through unchanged.
#   2. On extension_SiO2_anchor_mp-9258's strongest Si-O bond
#      (ICOHP=-5.599, strongly bonding by the ICOHP convention,
#      ICOBI=+0.513, a large bond order): the raw per-energy curves in the
#      deep bonding region (-10 to -3 eV, far below any frontier effects)
#      average COHP(E)=-0.338 (negative, bonding, as expected) but
#      COBI(E)=+0.041 (POSITIVE in that same bonding-dominated region).
# Consequently, the antibonding (destabilizing) part of an ICOBI curve is
# its NEGATIVE part, not its positive part -- the clip direction below is
# deliberately the mirror image of integrate_antibonding_in_window()
# above, not a copy-paste of it.
# ---------------------------------------------------------------------------


def load_cobicar(compound_dir: Path) -> Cohpcar:
    """Raw per-line COBICAR.lobster reader. Same file format/energy
    convention as load_cohpcar() (E - E_F already shifted by LOBSTER), just
    are_cobis=True so pymatgen parses the COBI/ICOBI columns instead of
    COHP/ICOHP -- see module-level note above on why the sign meaning of
    the resulting "positive"/"negative" values is NOT the same as for
    load_cohpcar()."""
    return Cohpcar(filename=str(compound_dir / "COBICAR.lobster"), are_cobis=True)


def integrate_icobi_antibonding_in_window(
    energies: np.ndarray, cobi_values: np.ndarray, e_ref: float, delta_e: float
) -> float:
    """ICOBI analog of integrate_antibonding_in_window(): integrates
    abs(min(cobi_values, 0)) -- the NEGATIVE (antibonding, per the sign
    note above) part of the COBI(E) trace -- over the one-sided window
    (e_ref - delta_e, e_ref], via the trapezoidal rule. Returns a
    non-negative magnitude, same convention as the ICOHP version, so the
    two are directly comparable in sign/interpretation despite integrating
    opposite-signed raw regions of their respective source curves.
    """
    energies = np.asarray(energies)
    cobi_values = np.asarray(cobi_values)
    mask = (energies > e_ref - delta_e) & (energies <= e_ref)
    if not mask.any():
        raise ValueError(
            f"Window ({e_ref - delta_e:.4f}, {e_ref:.4f}] contains no grid points "
            f"(grid spacing ~{energies[1]-energies[0]:.4f}) -- delta_e too small."
        )
    antibonding_only = np.clip(cobi_values[mask], a_min=None, a_max=0.0)
    return float(np.trapezoid(np.abs(antibonding_only), energies[mask]))


def icobi_antibonding_population_near_frontier(
    compound_dir: Path,
    is_metal: bool,
    delta_e: float = DEFAULT_DELTA_E,
    label: str = "average",
    vasprun_path: Optional[Path] = None,
) -> Dict[str, object]:
    """ICOBI analog of antibonding_population_near_frontier() -- same
    window, same E_ref logic (frontier_reference_energy() is metric-
    agnostic: it only depends on the electronic structure, not on which
    curve is being integrated), same one-sided-window physical rationale
    (only occupied states can destabilize the ground state), but
    integrating COBICAR.lobster's NEGATIVE part instead of COHPCAR.lobster's
    positive part (see the sign-convention note above
    integrate_icobi_antibonding_in_window()). Returns a raw magnitude
    (dimensionless bond-order units, not eV -- ICOBI itself is
    dimensionless) plus a version normalized by the total occupied |ICOBI|
    magnitude at E_ref on the same trace, mirroring the ICOHP version's
    raw/normalized pair.
    """
    cobicar = load_cobicar(compound_dir)
    e_ref = frontier_reference_energy(compound_dir, is_metal, vasprun_path)

    energies = np.array(cobicar.energies)
    spins = list(cobicar.cohp_data[label]["COHP"].keys())
    cobi_total = sum(cobicar.cohp_data[label]["COHP"][s] for s in spins)
    icobi_total = sum(cobicar.cohp_data[label]["ICOHP"][s] for s in spins)

    w_antibond = integrate_icobi_antibonding_in_window(energies, cobi_total, e_ref, delta_e)
    n_grid_points_in_window = int(((energies > e_ref - delta_e) & (energies <= e_ref)).sum())

    idx_ref = int(np.argmin(np.abs(energies - e_ref)))
    total_occupied_magnitude = abs(float(icobi_total[idx_ref]))
    w_antibond_normalized = (
        w_antibond / total_occupied_magnitude if total_occupied_magnitude > 0 else None
    )

    return {
        "label": label,
        "is_metal": is_metal,
        "e_ref": e_ref,
        "delta_e": delta_e,
        "window": [e_ref - delta_e, e_ref],
        "n_grid_points_in_window": n_grid_points_in_window,
        "w_antibond_icobi_raw": w_antibond,
        "total_occupied_icobi_magnitude": total_occupied_magnitude,
        "w_antibond_icobi_normalized": w_antibond_normalized,
    }
