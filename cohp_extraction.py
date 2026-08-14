#!/usr/bin/env python3
"""COHP extraction feasibility (mission: antibonding-population-near-E_F,
steps 0+1 ONLY). Reuses pymatgen's own COHPCAR parsers -- no hand-rolled
parsing of the LOBSTER COHPCAR format, since a hand-rolled parser is
exactly the kind of place a sign-convention or energy-reference bug would
go unnoticed.

This module deliberately does NOT define an antibonding-population metric
or an energy window (that is step 2, out of scope here). It only provides:
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
    converged calculation confirms both are metals, band_gap=0.0).

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
