"""Endobondic/exobondic bonding classification and the (deliberately more
cautious) viability read built on top of it, per Reitz & Dronskowski
(ic-2026-04181q) -- see reaction_analysis package docs / mission notes for
the manuscript's own definitions. This module never imports pymatgen or
touches real LOBSTER/DFT data; it operates purely on a delta_icohp (and,
for the viability label, a delta_energy) already computed elsewhere
(delta.py).

Sign convention: delta_icohp = ICOHP_sum(products) - ICOHP_sum(reactants)
(the manuscript's own convention, and exactly this package's delta.py
convention already -- no sign flip needed at this boundary, unlike the
reaction_icohp.py module elsewhere in this project).

  delta_icohp > 0 -> ENDOBONDIC: breaking the reactant's bonds costs more
    than the products' bonds recover -- a bonding-derived kinetic barrier
    exists, so an exothermic decomposition can still be kinetically
    blocked (the compound is metastable/viable despite delta_energy < 0).
  delta_icohp < 0 -> EXOBONDIC: no such barrier is visible from bonding
    alone.

Critical caveat repeated from the manuscript (their Mn2O7 discussion) and
enforced structurally below, NOT as a hardcoded chemistry rule: an
exobondic sign never proves a compound cannot exist. It only means this
static, bonding-only heuristic finds no barrier -- slow, non-bond-breaking
decomposition pathways (e.g. Mn2O7's gradual O2 loss) are invisible to it.
UNSTABLE_NONEXISTENT is therefore ALWAYS returned together with a warning
saying exactly this; classify_viability() never hardcodes Mn2O7 or any
other specific compound.
"""

from __future__ import annotations

from enum import Enum
from typing import NamedTuple, Optional


class BondingLabel(str, Enum):
    ENDOBONDIC = "endobondic"
    EXOBONDIC = "exobondic"


class ViabilityLabel(str, Enum):
    METASTABLE_VIABLE = "metastable_viable"
    UNSTABLE_NONEXISTENT = "unstable_nonexistent"
    STABLE_ON_HULL = "stable_on_hull"
    AMBIGUOUS_CHECK_KINETICS = "ambiguous_check_kinetics"


_KINETICS_CAVEAT = (
    "exobondic sign means no bonding-derived kinetic barrier was found, "
    "but per Reitz & Dronskowski (ic-2026-04181q) this never proves the "
    "compound cannot exist -- slow, non-bond-breaking decomposition "
    "pathways (e.g. Mn2O7's gradual O2 loss) are not captured by this "
    "static sign criterion; check the literature before concluding "
    "UNSTABLE_NONEXISTENT."
)


def classify_bonding(delta_icohp: float) -> BondingLabel:
    """delta_icohp == 0.0 exactly is classified ENDOBONDIC (no bonding
    driving force for decomposition found), the same side as > 0 -- there
    is no third label for the zero case in the manuscript."""
    return BondingLabel.ENDOBONDIC if delta_icohp >= 0 else BondingLabel.EXOBONDIC


class ViabilityResult(NamedTuple):
    label: ViabilityLabel
    bonding_label: BondingLabel
    warnings: list[str]


def classify_viability(
    delta_energy: float,
    delta_icohp: float,
    *,
    ambiguous_ratio_threshold: Optional[float] = None,
    exobondic_reference_magnitude: Optional[float] = None,
) -> ViabilityResult:
    """Combine a thermodynamic delta (delta_energy: products - reactants,
    same sign convention as delta_icohp; e.g. a DFT reaction energy or
    formation-energy difference computed by the caller -- this module
    never computes energies itself) with delta_icohp to produce a
    ViabilityLabel.

      delta_energy >= 0  -> STABLE_ON_HULL (does not spontaneously
        decompose in the first place; bonding is not even consulted).
      delta_energy < 0 and endobondic -> METASTABLE_VIABLE.
      delta_energy < 0 and exobondic  -> UNSTABLE_NONEXISTENT, always with
        the kinetics caveat warning (see module docstring).

    AMBIGUOUS_CHECK_KINETICS is never triggered from delta_energy/
    delta_icohp alone (the manuscript gives no static rule for it -- their
    own Mn2O7 call was a literature judgment, not derived from the sign or
    magnitude of ΔICOHP). It can only be requested explicitly by the
    caller supplying BOTH `ambiguous_ratio_threshold` and
    `exobondic_reference_magnitude` (e.g. the typical |delta_icohp| across
    a caller-chosen comparison set) -- if |delta_icohp| falls below
    ambiguous_ratio_threshold * exobondic_reference_magnitude, the result
    is downgraded from UNSTABLE_NONEXISTENT to AMBIGUOUS_CHECK_KINETICS
    with an additional warning. Neither parameter is ever guessed or
    defaulted to a nonzero value by this module.
    """
    warnings: list[str] = []

    if delta_energy >= 0:
        bonding = classify_bonding(delta_icohp)
        return ViabilityResult(ViabilityLabel.STABLE_ON_HULL, bonding, warnings)

    bonding = classify_bonding(delta_icohp)
    if bonding is BondingLabel.ENDOBONDIC:
        return ViabilityResult(ViabilityLabel.METASTABLE_VIABLE, bonding, warnings)

    warnings.append(_KINETICS_CAVEAT)

    if (
        ambiguous_ratio_threshold is not None
        and exobondic_reference_magnitude is not None
        and exobondic_reference_magnitude > 0
        and abs(delta_icohp) < ambiguous_ratio_threshold * exobondic_reference_magnitude
    ):
        warnings.append(
            f"|delta_icohp|={abs(delta_icohp):.4g} is below "
            f"{ambiguous_ratio_threshold:g} x the supplied reference exobondic "
            f"magnitude ({exobondic_reference_magnitude:.4g}) -- flagged "
            "ambiguous rather than confidently non-existent."
        )
        return ViabilityResult(ViabilityLabel.AMBIGUOUS_CHECK_KINETICS, bonding, warnings)

    return ViabilityResult(ViabilityLabel.UNSTABLE_NONEXISTENT, bonding, warnings)
