"""The three DeltaICOHP/DeltaICOBI normalizations, computed together
(never one in isolation) for a balanced Reaction.

All three follow the "products minus reactants" sign convention (like a
formation-energy/reaction-energy convention: Delta_rxn = E_products -
E_reactants) -- note this is the OPPOSITE sign convention from the
existing ad hoc `reaction_icohp.py` module (mission #5), which computes
reactant-minus-predicted-products. The two modules are not
interchangeable and their outputs must not be mixed without accounting
for the sign flip.

a) delta_per_formula_unit_eV (extensive, strict mass balance): the
   reaction's total Delta as balanced and written (i.e. with whatever
   coefficients the Reaction actually specifies -- typically smallest
   whole/rational numbers, one formula unit of the primary reactant/
   target, matching how formation energies are conventionally reported).
   This is the ICOHP/ICOBI analog of Delta_E_formation or Delta_E_hull.

b) delta_per_atom_eV = delta_per_formula_unit_eV / N_atoms_transferred,
   where N_atoms_transferred is the total atom count on either side of
   the balanced reaction (equal on both sides by construction). Lets
   reactions of very different formula-unit size be compared directly.

c) delta_per_bond_eV (intensive, DIAGNOSTIC ONLY, NOT a conservative
   reaction balance -- bonds are not conserved between differently-
   bonded structures): coefficient*n_bonds-weighted mean of
   mean_per_bond_eV across products, minus the same across reactants.
   Always returned with delta_per_bond_conservative=False.

If `reaction.bond_pair` is set and any involved CompoundEntry lacks that
bond type, the affected normalization(s) come back as NaN with an
explicit warning naming the compound and bond pair -- never a fatal
exception (same "batch continues" convention as percolation_path.py).
A malformed/unbalanced reaction likewise never raises out of
compute_delta(): the ReactionResult.error field carries the reason.

ReactionResult.bonding_label is classify.classify_bonding() applied to
delta_per_formula_unit_eV (Reitz & Dronskowski, ic-2026-04181q -- see
classify.py) -- endobondic/exobondic uses this module's own sign
convention directly (products - reactants), with no flip, since it
matches the manuscript's own ΔICOHP definition. None when
delta_per_formula_unit_eV is NaN (unbalanced/missing-data reaction).
"""

from __future__ import annotations

import math

from reaction_analysis.balance import ReactionBalanceError, check_balance
from reaction_analysis.classify import classify_bonding
from reaction_analysis.schema import CompoundEntry, IcohpSummary, Reaction, ReactionMember, ReactionResult

NAN = float("nan")


def _summary(entry: CompoundEntry, metric: str) -> IcohpSummary:
    if metric == "icohp":
        return entry.icohp
    if metric == "icobi":
        if entry.icobi is None:
            raise ValueError(f"{entry.compound_id!r} has no ICOBI data")
        return entry.icobi
    raise ValueError(f"Unknown metric {metric!r} (expected 'icohp' or 'icobi')")


def _value_per_fu(entry: CompoundEntry, metric: str, bond_pair: str | None) -> tuple[float | None, bool]:
    """(value per formula unit, availability). availability is False iff
    bond_pair was requested and this entry has no data for it."""
    summary = _summary(entry, metric)
    if bond_pair is None:
        return summary.sum_per_formula_unit_eV, True
    bt = summary.by_bond_type.get(bond_pair)
    if bt is None:
        return None, False
    return bt.sum_eV / entry.Z, True


def _mean_and_n_bonds(entry: CompoundEntry, metric: str, bond_pair: str | None) -> tuple[float | None, int | None, bool]:
    summary = _summary(entry, metric)
    if bond_pair is None:
        return summary.mean_per_bond_eV, summary.n_bonds, True
    bt = summary.by_bond_type.get(bond_pair)
    if bt is None:
        return None, None, False
    return bt.mean_eV, bt.n_bonds, True


def _side_total_per_fu(
    members: list[ReactionMember], entries: dict[str, CompoundEntry], metric: str, bond_pair: str | None,
    warnings: list[str],
) -> tuple[float, bool]:
    total = 0.0
    available = True
    for member in members:
        entry = entries[member.compound_id]
        value, ok = _value_per_fu(entry, metric, bond_pair)
        if not ok:
            available = False
            warnings.append(f"{entry.compound_id} has no bond_pair {bond_pair!r} data for metric {metric!r}")
            continue
        total += member.coefficient * value
    return total, available


def _side_weighted_mean_per_bond(
    members: list[ReactionMember], entries: dict[str, CompoundEntry], metric: str, bond_pair: str | None,
) -> float | None:
    total = 0.0
    total_weight = 0.0
    for member in members:
        entry = entries[member.compound_id]
        mean_val, n_bonds, ok = _mean_and_n_bonds(entry, metric, bond_pair)
        if not ok or not n_bonds:
            return None
        weight = member.coefficient * n_bonds
        total += weight * mean_val
        total_weight += weight
    if total_weight == 0:
        return None
    return total / total_weight


def _n_atoms_transferred(members: list[ReactionMember], entries: dict[str, CompoundEntry]) -> float:
    return sum(
        member.coefficient * sum(entries[member.compound_id].composition_per_formula_unit().values())
        for member in members
    )


def _compute_single_metric(reaction: Reaction, entries: dict[str, CompoundEntry], metric: str) -> ReactionResult:
    warnings: list[str] = []

    try:
        check_balance(reaction, entries)
    except ReactionBalanceError as exc:
        return ReactionResult(
            reaction_id=reaction.reaction_id, type=reaction.type, metric=metric, bond_pair=reaction.bond_pair,
            delta_per_formula_unit_eV=NAN, delta_per_atom_eV=NAN, delta_per_bond_eV=NAN,
            delta_per_bond_conservative=False, warnings=warnings, error=str(exc),
        )

    try:
        fu_reactants, ok_r = _side_total_per_fu(reaction.reactants, entries, metric, reaction.bond_pair, warnings)
        fu_products, ok_p = _side_total_per_fu(reaction.products, entries, metric, reaction.bond_pair, warnings)
        delta_fu = (fu_products - fu_reactants) if (ok_r and ok_p) else NAN

        n_atoms = _n_atoms_transferred(reaction.reactants, entries)
        delta_atom = (delta_fu / n_atoms) if (n_atoms > 0 and not math.isnan(delta_fu)) else NAN

        mean_reactants = _side_weighted_mean_per_bond(reaction.reactants, entries, metric, reaction.bond_pair)
        mean_products = _side_weighted_mean_per_bond(reaction.products, entries, metric, reaction.bond_pair)
        if mean_reactants is None or mean_products is None:
            delta_bond = NAN
            if reaction.bond_pair is not None:
                warnings.append(
                    f"delta_per_bond_eV undefined: bond_pair {reaction.bond_pair!r} missing on one or more members"
                )
        else:
            delta_bond = mean_products - mean_reactants

        bonding_label = None if math.isnan(delta_fu) else classify_bonding(delta_fu).value

        return ReactionResult(
            reaction_id=reaction.reaction_id, type=reaction.type, metric=metric, bond_pair=reaction.bond_pair,
            delta_per_formula_unit_eV=delta_fu, delta_per_atom_eV=delta_atom, delta_per_bond_eV=delta_bond,
            delta_per_bond_conservative=False, bonding_label=bonding_label, warnings=warnings, error=None,
        )
    except Exception as exc:  # noqa: BLE001 - batch must never die on one bad reaction
        return ReactionResult(
            reaction_id=reaction.reaction_id, type=reaction.type, metric=metric, bond_pair=reaction.bond_pair,
            delta_per_formula_unit_eV=NAN, delta_per_atom_eV=NAN, delta_per_bond_eV=NAN,
            delta_per_bond_conservative=False, warnings=warnings, error=f"{type(exc).__name__}: {exc}",
        )


def compute_delta(reaction: Reaction, entries: dict[str, CompoundEntry]) -> list[ReactionResult]:
    """One ReactionResult per requested metric (reaction.metric is
    "icohp", "icobi", or "icohp,icobi" -- the last one yields two
    results, never a single result with mixed fields)."""
    return [_compute_single_metric(reaction, entries, m) for m in reaction.metric.split(",")]
