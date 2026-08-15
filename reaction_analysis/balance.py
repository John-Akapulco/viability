"""Stoichiometric balance checking for a Reaction, and coefficient
auto-derivation for the decomposition-into-elements case.

Balance is checked on PER-FORMULA-UNIT composition
(CompoundEntry.composition_per_formula_unit(), i.e. composition/Z), not
on raw composition (atoms in the calculated cell). Raw cell atom counts
are an artifact of how large a supercell VASP/LOBSTER happened to run on
and are not comparable across compounds -- a ReactionMember's coefficient
is, as in a normal written chemical equation, a multiplier on formula
units. This is not spelled out explicitly in the original data-format
request; it is the only interpretation under which "coefficient x
composition" is dimensionally meaningful across compounds of different
cell sizes, so it is made explicit and tested here rather than left
ambiguous.
"""

from __future__ import annotations

from reaction_analysis.schema import CompoundEntry, Reaction, ReactionMember

_TOL = 1e-6


class ReactionBalanceError(ValueError):
    """Raised when a Reaction's reactant and product sides do not balance
    element-by-element. Always names the specific unbalanced element(s)
    and the imbalance -- never a generic failure."""


def element_balance(reaction: Reaction, entries: dict[str, CompoundEntry]) -> dict[str, float]:
    """Per-element (reactant total - product total) atom count, in
    formula-unit-normalized units. A balanced reaction has every value
    at (approximately) zero. Does not raise -- callers that want a hard
    check should use check_balance()."""
    diff: dict[str, float] = {}

    def _accumulate(members: list[ReactionMember], sign: float) -> None:
        for member in members:
            entry = entries[member.compound_id]
            for el, n_per_fu in entry.composition_per_formula_unit().items():
                diff[el] = diff.get(el, 0.0) + sign * member.coefficient * n_per_fu

    _accumulate(reaction.reactants, +1.0)
    _accumulate(reaction.products, -1.0)
    return diff


def check_balance(reaction: Reaction, entries: dict[str, CompoundEntry]) -> None:
    """Raise ReactionBalanceError, naming every unbalanced element and
    its imbalance, if the reaction does not balance to within floating-
    point tolerance. Silent (returns None) if it does."""
    diff = element_balance(reaction, entries)
    unbalanced = {el: d for el, d in diff.items() if abs(d) > _TOL}
    if unbalanced:
        detail = ", ".join(f"{el}: {d:+.6f}" for el, d in sorted(unbalanced.items()))
        raise ReactionBalanceError(
            f"Reaction '{reaction.reaction_id}' is not balanced (reactants - products, "
            f"per formula unit, should be ~0 for every element): {detail}"
        )


def derive_element_coefficients(
    target: CompoundEntry, element_entries: dict[str, CompoundEntry]
) -> list[ReactionMember]:
    """Case 1 helper: given a target CompoundEntry and a dict of
    element-symbol -> pure-element CompoundEntry (role="element",
    single-species composition), derive the ReactionMember coefficients
    that balance `target -> sum(elements)` -- avoids manual entry of
    coefficients that are fully determined by the target's own
    composition.

    Each element entry's own cell may hold more than one atom (e.g. a
    diatomic N2/O2 reference, or any multi-atom elemental structure) --
    the coefficient is target_atoms_of_el / element_entry's
    atoms-per-formula-unit of that element, so that
    coefficient * element_entry's per-fu composition reproduces exactly
    target_atoms_of_el.
    """
    target_comp = target.composition_per_formula_unit()
    members = []
    for el, n_target in target_comp.items():
        if el not in element_entries:
            raise KeyError(f"No elemental reference supplied for element '{el}' (target {target.compound_id!r})")
        ref = element_entries[el]
        ref_comp = ref.composition_per_formula_unit()
        if set(ref_comp) != {el}:
            raise ValueError(
                f"Elemental reference for '{el}' ({ref.compound_id!r}) is not a pure "
                f"single-element compound (composition: {ref_comp!r})"
            )
        atoms_per_fu = ref_comp[el]
        members.append(ReactionMember(compound_id=ref.compound_id, coefficient=n_target / atoms_per_fu))
    return members
