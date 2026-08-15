"""Pydantic data models for the DeltaICOHP/DeltaICOBI reaction-analysis
schema. This module defines the *shape* of the data only -- no real
LOBSTER data is parsed here (see parse_lobster.py) and no reaction is
computed here (see balance.py / delta.py).

No chemical roles, thresholds, or decision rules are hardcoded here --
CompoundEntry.role is pure caller-supplied metadata; nothing in this
package branches on its value.
"""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import BaseModel, Field


class BondTypeSummary(BaseModel):
    """Aggregate over every bond between one pair of atomic species (e.g.
    "Fe-O", element symbols only, alphabetically ordered by
    parse_lobster.py) -- for --bond-pair-style filtering in delta.py."""

    mean_eV: float
    sum_eV: float
    n_bonds: int


class IcohpSummary(BaseModel):
    """One metric's (ICOHP or ICOBI) summary over a compound's whole
    cell, in every normalization delta.py needs.

    sum_total_eV is the raw sum over every symmetry-inequivalent bond
    label in ICOHPLIST.lobster/ICOBILIST.lobster as parsed by
    parse_lobster.py -- see that module's docstring for the empirical
    check (on tests/fixtures) that LOBSTER lists each periodic bond once,
    not once per direction, which is what makes an unfiltered sum over
    all labels the correct "no double-counting" total.
    """

    sum_total_eV: float
    sum_per_atom_eV: float
    sum_per_formula_unit_eV: float
    mean_per_bond_eV: float
    n_bonds: int
    by_bond_type: dict[str, BondTypeSummary] = Field(default_factory=dict)


class CompoundEntry(BaseModel):
    """One calculated structure: a reaction target, an elemental
    reference, a hull-neighbor compound, or a polymorph.

    `role` is descriptive metadata only -- nothing in schema.py,
    balance.py, or delta.py inspects it; a Reaction's reactants/products
    lists are what actually determine a compound's role in any given
    reaction, and the same CompoundEntry can appear as a "target" in one
    Reaction and a "product" in another.
    """

    compound_id: str
    formula: str
    composition: dict[str, int]
    Z: int = Field(gt=0)
    space_group_symbol: str
    space_group_number: int
    role: Literal["target", "element", "hull_neighbor", "polymorph"]
    energy_total_eV: Optional[float] = None
    energy_per_atom_eV: Optional[float] = None
    icohp: IcohpSummary
    icobi: Optional[IcohpSummary] = None
    source_path: str

    def composition_per_formula_unit(self) -> dict[str, float]:
        """Atoms per element per formula unit (composition / Z) -- the
        quantity balance.py actually checks equality of, since
        `composition` itself is cell-size-dependent (an artifact of how
        large a supercell VASP/LOBSTER happened to run on) and not
        directly comparable across compounds without dividing out Z
        first."""
        return {el: n / self.Z for el, n in self.composition.items()}


class ReactionMember(BaseModel):
    compound_id: str
    coefficient: float = Field(gt=0)


class Reaction(BaseModel):
    reaction_id: str
    type: Literal[
        "decomposition_to_elements",
        "decomposition_to_compound_and_elements",
        "polymorph_transition",
    ]
    reactants: list[ReactionMember]
    products: list[ReactionMember]
    metric: Literal["icohp", "icobi", "icohp,icobi"] = "icohp"
    bond_pair: Optional[str] = None


class ReactionResult(BaseModel):
    reaction_id: str
    type: str
    metric: str
    bond_pair: Optional[str] = None
    delta_per_formula_unit_eV: float
    delta_per_atom_eV: float
    delta_per_bond_eV: float
    delta_per_bond_conservative: bool = False
    bonding_label: Optional[Literal["endobondic", "exobondic"]] = None
    warnings: list[str] = Field(default_factory=list)
    error: Optional[str] = None
