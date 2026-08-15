"""Tests for reaction_analysis.schema -- shape validation only, no real
LOBSTER parsing (see test_parse_lobster.py) and no reaction computation
(see test_balance.py / test_delta.py)."""

import unittest

from pydantic import ValidationError

from reaction_analysis.schema import (
    BondTypeSummary,
    CompoundEntry,
    IcohpSummary,
    Reaction,
    ReactionMember,
    ReactionResult,
)


def _make_icohp_summary(**overrides) -> IcohpSummary:
    defaults = dict(
        sum_total_eV=-2.25,
        sum_per_atom_eV=-0.75,
        sum_per_formula_unit_eV=-2.25,
        mean_per_bond_eV=-0.375,
        n_bonds=6,
        by_bond_type={
            "Na-Na": BondTypeSummary(mean_eV=-0.45, sum_eV=-0.9, n_bonds=2),
            "K-Na": BondTypeSummary(mean_eV=-0.6, sum_eV=-1.2, n_bonds=2),
            "K-K": BondTypeSummary(mean_eV=-0.05, sum_eV=-0.05, n_bonds=1),
        },
    )
    defaults.update(overrides)
    return IcohpSummary(**defaults)


def _make_compound_entry(**overrides) -> CompoundEntry:
    defaults = dict(
        compound_id="test_Na2K",
        formula="Na2K",
        composition={"Na": 2, "K": 1},
        Z=1,
        space_group_symbol="Pm-3m",
        space_group_number=221,
        role="target",
        icohp=_make_icohp_summary(),
        source_path="tests/fixtures/compound_Na2K",
    )
    defaults.update(overrides)
    return CompoundEntry(**defaults)


class TestCompoundEntry(unittest.TestCase):
    def test_round_trips_json(self):
        entry = _make_compound_entry()
        restored = CompoundEntry.model_validate_json(entry.model_dump_json())
        self.assertEqual(entry, restored)

    def test_composition_per_formula_unit_divides_by_Z(self):
        entry = _make_compound_entry(Z=2, composition={"Na": 4, "K": 2})
        self.assertEqual(entry.composition_per_formula_unit(), {"Na": 2.0, "K": 1.0})

    def test_role_must_be_a_known_literal(self):
        with self.assertRaises(ValidationError):
            _make_compound_entry(role="not_a_real_role")

    def test_Z_must_be_positive(self):
        with self.assertRaises(ValidationError):
            _make_compound_entry(Z=0)

    def test_icobi_optional_defaults_to_none(self):
        entry = _make_compound_entry()
        self.assertIsNone(entry.icobi)


class TestReactionShape(unittest.TestCase):
    def test_reaction_member_coefficient_must_be_positive(self):
        with self.assertRaises(ValidationError):
            ReactionMember(compound_id="x", coefficient=0)
        with self.assertRaises(ValidationError):
            ReactionMember(compound_id="x", coefficient=-1)

    def test_reaction_type_must_be_a_known_literal(self):
        with self.assertRaises(ValidationError):
            Reaction(
                reaction_id="r1",
                type="not_a_real_type",
                reactants=[ReactionMember(compound_id="a", coefficient=1)],
                products=[ReactionMember(compound_id="b", coefficient=1)],
            )

    def test_reaction_metric_defaults_to_icohp(self):
        reaction = Reaction(
            reaction_id="r1",
            type="polymorph_transition",
            reactants=[ReactionMember(compound_id="a", coefficient=1)],
            products=[ReactionMember(compound_id="b", coefficient=1)],
        )
        self.assertEqual(reaction.metric, "icohp")
        self.assertIsNone(reaction.bond_pair)

    def test_reaction_result_round_trips_json(self):
        result = ReactionResult(
            reaction_id="r1", type="polymorph_transition", metric="icohp", bond_pair=None,
            delta_per_formula_unit_eV=0.2, delta_per_atom_eV=0.2, delta_per_bond_eV=-0.2667,
            delta_per_bond_conservative=False, warnings=[], error=None,
        )
        restored = ReactionResult.model_validate_json(result.model_dump_json())
        self.assertEqual(result, restored)


if __name__ == "__main__":
    unittest.main()
