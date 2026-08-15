"""Tests for reaction_analysis.balance, on the fixture compounds parsed
via parse_lobster.py (examples/dataset/ + tests/fixtures/, see that
directory's README for why each one is there)."""

import unittest
from pathlib import Path

from reaction_analysis.balance import (
    ReactionBalanceError,
    check_balance,
    derive_element_coefficients,
)
from reaction_analysis.parse_lobster import parse_compound_entry
from reaction_analysis.schema import Reaction, ReactionMember

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "dataset"
FIXTURES_DIR = Path(__file__).parent / "fixtures"


def _entries():
    return {
        "Na": parse_compound_entry(EXAMPLES_DIR / "compound_A", role="element", compound_id="Na"),
        "K": parse_compound_entry(FIXTURES_DIR / "compound_K", role="element", compound_id="K"),
        "NaK": parse_compound_entry(EXAMPLES_DIR / "compound_B", role="hull_neighbor", compound_id="NaK"),
        "Na2K": parse_compound_entry(FIXTURES_DIR / "compound_Na2K", role="target", compound_id="Na2K"),
        "Na_polyA": parse_compound_entry(EXAMPLES_DIR / "compound_A", role="polymorph", compound_id="Na_polyA"),
        "Na_polyC": parse_compound_entry(
            EXAMPLES_DIR / "compound_C_disconnected", role="polymorph", compound_id="Na_polyC"
        ),
    }


class TestCheckBalance(unittest.TestCase):
    def setUp(self):
        self.entries = _entries()

    def test_decomposition_to_elements_balances(self):
        reaction = Reaction(
            reaction_id="Na2K_to_elements",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[
                ReactionMember(compound_id="Na", coefficient=2),
                ReactionMember(compound_id="K", coefficient=1),
            ],
        )
        check_balance(reaction, self.entries)  # should not raise

    def test_decomposition_to_compound_and_elements_balances(self):
        reaction = Reaction(
            reaction_id="Na2K_to_NaK_plus_Na",
            type="decomposition_to_compound_and_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[
                ReactionMember(compound_id="NaK", coefficient=1),
                ReactionMember(compound_id="Na", coefficient=1),
            ],
        )
        check_balance(reaction, self.entries)  # should not raise

    def test_polymorph_transition_balances(self):
        reaction = Reaction(
            reaction_id="Na_polyA_to_polyC",
            type="polymorph_transition",
            reactants=[ReactionMember(compound_id="Na_polyA", coefficient=1)],
            products=[ReactionMember(compound_id="Na_polyC", coefficient=1)],
        )
        check_balance(reaction, self.entries)  # should not raise

    def test_unbalanced_reaction_names_the_offending_element(self):
        reaction = Reaction(
            reaction_id="Na2K_missing_one_Na",
            type="decomposition_to_compound_and_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[ReactionMember(compound_id="NaK", coefficient=1)],  # missing 1 Na
        )
        with self.assertRaises(ReactionBalanceError) as ctx:
            check_balance(reaction, self.entries)
        self.assertIn("Na", str(ctx.exception))

    def test_unbalanced_by_more_than_one_element(self):
        reaction = Reaction(
            reaction_id="bogus",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[ReactionMember(compound_id="Na", coefficient=1)],  # missing K entirely, wrong Na too
        )
        with self.assertRaises(ReactionBalanceError) as ctx:
            check_balance(reaction, self.entries)
        message = str(ctx.exception)
        self.assertIn("Na", message)
        self.assertIn("K", message)


class TestDeriveElementCoefficients(unittest.TestCase):
    def setUp(self):
        self.entries = _entries()

    def test_derives_correct_coefficients_for_Na2K(self):
        members = derive_element_coefficients(
            self.entries["Na2K"], {"Na": self.entries["Na"], "K": self.entries["K"]}
        )
        by_id = {m.compound_id: m.coefficient for m in members}
        self.assertAlmostEqual(by_id["Na"], 2.0)
        self.assertAlmostEqual(by_id["K"], 1.0)

    def test_derived_coefficients_actually_balance(self):
        members = derive_element_coefficients(
            self.entries["Na2K"], {"Na": self.entries["Na"], "K": self.entries["K"]}
        )
        reaction = Reaction(
            reaction_id="derived",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=members,
        )
        check_balance(reaction, self.entries)  # should not raise

    def test_missing_element_reference_raises(self):
        with self.assertRaises(KeyError):
            derive_element_coefficients(self.entries["Na2K"], {"Na": self.entries["Na"]})  # no K supplied

    def test_non_pure_element_reference_rejected(self):
        with self.assertRaises(ValueError):
            derive_element_coefficients(
                self.entries["Na2K"], {"Na": self.entries["NaK"], "K": self.entries["K"]}
            )


if __name__ == "__main__":
    unittest.main()
