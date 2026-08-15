"""Tests for reaction_analysis.delta, on the fixture compounds parsed via
parse_lobster.py. Expected numbers below are hand-computed from the raw
ICOHP values in examples/dataset/{compound_A,compound_B,compound_C_disconnected}
and tests/fixtures/{compound_K,compound_Na2K} -- see tests/fixtures/README.md.

Sign convention: delta = products - reactants (see delta.py docstring;
opposite of reaction_icohp.py's convention -- not interchangeable).
"""

import math
import unittest
from pathlib import Path

from reaction_analysis.delta import compute_delta
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


class TestCase1DecompositionToElements(unittest.TestCase):
    """Na2K -> 2 Na + K. Na2K sum_per_fu=-2.25 (n_bonds=6); Na sum_per_fu=
    -2.2 (n_bonds=3); K sum_per_fu=-2.8 (n_bonds=3)."""

    def setUp(self):
        self.entries = _entries()
        self.reaction = Reaction(
            reaction_id="Na2K_to_elements",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[
                ReactionMember(compound_id="Na", coefficient=2),
                ReactionMember(compound_id="K", coefficient=1),
            ],
        )

    def test_delta_per_formula_unit(self):
        [result] = compute_delta(self.reaction, self.entries)
        self.assertIsNone(result.error)
        self.assertAlmostEqual(result.delta_per_formula_unit_eV, -4.95, places=6)

    def test_delta_per_atom(self):
        [result] = compute_delta(self.reaction, self.entries)
        self.assertAlmostEqual(result.delta_per_atom_eV, -1.65, places=6)

    def test_delta_per_bond_not_conservative(self):
        [result] = compute_delta(self.reaction, self.entries)
        self.assertAlmostEqual(result.delta_per_bond_eV, -0.425, places=6)
        self.assertFalse(result.delta_per_bond_conservative)

    def test_dimensional_consistency_atom_vs_formula_unit(self):
        """delta_per_atom_eV * N_atoms_transferred == delta_per_formula_unit_eV
        exactly, by construction (delta.py's own definition of (b) as
        (a) / N_atoms) -- see the section-7-style consistency check the
        original request asked for."""
        [result] = compute_delta(self.reaction, self.entries)
        n_atoms = 3  # Na2K has 2 Na + 1 K per formula unit
        self.assertAlmostEqual(result.delta_per_atom_eV * n_atoms, result.delta_per_formula_unit_eV, places=6)


class TestCase2DecompositionToCompoundAndElements(unittest.TestCase):
    """Na2K -> NaK + Na. NaK sum_per_fu=-8.0 (n_bonds=5)."""

    def setUp(self):
        self.entries = _entries()
        self.reaction = Reaction(
            reaction_id="Na2K_to_NaK_plus_Na",
            type="decomposition_to_compound_and_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[
                ReactionMember(compound_id="NaK", coefficient=1),
                ReactionMember(compound_id="Na", coefficient=1),
            ],
        )

    def test_deltas(self):
        [result] = compute_delta(self.reaction, self.entries)
        self.assertIsNone(result.error)
        self.assertAlmostEqual(result.delta_per_formula_unit_eV, -7.95, places=6)
        self.assertAlmostEqual(result.delta_per_atom_eV, -2.65, places=6)
        self.assertAlmostEqual(result.delta_per_bond_eV, -0.9, places=6)


class TestCase3PolymorphTransition(unittest.TestCase):
    """compound_A (Na) -> compound_C_disconnected (Na). A sum_per_fu=-2.2
    (n_bonds=3); C sum_per_fu=-2.0 (n_bonds=2)."""

    def setUp(self):
        self.entries = _entries()
        self.reaction = Reaction(
            reaction_id="Na_polyA_to_polyC",
            type="polymorph_transition",
            reactants=[ReactionMember(compound_id="Na_polyA", coefficient=1)],
            products=[ReactionMember(compound_id="Na_polyC", coefficient=1)],
        )

    def test_deltas(self):
        [result] = compute_delta(self.reaction, self.entries)
        self.assertIsNone(result.error)
        self.assertAlmostEqual(result.delta_per_formula_unit_eV, 0.2, places=6)
        self.assertAlmostEqual(result.delta_per_atom_eV, 0.2, places=6)
        self.assertAlmostEqual(result.delta_per_bond_eV, -1.0 - (-2.2 / 3), places=6)


class TestUnbalancedReactionNeverRaises(unittest.TestCase):
    def test_unbalanced_reaction_returns_nan_and_error_instead_of_raising(self):
        entries = _entries()
        reaction = Reaction(
            reaction_id="bad",
            type="decomposition_to_compound_and_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[ReactionMember(compound_id="NaK", coefficient=1)],  # missing 1 Na
        )
        [result] = compute_delta(reaction, entries)  # must not raise
        self.assertIsNotNone(result.error)
        self.assertIn("Na", result.error)
        self.assertTrue(math.isnan(result.delta_per_formula_unit_eV))
        self.assertTrue(math.isnan(result.delta_per_atom_eV))
        self.assertTrue(math.isnan(result.delta_per_bond_eV))


class TestBondPairRestriction(unittest.TestCase):
    def setUp(self):
        self.entries = _entries()

    def test_bond_pair_missing_on_one_product_member(self):
        reaction = Reaction(
            reaction_id="Na2K_to_elements_KK",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[
                ReactionMember(compound_id="Na", coefficient=2),
                ReactionMember(compound_id="K", coefficient=1),
            ],
            bond_pair="K-K",
        )
        [result] = compute_delta(reaction, self.entries)
        # products side: Na has no K-K bonds -> unavailable -> NaN, warned
        self.assertTrue(math.isnan(result.delta_per_formula_unit_eV))
        self.assertTrue(result.warnings)

    def test_missing_bond_pair_warns_without_crashing(self):
        # K has no "Na-Na" bond type at all.
        reaction = Reaction(
            reaction_id="Na2K_to_elements_NaNa",
            type="decomposition_to_elements",
            reactants=[ReactionMember(compound_id="Na2K", coefficient=1)],
            products=[
                ReactionMember(compound_id="Na", coefficient=2),
                ReactionMember(compound_id="K", coefficient=1),
            ],
            bond_pair="Na-Na",
        )
        [result] = compute_delta(reaction, self.entries)  # must not raise
        self.assertIsNone(result.error)
        self.assertTrue(math.isnan(result.delta_per_formula_unit_eV))
        self.assertTrue(any("K" in w and "Na-Na" in w for w in result.warnings))


class TestMultiMetricRequest(unittest.TestCase):
    def test_icohp_icobi_returns_two_results(self):
        entries = _entries()
        reaction = Reaction(
            reaction_id="Na_polyA_to_polyC_both_metrics",
            type="polymorph_transition",
            reactants=[ReactionMember(compound_id="Na_polyA", coefficient=1)],
            products=[ReactionMember(compound_id="Na_polyC", coefficient=1)],
            metric="icohp,icobi",
        )
        results = compute_delta(reaction, entries)
        self.assertEqual([r.metric for r in results], ["icohp", "icobi"])
        # Neither fixture has ICOBILIST.lobster -> icobi metric errors out per-compound, not a crash.
        icobi_result = results[1]
        self.assertIsNotNone(icobi_result.error)


if __name__ == "__main__":
    unittest.main()
