"""Tests for reaction_analysis.parse_lobster.

Not in the originally requested test-file list (test_schema.py,
test_balance.py, test_delta.py) -- added separately because the
double-counting empirical check the request asked for ("confirmer sur
ICOHPLIST.lobster réel si chaque liaison i-j est listée une seule fois ou
si les deux sens apparaissent séparément") doesn't fit naturally into any
of those three, and parsing itself deserves direct coverage regardless.
"""

import unittest
from pathlib import Path

from reaction_analysis.parse_lobster import (
    check_no_reverse_duplicates,
    parse_compound_entry,
    raw_bond_records,
)

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "dataset"
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestNoReverseDuplicates(unittest.TestCase):
    """The empirical check: LOBSTER lists each periodic bond once, not
    once per direction. Verified on every examples/dataset/ and
    tests/fixtures/ compound -- if this ever fails on real production
    data, IcohpSummary.sum_total_eV's unfiltered-sum convention (shared
    with reaction_icohp.py and percolation_path.py's edge-adding logic)
    would need revisiting."""

    def test_no_reverse_duplicates_in_any_fixture(self):
        dirs = [
            EXAMPLES_DIR / "compound_A",
            EXAMPLES_DIR / "compound_B",
            EXAMPLES_DIR / "compound_C_disconnected",
            FIXTURES_DIR / "compound_K",
            FIXTURES_DIR / "compound_Na2K",
        ]
        for d in dirs:
            with self.subTest(compound=d.name):
                records = raw_bond_records(d, "ICOHPLIST.lobster", are_cobis=False)
                self.assertEqual(check_no_reverse_duplicates(records), [])

    def test_detects_a_planted_reverse_duplicate(self):
        # Sanity check on the detector itself: compound_B's bond #2
        # (Na1, K2, (0,0,0)) reversed would be (K2, Na1, (0,0,0)) -- not
        # naturally present, so plant it and confirm it's caught.
        records = raw_bond_records(EXAMPLES_DIR / "compound_B", "ICOHPLIST.lobster", are_cobis=False)
        planted = records + [{"atom1": "K2", "atom2": "Na1", "translation": (0, 0, 0), "value": -0.5}]
        problems = check_no_reverse_duplicates(planted)
        self.assertTrue(problems, "detector should flag the planted reverse-direction duplicate")


class TestParseCompoundEntry(unittest.TestCase):
    def test_compound_A_pure_element(self):
        entry = parse_compound_entry(EXAMPLES_DIR / "compound_A", role="element")
        self.assertEqual(entry.compound_id, "compound_A")
        self.assertEqual(entry.composition, {"Na": 1})
        self.assertEqual(entry.Z, 1)
        self.assertEqual(entry.icohp.n_bonds, 3)
        self.assertAlmostEqual(entry.icohp.sum_total_eV, -2.2)
        self.assertAlmostEqual(entry.icohp.sum_per_atom_eV, -2.2)
        self.assertAlmostEqual(entry.icohp.sum_per_formula_unit_eV, -2.2)
        self.assertAlmostEqual(entry.icohp.mean_per_bond_eV, -2.2 / 3)
        self.assertIsNone(entry.icobi)
        self.assertIn("Na-Na", entry.icohp.by_bond_type)
        self.assertEqual(entry.icohp.by_bond_type["Na-Na"].n_bonds, 3)

    def test_compound_B_binary_by_bond_type(self):
        entry = parse_compound_entry(EXAMPLES_DIR / "compound_B", role="target")
        self.assertEqual(entry.composition, {"Na": 1, "K": 1})
        self.assertEqual(entry.Z, 1)
        self.assertAlmostEqual(entry.icohp.sum_total_eV, -8.0)
        keys = set(entry.icohp.by_bond_type)
        self.assertEqual(keys, {"Na-Na", "K-Na"})
        self.assertEqual(entry.icohp.by_bond_type["K-Na"].n_bonds, 2)
        self.assertAlmostEqual(entry.icohp.by_bond_type["K-Na"].sum_eV, -1.0)

    def test_compound_Na2K_multi_element(self):
        entry = parse_compound_entry(FIXTURES_DIR / "compound_Na2K", role="target")
        self.assertEqual(entry.composition, {"Na": 2, "K": 1})
        self.assertEqual(entry.Z, 1)
        self.assertAlmostEqual(entry.icohp.sum_total_eV, -2.25)
        self.assertAlmostEqual(entry.icohp.sum_per_formula_unit_eV, -2.25)
        self.assertAlmostEqual(entry.icohp.mean_per_bond_eV, -2.25 / 6)
        self.assertAlmostEqual(entry.icohp.by_bond_type["Na-Na"].sum_eV, -1.0)
        self.assertAlmostEqual(entry.icohp.by_bond_type["K-Na"].sum_eV, -1.2)
        self.assertAlmostEqual(entry.icohp.by_bond_type["K-K"].sum_eV, -0.05)

    def test_role_is_never_inferred_and_is_required(self):
        with self.assertRaises(TypeError):
            parse_compound_entry(EXAMPLES_DIR / "compound_A")  # type: ignore[call-arg]

    def test_compound_id_defaults_to_directory_name(self):
        entry = parse_compound_entry(FIXTURES_DIR / "compound_K", role="element")
        self.assertEqual(entry.compound_id, "compound_K")


class TestNearestNeighborBondFilter(unittest.TestCase):
    """compound_Na2K's Na-Na pair has distances [2.0, 4.0, 4.0] Angstrom
    (bonds #1 Na1-Na2@2.0, #4 Na1-Na1@4.0, #5 Na2-Na2@4.0 in its
    ICOHPLIST.lobster) -- a real gap, unlike every other pair in the toy
    fixtures (identical or single distances), making it the one case in
    this project's existing fixtures that actually exercises
    nearest_neighbor.py's filtering end to end through parse_lobster.py."""

    def test_default_is_unfiltered_unchanged(self):
        entry = parse_compound_entry(FIXTURES_DIR / "compound_Na2K", role="target")
        self.assertEqual(entry.icohp.n_bonds, 6)
        self.assertAlmostEqual(entry.icohp.sum_total_eV, -2.25)

    def test_nearest_neighbor_filter_drops_the_far_Na_Na_bonds(self):
        entry = parse_compound_entry(
            FIXTURES_DIR / "compound_Na2K", role="target", bond_filter="nearest_neighbor"
        )
        # Na-Na: only the 2.0 A bond (-0.80) survives; the two 4.0 A bonds
        # (-0.10 each) are cut. Na-K: both bonds at 2.5 A (identical
        # distances, no gap) survive. K-K: single 4.0 A bond, trivially
        # its own shell.
        self.assertEqual(entry.icohp.n_bonds, 4)
        self.assertAlmostEqual(entry.icohp.sum_total_eV, -2.05)
        self.assertEqual(entry.icohp.by_bond_type["Na-Na"].n_bonds, 1)
        self.assertAlmostEqual(entry.icohp.by_bond_type["Na-Na"].sum_eV, -0.80)
        self.assertEqual(entry.icohp.by_bond_type["K-Na"].n_bonds, 2)

    def test_unknown_bond_filter_raises(self):
        with self.assertRaises(ValueError):
            parse_compound_entry(
                FIXTURES_DIR / "compound_Na2K", role="target", bond_filter="bogus"
            )


if __name__ == "__main__":
    unittest.main()
