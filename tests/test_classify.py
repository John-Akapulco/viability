"""Tests for reaction_analysis.classify -- pure delta_icohp/delta_energy
inputs, no LOBSTER/pymatgen dependency."""

import unittest

from reaction_analysis.classify import (
    BondingLabel,
    ViabilityLabel,
    classify_bonding,
    classify_viability,
)


class TestClassifyBonding(unittest.TestCase):
    def test_positive_is_endobondic(self):
        self.assertEqual(classify_bonding(1.5), BondingLabel.ENDOBONDIC)

    def test_negative_is_exobondic(self):
        self.assertEqual(classify_bonding(-1.5), BondingLabel.EXOBONDIC)

    def test_zero_is_endobondic_by_convention(self):
        self.assertEqual(classify_bonding(0.0), BondingLabel.ENDOBONDIC)


class TestClassifyViability(unittest.TestCase):
    def test_non_negative_delta_energy_is_stable_on_hull_regardless_of_bonding(self):
        result = classify_viability(delta_energy=0.2, delta_icohp=-5.0)
        self.assertEqual(result.label, ViabilityLabel.STABLE_ON_HULL)
        self.assertEqual(result.warnings, [])

    def test_exothermic_and_endobondic_is_metastable_viable(self):
        result = classify_viability(delta_energy=-0.3, delta_icohp=2.0)
        self.assertEqual(result.label, ViabilityLabel.METASTABLE_VIABLE)
        self.assertEqual(result.bonding_label, BondingLabel.ENDOBONDIC)
        self.assertEqual(result.warnings, [])

    def test_exothermic_and_exobondic_is_unstable_nonexistent_with_kinetics_caveat(self):
        result = classify_viability(delta_energy=-0.3, delta_icohp=-2.0)
        self.assertEqual(result.label, ViabilityLabel.UNSTABLE_NONEXISTENT)
        self.assertEqual(result.bonding_label, BondingLabel.EXOBONDIC)
        self.assertEqual(len(result.warnings), 1)
        self.assertIn("never proves the compound cannot exist", result.warnings[0])

    def test_ambiguous_requires_both_reference_params_explicitly(self):
        # Small |delta_icohp| relative to a caller-supplied reference, but
        # no ambiguous_ratio_threshold given -> must NOT auto-trigger.
        result = classify_viability(
            delta_energy=-0.3, delta_icohp=-0.01, exobondic_reference_magnitude=2.0
        )
        self.assertEqual(result.label, ViabilityLabel.UNSTABLE_NONEXISTENT)

        # Threshold given but no reference magnitude -> still must NOT
        # auto-trigger.
        result2 = classify_viability(
            delta_energy=-0.3, delta_icohp=-0.01, ambiguous_ratio_threshold=0.5
        )
        self.assertEqual(result2.label, ViabilityLabel.UNSTABLE_NONEXISTENT)

    def test_ambiguous_triggers_when_both_params_given_and_magnitude_is_small(self):
        result = classify_viability(
            delta_energy=-0.3,
            delta_icohp=-0.01,
            ambiguous_ratio_threshold=0.5,
            exobondic_reference_magnitude=2.0,
        )
        self.assertEqual(result.label, ViabilityLabel.AMBIGUOUS_CHECK_KINETICS)
        self.assertEqual(len(result.warnings), 2)

    def test_ambiguous_does_not_trigger_when_magnitude_is_not_small(self):
        result = classify_viability(
            delta_energy=-0.3,
            delta_icohp=-5.0,
            ambiguous_ratio_threshold=0.5,
            exobondic_reference_magnitude=2.0,
        )
        self.assertEqual(result.label, ViabilityLabel.UNSTABLE_NONEXISTENT)

    def test_mn2o7_does_not_auto_classify_as_confidently_nonexistent(self):
        # Mn2O7 -> 2 MnO2 + 3/2 O2 is the manuscript's own guard-rail case
        # (exobondic, delta_icohp = -186 kJ/mol, yet Mn2O7 is a real,
        # observed compound that just decomposes slowly). This module
        # hardcodes nothing about Mn2O7 -- the requirement being tested is
        # structural: ANY exobondic call must carry the kinetics caveat,
        # so a caller can never mistake UNSTABLE_NONEXISTENT for a
        # confident, unqualified verdict.
        from reaction_analysis.units import kj_per_mol_to_ev

        delta_icohp_eV = kj_per_mol_to_ev(-186.0)
        result = classify_viability(delta_energy=-0.05, delta_icohp=delta_icohp_eV)
        self.assertEqual(result.label, ViabilityLabel.UNSTABLE_NONEXISTENT)
        self.assertTrue(
            any("never proves the compound cannot exist" in w for w in result.warnings),
            "UNSTABLE_NONEXISTENT must always carry the kinetics caveat warning",
        )


if __name__ == "__main__":
    unittest.main()
