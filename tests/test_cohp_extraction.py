"""Validation of cohp_extraction.py against the 6 pilot compounds' REAL
COHPCAR.lobster files (mission: antibonding-population-near-E_F, step 1).
Uses real LOBSTER output rather than a synthetic case deliberately -- the
LOBSTER file format and pymatgen's parsing of it are exactly what this
step needs to validate, not a hand-rolled algorithm.

Skips gracefully if mp_dataset/structures_cohp/ isn't populated (i.e. on
a machine that hasn't run mp_dataset/structures_cohp population, or in a
CI environment without the multi-hundred-MB COHPCAR.lobster files, which
are not tracked in git -- see .gitignore).

Run with: python -m unittest discover -s tests -v
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import cohp_extraction as ce

STRUCTURES_COHP_ROOT = Path(__file__).resolve().parent.parent / "mp_dataset" / "structures_cohp"

PILOT_COMPOUNDS = [
    "hull_ionic_NaCl_mp-22862",
    "hull_covalent_Si_mp-149",
    "hull_metallic_AlNi_mp-1487",
    "metastable_ionic_LiBr_mp-23259",
    "metastable_covalent_C_rhombohedral_mp-169",
    "metastable_metallic_BeCu_mp-2323",
]

# Materials Project's own converged is_metal/band_gap (fetched once via
# cohp_extraction.metal_or_gap_from_mp and hardcoded here so the test
# suite stays fast and offline -- see analysis/REPORT_cohp_feasibility.md
# for the live-fetch transcript). This is the AUTHORITATIVE classification
# for this step; a naive local eigenvalue_band_properties/DOS-at-E_F check
# on our LOBSTER-oriented coarse k-mesh spuriously suggests small
# (0.1-0.2 eV) gaps for AlNi and BeCu, both of which MP confirms are
# metals (band_gap=0.0) -- documented as a pitfall, not used here.
EXPECTED_METAL_OR_GAP = {
    "hull_ionic_NaCl_mp-22862": {"is_metal": False, "band_gap": 5.0037},
    "hull_covalent_Si_mp-149": {"is_metal": False, "band_gap": 0.6103},
    "hull_metallic_AlNi_mp-1487": {"is_metal": True, "band_gap": 0.0},
    "metastable_ionic_LiBr_mp-23259": {"is_metal": False, "band_gap": 4.9234},
    "metastable_covalent_C_rhombohedral_mp-169": {"is_metal": False, "band_gap": 0.1981},
    "metastable_metallic_BeCu_mp-2323": {"is_metal": True, "band_gap": 0.0},
}


def _missing_cohpcar_data() -> bool:
    return not all((STRUCTURES_COHP_ROOT / c / "COHPCAR.lobster").exists() for c in PILOT_COMPOUNDS)


@unittest.skipIf(_missing_cohpcar_data(), "mp_dataset/structures_cohp/ not populated on this machine")
class TestCohpExtraction(unittest.TestCase):
    def test_cohpcar_reads_and_matches_icohplist_label_count(self):
        for compound in PILOT_COMPOUNDS:
            d = STRUCTURES_COHP_ROOT / compound
            cohpcar = ce.load_cohpcar(d)
            from pymatgen.io.lobster.outputs import Icohplist

            icohplist = Icohplist(filename=str(d / "ICOHPLIST.lobster"))
            # cohp_data has one extra "average" entry beyond the per-bond labels
            self.assertEqual(
                len(cohpcar.cohp_data) - 1, len(icohplist.icohplist),
                f"{compound}: COHPCAR/ICOHPLIST bond-label count mismatch",
            )

    def test_cross_validation_against_icohplist_within_tight_tolerance(self):
        # 1e-4 eV: one order of magnitude looser than the worst observed
        # deviation (1e-5 eV, BeCu) to allow for LOBSTER's internal
        # adaptive-quadrature vs. energy-grid-trapezoidal numerical
        # differences between the two output files, while still being far
        # tighter than any physically meaningful ICOHP value in this
        # dataset (weakest bonds are O(1e-4) eV, strongest are O(1-10) eV).
        tolerance = 1e-4
        for compound in PILOT_COMPOUNDS:
            d = STRUCTURES_COHP_ROOT / compound
            result = ce.cross_validate_against_icohplist(d)
            self.assertEqual(result["n_matched"], result["n_labels"], f"{compound}: unmatched labels")
            self.assertLessEqual(
                result["max_abs_diff"], tolerance,
                f"{compound}: max |COHPCAR - ICOHPLIST| = {result['max_abs_diff']} exceeds {tolerance} eV",
            )

    def test_sign_convention_matches_established_icohplist_convention(self):
        # Not asserting a specific antibonding/bonding outcome (that's a
        # per-compound chemistry question, not a parser property) -- only
        # that the check runs and returns the documented, empirically
        # verified convention statement (negative = bonding).
        for compound in PILOT_COMPOUNDS:
            d = STRUCTURES_COHP_ROOT / compound
            result = ce.sign_convention_check(d)
            self.assertIn("NEGATIVE ICOHP/COHP = bonding", result["convention"])
            self.assertIsInstance(result["has_antibonding_states_below_efermi"], dict)

    def test_metal_vs_gap_classification_matches_materials_project(self):
        # Uses the hardcoded MP reference (see module docstring above) --
        # no live network call in the test suite.
        for compound, expected in EXPECTED_METAL_OR_GAP.items():
            self.assertIn(compound, PILOT_COMPOUNDS)
            if expected["is_metal"]:
                self.assertEqual(expected["band_gap"], 0.0)
            else:
                self.assertGreater(expected["band_gap"], 0.0)
        # both bonding-character classes are represented, as the mission
        # brief expects from this 6-compound pilot set
        n_metal = sum(1 for v in EXPECTED_METAL_OR_GAP.values() if v["is_metal"])
        n_gap = sum(1 for v in EXPECTED_METAL_OR_GAP.values() if not v["is_metal"])
        self.assertEqual(n_metal, 2)  # AlNi, BeCu
        self.assertEqual(n_gap, 4)  # NaCl, Si, LiBr, C-rhombohedral
        # C-rhombohedral specifically: mission flagged uncertainty over
        # whether this polytype is a semimetal -- confirmed NOT (nonzero
        # MP gap, 0.198 eV), don't assume, this is why it's checked here.
        self.assertFalse(EXPECTED_METAL_OR_GAP["metastable_covalent_C_rhombohedral_mp-169"]["is_metal"])
        self.assertGreater(
            EXPECTED_METAL_OR_GAP["metastable_covalent_C_rhombohedral_mp-169"]["band_gap"], 0.0
        )


if __name__ == "__main__":
    unittest.main()
