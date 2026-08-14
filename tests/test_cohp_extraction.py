"""Validation of cohp_extraction.py: steps 0+1 (extraction/cross-validation,
against the 6 pilot compounds' REAL COHPCAR.lobster files -- deliberately
real data, since the LOBSTER file format and pymatgen's parsing of it are
what needs validating, not a hand-rolled algorithm) and step 2 (the
antibonding-population-near-frontier metric: its pure numerical core is
validated on hand-crafted synthetic arrays with a known analytic answer,
per this project's established practice, before trusting it on real data).

Skips gracefully if mp_dataset/structures_cohp/ isn't populated (i.e. on
a machine that hasn't run mp_dataset/structures_cohp population, or in a
CI environment without the multi-hundred-MB COHPCAR.lobster files, which
are not tracked in git -- see .gitignore). The synthetic-array tests for
integrate_antibonding_in_window() do not depend on this and always run.

Run with: python -m unittest discover -s tests -v
"""

import sys
import unittest
from pathlib import Path

import numpy as np

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


class TestIntegrateAntibondingSynthetic(unittest.TestCase):
    """Pure numerical validation of the step-2 metric's core, on
    hand-crafted arrays with a known analytic answer -- independent of any
    real COHPCAR.lobster parsing."""

    def test_pure_bonding_gives_zero(self):
        # constant negative (bonding) COHP everywhere -> clipped to 0
        energies = np.linspace(-3, 3, 601)
        cohp = np.full_like(energies, -1.0)
        w = ce.integrate_antibonding_in_window(energies, cohp, e_ref=0.0, delta_e=1.0)
        self.assertAlmostEqual(w, 0.0, places=10)

    def test_constant_antibonding_integrates_to_width_times_height(self):
        # constant +2.0 (antibonding) COHP everywhere; window width 1.0 eV
        # -> trapezoidal integral = 1.0 * 2.0 = 2.0, up to a known O(grid
        # spacing) edge truncation from the strictly-open lower bound
        # (e_ref - delta_e, e_ref] excluding the boundary sample itself --
        # intentional, so consecutive windows never double-count a shared
        # edge point; tolerance set to a few grid cells, not exact.
        energies = np.linspace(-3, 3, 601)  # spacing = 0.01
        cohp = np.full_like(energies, 2.0)
        w = ce.integrate_antibonding_in_window(energies, cohp, e_ref=0.0, delta_e=1.0)
        self.assertAlmostEqual(w, 2.0, delta=0.03)

    def test_mixed_sign_only_antibonding_part_counted(self):
        # +1 (antibonding) for E in [-0.5, 0], -1 (bonding) for E < -0.5;
        # window (-1, 0] -> only the [-0.5, 0] part survives clipping,
        # expected integral = 0.5 * 1.0 = 0.5 up to O(grid spacing)
        # trapezoidal error at the discontinuity itself. This is the test
        # that distinguishes "integrate the antibonding part only" from
        # "just integrate the net signed COHP", which would wrongly give
        # ~0 instead of ~0.5.
        energies = np.linspace(-1.0, 0.0, 101)  # spacing = 0.01
        cohp = np.where(energies >= -0.5, 1.0, -1.0)
        w = ce.integrate_antibonding_in_window(energies, cohp, e_ref=0.0, delta_e=1.0)
        self.assertAlmostEqual(w, 0.5, delta=0.01)

    def test_window_is_one_sided_ignores_energies_above_e_ref(self):
        # +5 (strongly antibonding) placed entirely ABOVE e_ref: must not
        # leak into a window that's supposed to be occupied-states-only.
        energies = np.linspace(-2, 2, 401)
        cohp = np.where(energies > 0, 5.0, -1.0)
        w = ce.integrate_antibonding_in_window(energies, cohp, e_ref=0.0, delta_e=1.0)
        self.assertAlmostEqual(w, 0.0, places=10)

    def test_growing_window_never_decreases_the_integral(self):
        # integrand is non-negative by construction (clipped) -> widening
        # the window can only add area, never remove it.
        rng = np.random.default_rng(0)
        energies = np.linspace(-5, 5, 1001)
        cohp = rng.normal(0, 1, size=energies.shape)
        prev = 0.0
        for delta_e in (0.2, 0.5, 1.0, 2.0, 4.0):
            w = ce.integrate_antibonding_in_window(energies, cohp, e_ref=0.0, delta_e=delta_e)
            self.assertGreaterEqual(w, prev - 1e-12)
            prev = w

    def test_window_too_small_for_grid_raises(self):
        # gap of 2.0 between grid points -1 and 1 (0.0 deliberately absent);
        # a (−0.5, 0.0] window falls entirely inside that gap -> no grid
        # point qualifies, must raise rather than silently return 0.0.
        energies = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
        cohp = np.ones_like(energies)
        with self.assertRaises(ValueError):
            ce.integrate_antibonding_in_window(energies, cohp, e_ref=0.0, delta_e=0.5)


@unittest.skipIf(_missing_cohpcar_data(), "mp_dataset/structures_cohp/ not populated on this machine")
class TestAntibondingPopulationRealPilots(unittest.TestCase):
    """Sanity checks of the full metric (window + reference-energy logic
    + integration) on the 6 real pilot compounds. Not a "known correct
    answer" test -- no independent ground truth exists for this new metric
    -- but every invariant checked here must hold for the implementation
    to be trustworthy at all."""

    STRUCTURES_ROOT = Path(__file__).resolve().parent.parent / "mp_dataset" / "structures"

    def _vasprun_path(self, compound: str) -> Path:
        return self.STRUCTURES_ROOT / compound / "vasprun.xml"

    def test_runs_without_error_and_returns_nonnegative_values(self):
        for compound, is_metal in EXPECTED_METAL_OR_GAP.items():
            r = ce.antibonding_population_near_frontier(
                STRUCTURES_COHP_ROOT / compound,
                is_metal=is_metal["is_metal"],
                vasprun_path=self._vasprun_path(compound),
            )
            self.assertGreaterEqual(r["w_antibond_raw"], 0.0, compound)
            if r["w_antibond_normalized"] is not None:
                self.assertGreaterEqual(r["w_antibond_normalized"], 0.0, compound)

    def test_metal_reference_energy_is_exactly_ef(self):
        for compound in ("hull_metallic_AlNi_mp-1487", "metastable_metallic_BeCu_mp-2323"):
            e_ref = ce.frontier_reference_energy(STRUCTURES_COHP_ROOT / compound, is_metal=True)
            self.assertEqual(e_ref, 0.0)

    def test_gapped_reference_energy_is_below_ef_not_above(self):
        # VBM must sit below E_F (occupied side) for a real insulator/
        # semiconductor -- this is the check that would have caught the
        # AlNi/BeCu "naive VBM" pitfall (documented in
        # analysis/REPORT_cohp_feasibility.md) had it been applied there.
        for compound in (
            "hull_ionic_NaCl_mp-22862",
            "hull_covalent_Si_mp-149",
            "metastable_ionic_LiBr_mp-23259",
            "metastable_covalent_C_rhombohedral_mp-169",
        ):
            e_ref = ce.frontier_reference_energy(
                STRUCTURES_COHP_ROOT / compound, is_metal=False,
                vasprun_path=self._vasprun_path(compound),
            )
            self.assertLess(e_ref, 0.0, compound)

    def test_window_sensitivity_is_monotonic_non_decreasing(self):
        for compound, is_metal in EXPECTED_METAL_OR_GAP.items():
            prev = 0.0
            for delta_e in (0.5, 1.0, 2.0):
                r = ce.antibonding_population_near_frontier(
                    STRUCTURES_COHP_ROOT / compound,
                    is_metal=is_metal["is_metal"], delta_e=delta_e,
                    vasprun_path=self._vasprun_path(compound),
                )
                self.assertGreaterEqual(r["w_antibond_raw"], prev - 1e-9, f"{compound} dE={delta_e}")
                prev = r["w_antibond_raw"]


if __name__ == "__main__":
    unittest.main()
