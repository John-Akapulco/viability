"""Tests for reaction_analysis.nearest_neighbor -- synthetic distance
lists only, no LOBSTER/pymatgen dependency (see module docstring for why
this stays generic)."""

import unittest

from reaction_analysis.nearest_neighbor import detect_first_shell, first_shell_records


class TestDetectFirstShell(unittest.TestCase):
    def test_manuscript_style_cluster_with_growing_in_shell_gaps(self):
        # The exact failure mode the "running max gap" approach (tried
        # first, discarded) got wrong: in-shell gaps that themselves grow
        # (0.02 -> 0.05 -> 0.18) must not trigger a cut before the real
        # jump to the next shell (1.67).
        d_min, d_cut = detect_first_shell([2.58, 2.60, 2.65, 2.83, 4.50, 4.60])
        self.assertEqual(d_min, 2.58)
        self.assertEqual(d_cut, 2.83)

    def test_no_gap_at_all_returns_whole_list_as_one_shell(self):
        d_min, d_cut = detect_first_shell([2.0, 2.1, 2.2, 2.3])
        self.assertEqual((d_min, d_cut), (2.0, 2.3))

    def test_identical_distances_returns_whole_list_as_one_shell(self):
        d_min, d_cut = detect_first_shell([3.0, 3.0, 3.0])
        self.assertEqual((d_min, d_cut), (3.0, 3.0))

    def test_single_distance(self):
        self.assertEqual(detect_first_shell([2.0]), (2.0, 2.0))

    def test_two_distances_never_cut_insufficient_data(self):
        # Fewer than 2 gaps: no baseline to call a gap an outlier against
        # -- conservatively kept as one shell, even when clearly far apart.
        self.assertEqual(detect_first_shell([2.0, 4.0]), (2.0, 4.0))

    def test_three_distances_with_a_clear_outlier_gap_does_cut(self):
        d_min, d_cut = detect_first_shell([2.0, 2.02, 2.04, 5.0])
        self.assertEqual((d_min, d_cut), (2.0, 2.04))

    def test_empty_list_raises(self):
        with self.assertRaises(ValueError):
            detect_first_shell([])

    def test_gap_ratio_is_configurable(self):
        # A gap that qualifies at the default ratio can be suppressed by
        # raising gap_ratio -- confirms the threshold is a real knob, not
        # a hardcoded distance.
        distances = [2.0, 2.02, 2.04, 2.5]
        _, d_cut_default = detect_first_shell(distances, gap_ratio=2.5)
        _, d_cut_strict = detect_first_shell(distances, gap_ratio=100.0)
        self.assertEqual(d_cut_default, 2.04)
        self.assertEqual(d_cut_strict, 2.5)


class TestFirstShellRecords(unittest.TestCase):
    def test_groups_by_pair_key_and_filters_independently(self):
        records = [
            {"pair": "A-B", "d": 2.0}, {"pair": "A-B", "d": 2.02}, {"pair": "A-B", "d": 2.04},
            {"pair": "A-B", "d": 5.0},
            {"pair": "C-D", "d": 3.0}, {"pair": "C-D", "d": 3.1}, {"pair": "C-D", "d": 3.2},
        ]
        kept = first_shell_records(records, pair_key=lambda r: r["pair"], distance=lambda r: r["d"])
        kept_ab = sorted(r["d"] for r in kept if r["pair"] == "A-B")
        kept_cd = sorted(r["d"] for r in kept if r["pair"] == "C-D")
        self.assertEqual(kept_ab, [2.0, 2.02, 2.04])
        self.assertEqual(kept_cd, [3.0, 3.1, 3.2])

    def test_empty_input_returns_empty(self):
        self.assertEqual(first_shell_records([], pair_key=lambda r: r, distance=lambda r: 0.0), [])


if __name__ == "__main__":
    unittest.main()
