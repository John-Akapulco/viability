"""Automatic nearest-neighbor (first coordination shell) detection, used to
restrict ICOHPLIST.lobster bond summation to first-shell bonds only before
summing per species-pair -- the convention Reitz & Dronskowski
(ic-2026-04181q, "Calculus of Bonding Energetics" section) use for their
endobondic/exobondic ICOHP totals, and a deliberate DEPARTURE from every
other summation in this project (percolation_path.py, reaction_icohp.py,
and this package's own existing sum_total_eV), which are unfiltered by
design (see parse_lobster.py's module docstring for why unfiltered is
correct there). See parse_lobster.py's `bond_filter` parameter for where
these two conventions coexist.

No absolute distance cutoff is hardcoded anywhere in this module -- the
shell boundary is always the first point where the sorted distance list
jumps by more than `gap_ratio` times the largest gap seen so far within
the putative shell: a relative, self-calibrating rule matching the
manuscript's "first gap in the bond-distance spectrum" description without
assuming any particular bond length range or chemistry.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Callable, Hashable, TypeVar

# Numerical floor (Angstrom), not a chemical cutoff: guards against a
# degenerate cut when the first two distances are floating-point-identical
# (gap ~0), which would otherwise make gap_ratio * 0 == 0 and trigger a
# false cut on the very next, arbitrarily small, nonzero gap.
_NUMERICAL_FLOOR_ANGSTROM = 1e-3

T = TypeVar("T")


def detect_first_shell(distances: list[float], gap_ratio: float = 2.5) -> tuple[float, float]:
    """Given every bond distance observed for one species pair (e.g. every
    Pb-N distance in a cell), return (d_min, d_cut): distances in
    [d_min, d_cut] belong to the first coordination shell.

    Algorithm: sort ascending and look at the consecutive-distance gaps.
    The shell boundary is placed at the single LARGEST gap, provided it is
    an outlier -- at least `gap_ratio` times the mean of every OTHER gap
    in the list -- otherwise every distance is called one shell (no
    boundary was warranted).

    Deliberately not a naive left-to-right "first gap bigger than the
    running max" walk: that approach was tried first and produces false
    cuts whenever in-shell gaps merely grow monotonically (e.g. distances
    2.58, 2.60, 2.65, 2.83 -- a real single first-shell cluster where each
    successive gap is itself bigger than the last, 0.02 -> 0.05 -> 0.18 --
    a running-max comparison cuts inside this cluster; comparing the one
    largest gap [here, the jump to the next real shell] against the mean
    of the rest does not). This assumes a single dominant shell-defining
    gap, which is the case this module is built for (finding the first
    shell boundary, not segmenting every shell) -- it is not a general
    multi-shell clustering algorithm.

    `gap_ratio` is the one tunable knob (default 2.5) -- never an absolute
    Angstrom cutoff; it is always applied relative to the list's own other
    gaps.
    """
    if not distances:
        raise ValueError("detect_first_shell() called with an empty distance list")
    d = sorted(distances)
    if len(d) < 3:
        # Fewer than 2 gaps: no "mean of the other gaps" baseline can be
        # formed, so there is no basis to call the lone gap (if any) an
        # outlier -- conservatively treat every distance as one shell
        # rather than guess.
        return (d[0], d[-1])

    gaps = [d[i + 1] - d[i] for i in range(len(d) - 1)]
    max_gap = max(gaps)
    idx = gaps.index(max_gap)
    other_gaps = gaps[:idx] + gaps[idx + 1:]
    baseline = max(sum(other_gaps) / len(other_gaps), _NUMERICAL_FLOOR_ANGSTROM)
    if max_gap >= gap_ratio * baseline:
        return (d[0], d[idx])
    return (d[0], d[-1])


def first_shell_records(
    records: list[T],
    *,
    pair_key: Callable[[T], Hashable],
    distance: Callable[[T], float],
    gap_ratio: float = 2.5,
) -> list[T]:
    """Group `records` (any dict/object list -- e.g. parse_lobster.py's raw
    bond records) by `pair_key(record)` (e.g. the sorted species-pair
    string "Pb-N"), apply detect_first_shell() to each group's distances,
    and return only the records within their own pair's first shell.

    Grouping key and distance extraction are both caller-supplied so this
    module stays generic (no LOBSTER/pymatgen-specific field names) and
    independently testable on plain synthetic records.
    """
    by_pair: dict[Hashable, list[T]] = defaultdict(list)
    for r in records:
        by_pair[pair_key(r)].append(r)

    kept: list[T] = []
    for pair, group in by_pair.items():
        d_min, d_cut = detect_first_shell([distance(r) for r in group], gap_ratio=gap_ratio)
        kept.extend(r for r in group if d_min - 1e-9 <= distance(r) <= d_cut + 1e-9)
    return kept
