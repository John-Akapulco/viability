# Antibonding-population-near-frontier: window and metric definition (step 2)

**Status: defined, implemented, and validated on the 6 pilot compounds only.** Extension to the 186-compound dataset is a separate, not-yet-authorized mission, per the same validate-small-then-decide pattern used for every other descriptor in this project (network dimensionality, periodic min-cut).

## 1. Physical motivation (recap)

Hypothesis under test (Dronskowski/Deringer/Tchougréeff-style analysis, by analogy with Peierls/Jahn-Teller distortions): a crystal whose highest-energy *occupied* electronic states have significant antibonding COHP character sits close to an electronic instability — a symmetry-lowering distortion that depopulates those antibonding states (pushing them above the new gap) while stabilizing bonding combinations below can lower the total energy. Only **occupied** antibonding character matters for this argument — unoccupied antibonding states don't destabilize the ground state, they're simply empty.

This is a genuinely different question from everything computed so far in this project: `percolation_path.py`, `network_dimensionality.py`, and `periodic_mincut.py` all operate on **integrated** ICOHP/ICOBI (a single number per bond, already summed over all energies up to E_F) — none of them look at *how* COHP is distributed in energy, i.e. none can see whether the frontier states specifically are more or less antibonding than the average.

## 2. Reference energy: E_F for metals, VBM for gapped compounds

- **Metal**: E_ref = E_F. LOBSTER's own `COHPCAR.lobster` file is written with energies already shifted so E_F = 0 (confirmed empirically in step 1, not merely assumed from the pymatgen docstring) — so for metals E_ref = 0.0 directly on the file's own energy axis, no extra computation needed.
- **Gapped compound**: no states exist exactly at E_F (it sits in the gap), so E_ref = VBM instead — the top of the *occupied* manifold, the correct analogue of "the frontier" for an insulator/semiconductor in a second-order Jahn-Teller argument (states just below VBM are the ones whose antibonding character could drive a gap-opening distortion).

**Alignment, verified not assumed**: VBM comes from `pymatgen.io.vasp.Vasprun.eigenvalue_band_properties` on the *same* compound's `vasprun.xml`, which reports VBM on VASP's own absolute eigenvalue scale — a different, per-calculation-arbitrary zero than LOBSTER's E-E_F-shifted `COHPCAR.lobster` axis. Before trusting `VBM - E_F` as a meaningful shift, VASP's own reported `efermi` (from `vasprun.xml`) was checked against LOBSTER's reported `efermi` (from `COHPCAR.lobster`) for all 6 pilots:

| Compound | VASP E_F (eV) | LOBSTER E_F (eV) | diff |
|---|---:|---:|---:|
| NaCl | −0.1859 | −0.1859 | <10⁻⁴ |
| Si | 6.0139 | 6.0139 | <10⁻⁴ |
| AlNi | 10.2672 | 10.2672 | <10⁻⁴ |
| LiBr | 1.0058 | 1.0057 | <10⁻⁴ |
| C rhombohedral | 4.3681 | 4.3681 | <10⁻⁴ |
| BeCu | 7.8147 | 7.8147 | <10⁻⁴ |

Agreement to <10⁻⁴ eV across the board confirms the two energy scales are safely interchangeable, so `E_ref = VBM_absolute − E_F,LOBSTER_absolute` correctly lands VBM on `COHPCAR.lobster`'s own E-E_F axis. This check is now a hard runtime assertion in `cohp_extraction.frontier_reference_energy()` (raises if the two E_F values disagree by more than 0.01 eV) — should this ever fail on a new compound (e.g. one where LOBSTER's basis-fitted Fermi level diverges more from VASP's), that compound must be investigated before its metric is trusted, not silently computed anyway.

**Known pitfall avoided**: `is_metal` for this classification must come from Materials Project's own converged calculation (`cohp_extraction.metal_or_gap_from_mp`), not from a naive local check on this project's LOBSTER-oriented coarse k-mesh — confirmed in step 1 that such a local check spuriously suggests small gaps (0.14–0.17 eV) for AlNi and BeCu, both unambiguous metals. Using the wrong reference (a spurious local "VBM" that, for these two, actually sits *above* E_F — confirmed: +0.12 eV for AlNi, +0.02 eV for BeCu) would silently window the wrong energy region entirely.

## 3. Window: one-sided, `(E_ref − ΔE, E_ref]`

Not a window straddling E_ref: only occupied states (at or below the reference) can be antibonding-destabilizing in the ground state, per §1. Default **ΔE = 1.0 eV**, with sensitivity checked at ΔE ∈ {0.5, 1.0, 2.0} eV (same practice as the θ sensitivity check for network dimensionality). 1.0 eV is a round, defensible choice: comfortably wider than the calculation's own Gaussian smearing (`SIGMA = 0.05` eV) and the `COHPCAR.lobster` energy-grid spacing (~0.05 eV, so ~20 grid points fall in a 1 eV window) so the result isn't dominated by smearing/grid noise, while still being narrow enough to probe the *frontier* specifically rather than recovering the whole-compound average that `icohp_sum`/`icohp_mean` already capture.

## 4. Metric: integrated antibonding (positive) COHP in the window

$$W_\text{antibond} = \int_{E_\text{ref}-\Delta E}^{E_\text{ref}} \max(\text{COHP}(E),\, 0)\, dE$$

on LOBSTER's own **"average"** trace (summed over all bond pairs — the natural whole-compound scalar, matching how `icohp_sum`/`icohp_mean` are already whole-compound aggregates elsewhere in this project), trapezoidal integration on the native `COHPCAR.lobster` energy grid. Only the antibonding (positive) part is kept — a net/signed integral over the window would let bonding and antibonding character in the same window cancel, hiding exactly the "large simultaneous bonding + antibonding contributions near the frontier" case that's most diagnostic of an electronic instability.

Two variants reported, following the same raw/normalized pattern already used for the percolation weight and the min-cut descriptor:

- **Raw**: $W_\text{antibond}$ itself, in LOBSTER's ICOHP convention units (labelled "eV" by LOBSTER's own convention, though it is a Hamilton-population integral, not literally an energy).
- **Normalized**: $W_\text{antibond} \,/\, |\text{ICOHP}_\text{total occupied}|$, where the denominator is the "average" trace's own integrated ICOHP magnitude at E_ref — giving a dimensionless fraction ("what share of the total occupied bonding population near the frontier is antibonding"), more comparable across compounds of very different overall bond-strength scale.

Sign convention (empirically confirmed in step 1, reused here without re-deriving): **negative COHP/ICOHP = bonding, positive = antibonding.**

## 5. Validation on the 6 pilots

No independent ground truth exists for this new metric (none of the 6 pilots is a documented Peierls/JT-distorting compound), so validation here means: the pure numerical core passes hand-computable synthetic tests with known analytic answers, and the full pipeline behaves consistently (no crashes, sensible signs/magnitudes, expected invariants hold) on real data.

**Synthetic tests** (`tests/test_cohp_extraction.py::TestIntegrateAntibondingSynthetic`, 6 tests): pure bonding → 0; constant antibonding → width×height (up to a small, expected O(grid-spacing) edge-truncation from the strictly-one-sided window); a mixed bonding/antibonding case where only the antibonding part is counted (the test that would catch a "net signed integral" bug); window correctly ignores energies above E_ref; growing the window never decreases the integral (integrand is non-negative by construction); a too-small window with no grid point inside raises rather than silently returning 0.

**Real-pilot results** (ΔE = 1.0 eV, "average" trace):

| Compound | E_ref (eV, rel. to E_F) | $W_\text{antibond}$ raw | $W_\text{antibond}$ normalized |
|---|---:|---:|---:|
| NaCl (gapped) | −0.383 | 0.00991 | **0.111** |
| Si (gapped) | −0.259 | ~0.0 | ~10⁻⁵ |
| AlNi (metal) | 0.000 | 0.00052 | 0.00163 |
| LiBr (gapped) | −0.476 | ~0.0 | ~10⁻⁶ |
| C rhombohedral (gapped) | −0.361 | 0.0 | 0.0 |
| BeCu (metal) | 0.000 | 0.0 | 0.0 |

All 6 pilots — none documented as prone to a Peierls/JT-type distortion — show low-to-negligible antibonding population near the frontier, as expected for a first sanity pass. NaCl stands out as the highest of the 6 (normalized ≈ 0.11): this is **consistent with, not contradicted by**, step 1's independent qualitative finding that NaCl's "average" COHP trace was the one pilot compound flagged as having antibonding states below E_F at all (`has_antibnd_states_below_efermi`) — a reassuring cross-check between the two steps rather than a new, unrelated result. Physically plausible explanation: the "average" trace mixes in same-species (Na–Na, Cl–Cl) long-range contacts, which are closed-shell/Pauli-repulsion-dominated and often carry some antibonding character even in an otherwise cleanly-bonded ionic solid — worth revisiting with a per-bond-pair-type breakdown rather than only the all-pairs "average" if this direction is pursued further (see Limits).

AlNi shows a real window-sensitivity effect: normalized value grows from 0.00004 (ΔE=0.5) to 0.0016 (ΔE=1.0), i.e. its antibonding character sits mostly in the 0.5–1.0 eV range below E_F rather than immediately at the frontier — the kind of window-dependent detail this metric is specifically built to surface, that an integrated-over-everything ICOHP could never show.

## 6. What this does NOT establish yet

- No claim that this metric predicts anything (stability, persistence, or otherwise) — that requires the same kind of correlation testing already applied to the percolation weight and min-cut, against `energy_above_hull`/`formation_energy_per_atom` or (better, per mission #3's own conclusion) real persistence labels. Not done here.
- Only the "average" (all-bonds) trace was used; a per-bond-pair or per-orbital breakdown might be more diagnostic and is a natural next refinement, not attempted here.
- 6 compounds, none independently known to be electronically unstable — this validates the *pipeline*, not the *hypothesis*. A dataset including at least one well-documented Peierls/JT-distorting compound would be needed to test whether the metric actually discriminates, before any extension to the 186-compound dataset is worth the compute.

---

*Code*: `cohp_extraction.py` (`frontier_reference_energy`, `integrate_antibonding_in_window`, `antibonding_population_near_frontier`), `percolation_path.py` untouched. *Tests*: `tests/test_cohp_extraction.py` (14/14 passing; full project suite 34/34).
