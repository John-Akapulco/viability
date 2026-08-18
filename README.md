# viability

**Central research question: does Δ(ICOHP) —
the reaction-ICOHP descriptor for a compound's decomposition into
elements, and specifically the sign-based endobondic/exobondic
classification it supports —
distinguish thermodynamically stable, metastable, and unstable
compounds?** Tested over every element and compound with a computed
case-1 decomposition reaction (`mp_dataset/structures/`, 597 compounds,
518 case-1 reactions), pooled without regard to which batch a compound
came from — see
**[`analysis/REPORT_delta_icohp_viability.md`](analysis/REPORT_delta_icohp_viability.md)**.

**Headline**: concordance against 7 reference validation reactions with
independently known ΔICOHP values is essentially exact (Lin's Concordance Correlation
Coefficient = 0.999998, sign agreement 7/7, mean residual 0.76 kJ/mol —
this is a genuine agreement test against the identity line, not just a
correlation, and the implementation is validated by it). Extended to the
project's full 518-reaction history, the sign-based endobondic/exobondic
split **now does significantly discriminate** experimentally-realized
compounds from theoretical-only ones (Fisher exact p=0.0053, odds ratio
1.82; Mann-Whitney on the continuous Δ(ICOHP)/atom itself, same split,
p=0.0088) — a result the marginal-formation-energy and max-hull-distance
batches were built specifically to test. The finer three-way split
(exp_metastable/exp_stable/theo_metastable) is unchanged from before
(Kruskal-Wallis p=0.278, n=177 — those newer batches carry a
`theoretical` flag but not yet a `family` label, so this specific test
has not grown). Against continuous stability targets, Δ(ICOHP)/atom
still correlates with `formation_energy_per_atom` (ρ=−0.297,
p=7.3×10⁻¹²) but more weakly than on the earlier, smaller population
(ρ=−0.50) and no longer even borderline with `energy_above_hull`
(ρ=−0.058, p=0.19) — the added batches' deliberately marginal/extreme
thermodynamics dilute a rank correlation partly carried by the earlier
population's narrower stability range. See the report for the full
picture, and its §4 for the identified next step (extend `family`
labels to the newer batches so the three-way split can be tested at the
same scale as the binary one).

This question is answered by two modules working together:
`reaction_icohp.py` (the original ICOHP-analog-of-formation-energy
metric, mission #5) and the schema-driven `reaction_analysis/` package
(meant to eventually supersede it) with its endobondic/exobondic
classifier. A second, independent line of evidence —
`cohp_extraction.py`'s antibonding-population-near-the-frontier metric,
mission #4 — asks a related but distinct question (how COHP is
distributed in energy, not the total reaction-level ICOHP balance) and
is documented alongside it below. The project's original graph-based
methodology (`Percolation_viability/percolation_path.py`, mission #1;
`Percolation_viability/network_dimensionality.py`/
`Percolation_viability/periodic_mincut.py`, mission #3) established the
periodic-graph representation of ICOHP/ICOBI data and the statistical
conventions every later mission still reuses (Spearman correlation only,
mandatory `bond_type`/`is_metal` stratification, n<15 subgroups always
flagged, no SISSO until a signal is clean and well-powered per stratum)
— but its own headline correlations did not hold up at the current
dataset scale (see the Appendix) and it is no longer this project's
central question, so it now lives in its own `Percolation_viability/`
directory as background/appendix material rather than a current
descriptor.

## Descriptors (ordered by relevance to the central question, not by mission number)

1. **[Reaction-ICOHP (Δ(ICOHP)) and endobondic/exobondic classification](#reaction-icohp-δicohp-and-endobondicexobondic-classification-mission-5--central-question)** — mission #5, `reaction_icohp.py` + `reaction_analysis/`. The project's central question as of 2026-08-16. Numerically the strongest correlation in the project (ρ=0.502 vs. `formation_energy_per_atom`) and the first descriptor whose signal survives an `is_metal` stratification; near-exact concordance with 7 reference validation reactions (CCC=0.999998); sign-based viability discrimination shows a real but not-yet-significant trend.
2. **[Antibonding population near the frontier (E_F/VBM)](#antibonding-population-near-the-frontier-ef-vbm-mission-4)** — mission #4, `cohp_extraction.py`. A related but distinct energy-resolved-COHP question. Its `bond_type=covalent` subgroup is the most robust bond-type-stratified result in the project.
3. **[Appendix: original graph-based methodology](#appendix-original-graph-based-methodology-percolation-dimensionality-min-cut--missions-1-and-3)** — `Percolation_viability/percolation_path.py` (mission #1) and `Percolation_viability/network_dimensionality.py`/`Percolation_viability/periodic_mincut.py` (mission #3). Historical foundation (periodic graph construction, statistical conventions); no longer the project's central question, and min-cut's own headline correlation did not survive the dataset's growth to 349 compounds (diagnosed, not just noted — see the appendix).

---

## Reaction-ICOHP (Δ(ICOHP)) and endobondic/exobondic classification (mission #5 — central question)

A thermochemistry-flavored question: not a compound's own bonding
topology or energy distribution in isolation, but whether its total ICOHP
is "worth more" than the same atoms would have in a reference
configuration — the ICOHP analog of `formation_energy_per_atom`. Module
`reaction_icohp.py`, three reaction types (decomposition into elements,
polymorph comparison, decomposition into a compound + elements), balanced
via `pymatgen.analysis.reaction_calculator.Reaction`. Defined and
validated on n=1–5 hand-worked real examples (Ca3N2, Mn2O7, carbon
allotropes, TiO2 high-pressure polymorphs) in
**[`analysis/METRIC_DEFINITION_reaction_icohp.md`](analysis/METRIC_DEFINITION_reaction_icohp.md)**,
which already flagged the key caveat before any statistics were run:
ICOHP sees orbital-overlap bond population, not electrostatic/Madelung or
van der Waals energy, so strongly ionic compounds and van-der-Waals-bound
polymorphs (e.g. graphite) are expected to misbehave.

Case 1 (decomposition into elements), using 62 elemental reference
calculations (`mp_dataset/download_elements_reference.py` + hand-picked
extension compounds), computed for **518 case-1 reactions** pooled
across every batch — 0 missing elemental references. Tested in
`analysis/stats_analysis_reaction_icohp.py`. See
**[`analysis/REPORT_reaction_icohp.md`](analysis/REPORT_reaction_icohp.md)**.
Headline: ρ=0.297, p≈0 (n=511) against `formation_energy_per_atom` — a
real, highly significant correlation, though weaker than it looked at
smaller scale. The between-group picture has shifted shape rather than
resolved: `bond_type=ionic` (ρ=0.360, p=0.0004) and the newer
`bond_type=mixed`/Zintl category (ρ=−0.525, opposite sign, p≈0) both
carry their own significant signal, `bond_type=covalent` too but flips
sign (ρ=−0.354, p=0.023), and only `bond_type=metallic` (the majority
class, 316/511 rows) stays non-significant. **The coarser `is_metal`
split — the one that used to survive stratification when `bond_type`
didn't — has since converged onto `bond_type=metallic`'s exact
population** (identical n, ρ, p) and now inherits its non-significant
result; `bond_type`, not `is_metal`, is the split that carries
information at the current dataset composition. Against
`energy_above_hull` the correlation is weak and no longer even
borderline at the pooled level (ρ=0.058, p=0.19), though
`bond_type=mixed` again survives on its own (ρ=−0.307, p=0.016). On a
288-compound subset shared with every other descriptor in the project,
antibonding population (mission #4) edges out reaction ICOHP on both
targets — reaction ICOHP is a strong descriptor here, not unambiguously
*the* strongest anymore. Case 2 (polymorph comparison), at **74
polymorph groups** (up from an original 8): 47.3% agreement between
most-bonding and most-stable member, statistically indistinguishable
from the 42.4% chance baseline (P=0.22) — bonding does not track
polymorph stability, the same conclusion as at every earlier scale.
Case 3 (decomposition into a compound + elements) remains not
attempted; see the report for why.

**`reaction_analysis/` (schema-driven redesign of this axis)**: a
from-scratch Pydantic schema (`CompoundEntry`, `Reaction`,
`ReactionResult`) meant to eventually cover all three reaction types
above through one common, testable data model, rather than the ad hoc
`reaction_icohp.py` functions above. Ships with `parse_lobster.py`
(builds a `CompoundEntry` from `ICOHPLIST.lobster`/`ICOBILIST.lobster` +
structure, with an explicit regression test confirming LOBSTER lists each
periodic bond once, not once per direction — the assumption
`sum_total_eV`'s unfiltered summation depends on), `balance.py`
(element-by-element stoichiometric balance checking, coefficient
auto-derivation for decomposition-into-elements), and `delta.py` (the
three ΔICOHP/ΔICOBI normalizations — per formula unit, per atom, and a
non-conservative per-bond diagnostic — computed together, never one in
isolation). Populated with a first real-data batch (6 compounds, 3
`decomposition_to_elements` reactions), cross-checked 3/3 against
`reaction_icohp.py`'s numbers, then extended to cover the full case-1
population (`analysis/populate_reaction_analysis_case1_full.py`,
`reactions_dataset/` at the repo root) — matches `reaction_icohp_case1.csv`
row-for-row (518/518) after the known sign flip. Its sign convention
(products − reactants) is the *opposite* of `reaction_icohp.py`'s — the
two are not interchangeable.

**Endobondic/exobondic classification**: `nearest_neighbor.py`
(first-coordination-shell bond filtering — a relative, self-calibrating
gap detector, no hardcoded Angstrom cutoff), `classify.py` (`BondingLabel`
endobondic/exobondic from ΔICOHP's sign, and a deliberately more cautious
`ViabilityLabel` — `UNSTABLE_NONEXISTENT` always carries a warning that an
exobondic sign never proves non-existence, e.g. Mn2O7, a real compound
that decomposes by slow, gradual O2 loss rather than abrupt bond
rupture), and `units.py` (eV ↔ kJ/mol). `parse_lobster.py` gained an
opt-in `bond_filter="nearest_neighbor"` (default stays `"unfiltered"`, to
avoid silently changing the already-validated case-1 population above),
and `delta.py`'s `ReactionResult` now carries a `bonding_label`.

**Validation against 7 reference reactions and extension to the full
dataset** (**[`analysis/REPORT_delta_icohp_viability.md`](analysis/REPORT_delta_icohp_viability.md)**,
`analysis/test_delta_icohp_viability.py`): all 7 reference validation
cases reproduce to a Lin's CCC of 0.999998 (sign agreement 7/7, mean
residual 0.76 kJ/mol) — independent of this project's own CSP data.
Extended to every case-1 reaction with a computed result (518, pooled
across all batches — not split by which one a compound happened to
arrive in), the endobondic/exobondic sign **now significantly**
distinguishes experimentally-realized from theoretical-only compounds
(Fisher exact p=0.0053, odds ratio 1.82; Mann-Whitney on the continuous
Δ(ICOHP)/atom, p=0.0088) — the marginal-formation-energy and
max-hull-distance batches supplied the borderline-stability contrast
this test needed. The finer three-way split
(exp_metastable/exp_stable/theo_metastable, still n=177) remains
not-yet-significant (Kruskal-Wallis p=0.278) because those batches
carry a `theoretical` flag but not yet a `family` label. **Not yet
connected to real LOBSTER data for a full `ViabilityLabel`
classification** (needs a population deliberately chosen to sit near a
decomposition edge, closer to the 7 reference cases than to this
project's general-purpose dataset) — the concrete next step for this
central question, alongside extending `family` labels to the newer
batches.

## Antibonding population near the frontier (E_F/VBM) (mission #4)

A distinct question from integrated ICOHP/ICOBI (a single number per
bond): not *how much* bonding a compound has in total, but *how COHP is
distributed in energy* — specifically, whether the highest-energy
occupied states carry antibonding character, by analogy with
Peierls/Jahn-Teller electronic instabilities. `cohp_extraction.py`, built
on `pymatgen.io.lobster.outputs.Cohpcar` / `pymatgen.electronic_structure.cohp.CompleteCohp`
(no hand-rolled COHPCAR parsing). Cross-validated against the
already-validated `ICOHPLIST.lobster` across 558 bond labels in the 6
pilot compounds (exact match for 5/6, 1e-5 eV for the 6th); metal/gap
classification cross-checked against Materials Project rather than
derived locally, which caught a real pitfall (our LOBSTER-oriented coarse
k-mesh spuriously suggests small gaps for two known metals). See
**[`analysis/REPORT_cohp_feasibility.md`](analysis/REPORT_cohp_feasibility.md)**
(extraction pipeline validation) and
**[`analysis/METRIC_DEFINITION_antibonding.md`](analysis/METRIC_DEFINITION_antibonding.md)**
(the window/metric definition itself: one-sided window below E_F/VBM,
integrated antibonding-only COHP, raw + normalized). Validated on the 6
pilots only (synthetic numerical tests + real-data sanity checks,
`tests/test_cohp_extraction.py`) before any extension.

Extended to the full 186-compound dataset in `analysis/compute_antibonding_all.py`
and tested against `formation_energy_per_atom` in
`analysis/stats_analysis_antibonding.py`; both rerun unmodified over the
349-compound dataset once the `extension4` batch (89 alkali/alkaline-earth
binaries) was added (347/349 succeeded, 2 pre-existing COD-sourced gaps).
See **[`analysis/REPORT_antibonding.md`](analysis/REPORT_antibonding.md)**.
Headline: the normalized metric reaches ρ=−0.282, p=9.2×10⁻⁸ (n=346) —
weaker than the ρ=−0.328 (n=186) originally reported, but still
comfortably ahead of every descriptor except reaction ICOHP. The sign
(more antibonding population near the frontier associates with *more
negative*, i.e. more stable, formation energy) is the opposite of the
naive Peierls/Jahn-Teller reading, and diagnostics show the global number
is substantially a between-group effect — but `is_metal=False` now also
survives on its own (n=157, p=0.001), which it did not at n=186.
**`bond_type=covalent` (n=33) not only still holds up under
stratification but got *stronger* with more data** (ρ=−0.551, p=0.0009,
vs. the original ρ=−0.490 at n=23) — still the single most robust
bond-type-stratified result in the project. The old `bond_type=ionic`
result (n=9, ρ=−0.683) did not replicate once the sample grew to n=53
(ρ=−0.184, not significant) — a clean demonstration of why n<15 rows are
flagged rather than trusted. See the report for the full within-group
diagnostics, the ΔE-sensitivity check (still robust across 0.5–2.0 eV),
and why the sign is not yet interpretable as confirming or refuting the
instability hypothesis.

`ICOBI`-based near-frontier windowing (as opposed to `ICOHP`/COHP) is a
natural extension not yet implemented — `percolation_path.py` already
treats ICOHP and ICOBI symmetrically as alternative edge weights, but
`cohp_extraction.py` is COHP/ICOHP-only for now.

---

## Appendix: original graph-based methodology (percolation, dimensionality, min-cut — missions #1 and #3)

The project's original angle, and the reason every descriptor above
shares a common data model — but no longer the project's central
question (see the top of this file), and demoted here after its own
headline correlation (min-cut, mission #3) turned out not to survive the
dataset's growth to 349 compounds (see below). Kept for the periodic
graph representation and statistical conventions every later mission
still reuses, and because `percolation_path.py`'s CLI remains a real,
usable tool independent of which mission is currently the project's
focus.

### Network dimensionality + periodic min-cut (mission #3)

Two descriptors, `percolation_path.py` untouched:
`network_dimensionality.py` (0D-3D classification via a relative bond-
strength threshold + BFS connectivity, since the existing graph's ~6 Å
cutoff makes almost everything look 3D) and `periodic_mincut.py` (minimum
total bond strength separating the crystal into two halves along a
direction, via `networkx` max-flow/min-cut on a finite ribbon graph — a
different physical question from the percolation weight: separability, not
traversability). Both validated on synthetic cases with known answers
(`tests/test_network_dimensionality.py`, `tests/test_periodic_mincut.py`)
before running on real data. See
**[`analysis/REPORT_dimensionality_mincut.md`](analysis/REPORT_dimensionality_mincut.md)**:
tested against `formation_energy_per_atom` (`mp_dataset/
fetch_formation_energy.py`) rather than `energy_above_hull`. Headline (on
the original 186-compound population, still exactly reproducible today):
min-cut (normalized) was the first descriptor in the project to reach a
significant *global* correlation (ρ=0.285, p=0.0001, n=186) — stronger
than the original percolation weight ever achieved — though likely driven
partly by between-bond-type clustering rather than a clean within-type
relationship; dimensionality alone does not separate formation energy
(Kruskal-Wallis p=0.19). **This pooled number is no longer trustworthy at
the current 349-compound scale (ρ=0.089, p=0.098, not significant) — not
because the original result was wrong, but because later, unrelated
extension batches diluted it: ~65 elemental-reference compounds (added
for `reaction_icohp.py`, sitting at `formation_energy_per_atom`≈0 by
construction) and `extension4`'s deliberately far-above-hull half (which
shows a *significant correlation of the opposite sign*, ρ=−0.40, p=0.011,
n=39) got pooled in without re-checking this descriptor specifically. See
`REPORT_dimensionality_mincut.md` §5 for the full breakdown. **That
opposite-sign result is itself now explained (§5.5): a between-anion-
chemistry confound (F- vs. N-azide-vs-O/P/S/Cl-containing systems), not
a cell-size/symmetry artifact or a real far-from-hull effect — it
collapses to ρ≈−0.07, not significant, once anion identity or exact
chemsys composition is controlled for.**

### `percolation_path.py` — the original descriptor (mission #1)

The project's first descriptor and the reason everything above shares a
common data model. Post-processing on a relaxed crystal structure and its
ICOHP/ICOBI (LOBSTER) data: no molecular dynamics, no NEB calculation, no
physical supercell — periodicity is handled by vector-labeling the bond
graph's edges (a voltage graph / gain graph), and the minimum-weight
non-contractile path is found by Dijkstra search over the extended state
space `(atom, cumulative translation)`.

Tested against `energy_above_hull` (mission #1) and, on the primitive-vs-
conventional-cell question, revisited on the 6 pilot compounds (mission
#2, see below): **no significant correlation was found on its own** — see
[Analysis](#analysis-of-percolation_pathpy-missions-1-2) — but the graph
construction and Dijkstra machinery here are exactly what every later
descriptor's own graph (min-cut's ribbon graph, the antibonding metric's
energy window) still builds on or was designed in explicit contrast to.

#### Installation

```bash
pip install -r requirements.txt   # pymatgen
```

Requires Python ≥ 3.9 (uses `pymatgen.io.lobster.outputs.Icohplist`).

#### Expected input layout

```
dataset/
  compound_A/
    POSCAR (or CONTCAR, or *.cif)      # relaxed structure, primitive cell
    ICOHPLIST.lobster                   # LOBSTER output, with translation columns
    ICOBILIST.lobster                   # optional
  compound_B/
    ...
```

Each subdirectory of `--root` is treated as an independent compound. The
`ICOHPLIST.lobster`/`ICOBILIST.lobster` file must include the per-bond
lattice translation column (`tx ty tz`) — the standard case since LOBSTER
≥ 3. This is the information that lets the labeled periodic graph be
reconstructed without ever physically duplicating atoms.

#### Usage

```bash
python Percolation_viability/percolation_path.py --root dataset --metric icohp --output results.csv
python Percolation_viability/percolation_path.py --root dataset --metric icohp,icobi --output results.csv --also-json results.json
python Percolation_viability/percolation_path.py --root dataset --metric icohp --bond-pair Fe-O --output results.csv
```

Main options:
- `--metric icohp|icobi|icohp,icobi`: which quantity/quantities to use as
  edge weight (|ICOHP| or |ICOBI|). Both can be computed side by side for
  comparison.
- `--bond-pair Fe-O`: restricts the graph and the aggregates to bonds
  between these two species (otherwise every bond in the file is used).
- `--coord-bound N`: bound on the cumulative-translation search space
  explored during pathfinding. Defaults to a value derived automatically
  from the largest translation vector present in the input data (so it's
  directly tied to the upstream ICOHP/ICOBI calculation's cutoff radius,
  not a convergence parameter to tune by hand).

#### Output

One record per compound (CSV, one row per compound; or JSON, nested
structure with per-direction detail). Main columns (prefixed by the
metric name, e.g. `icohp_...`):

- `*_percolation_weight_min`: weight of the weakest percolation path,
  across all directions — the main descriptor.
- `*_percolation_direction`: the corresponding direction (`a`, `b`, or `c`).
- `*_percolation_weight_a/b/c` and `*_percolation_status_a/b/c`: minimum
  weight per direction, with an explicit status (`ok` or `disconnected`
  if no non-contractile path exists in that direction — never a silent
  infinite value).
- `*_sum`, `*_mean`, `*_min`, `*_max`: classic aggregates over the same
  bonds, for direct comparison against the percolation descriptor.
- `error`, `warnings`: per-compound diagnostics (batch processing never
  stops on a failed compound — the error is logged in the corresponding
  row and the batch continues).

#### Algorithm (summary)

1. Each ICOHP/ICOBI bond becomes an edge `(atom_i, atom_j, (nx,ny,nz))`
   of the periodic graph, weighted `|value|`. Both directions of the edge
   are added (opposite translation the other way).
2. For each lattice direction `a`, `b`, `c`, the minimum-weight path in
   the state space `(atom, cumulative translation)` connecting
   `(atom, (0,0,0))` to `(atom, direction)` is searched for each possible
   starting atom. This is plain Dijkstra (weights ≥ 0 since they're
   absolute values), valid because the state space is finite once bounded
   by `--coord-bound`.
3. The minimum over starting atoms gives that direction's percolation
   weight; the minimum over the three directions gives the main
   descriptor.
4. If no target state is reached in a direction (bond network
   disconnected along that direction given the ICOHP/ICOBI cutoff
   radius), this is reported explicitly (`status: disconnected`), never
   silently as an infinite or zero weight.

#### Tests

```bash
python -m unittest discover -s tests -v
```

Covers: isotropic case (identical weight in all 3 directions), anisotropic
case (weakest direction correctly identified), disconnected case (reported
explicitly, no crash or wrong value), a case where a cheaper multi-hop
indirect path beats a more expensive direct bond (validates that the
algorithm actually solves a global shortest-path problem, not just picking
the weakest bond), and an end-to-end integration test with a synthetic
`ICOHPLIST.lobster` parsed through pymatgen.

#### Examples

`examples/dataset/` contains three toy compounds illustrating the three
cases above (anisotropic, indirect path cheaper than the direct bond,
disconnected direction):

```bash
python Percolation_viability/percolation_path.py --root examples/dataset --metric icohp --output /tmp/example_results.csv -v
```

#### Analysis of `percolation_path.py` (missions #1-#2)

`analysis/` contains a statistical study of the percolation descriptor
against thermodynamic stability (`energy_above_hull`) over the 186-compound
dataset in `mp_dataset/structures/` (60 experimental-stable + 60
experimental-metastable + 60 theoretical-metastable, from `mp_dataset/
select_campaign.py`, plus the 6-compound pilot). See
**[`analysis/REPORT.md`](analysis/REPORT.md)** for the full write-up:
Spearman correlations (overall and per bond_type) for the raw and
normalized percolation weight vs. the classic ICOHP aggregates, a reference
logistic regression (stable vs. metastable) with cross-validated AUC,
figures, compute-time table, and a documented limits/next-steps section.
Headline result: no significant overall correlation with `energy_above_hull`
at this sample size, and no clear predictive edge over the classic
aggregates in the one subgroup (metallic) that shows a nominally
significant signal — see the report for why, and what would need to change
before this justifies a symbolic-regression (SISSO) pass.

Pipeline: `analysis/build_dataset.py` (join `percolation_path.py` output +
MP metadata → `analysis/percolation_vs_hull.csv`) then
`analysis/stats_analysis.py` (correlations, logistic regression, figures →
`analysis/stats_summary.json` + `analysis/figures/`).

`analysis/REPORT_conventional_pilot.md` (mission #2): tested the pilot
report's hypothesis that the primitive/conventional cell choice explains
the near-zero correlation above, on the 6 pilot compounds. Verdict: the
choice does change the percolation weight substantially for compounds
whose cell actually differs (3x-55x, well past DFT noise), but the
strongest bond still never participates in the winning path, and the
weight gets *smaller* rather than larger as hypothesized — a bigger cell
just gives the minimum-weight search more long-range weak bonds to
exploit. Extension to the full dataset was not pursued based on this
result.

#### Known limitations (`percolation_path.py` CLI)

- Assumes `ICOHPLIST.lobster`/`ICOBILIST.lobster` come from the same
  LOBSTER calculation as the supplied structure (same atom order, same
  cell); no cross-check of structure/POSCAR-LOBSTER consistency beyond
  atom-index validation.
- Per-direction Dijkstra complexity is `O(atoms)` runs over a state space
  of size `atoms × (2·coord_bound+1)³`; suited to cells of a few dozen to
  a few hundred atoms per compound, not to giant cells.
