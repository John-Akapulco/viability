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
p=0.0088), and the full `ViabilityLabel` (sign combined with
thermodynamics, `classify_viability()`) discriminates far more strongly
still (Fisher exact **p=2.8×10⁻¹⁵**, odds ratio **6.0** — experimentally-
realized compounds are ~6× less likely to be predicted
`UNSTABLE_NONEXISTENT` than theoretical-only ones) — a result the
marginal-formation-energy and max-hull-distance batches were built
specifically to test. The finer three-way split
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
picture, and its §4 for the identified next steps (extend `family`
labels to the newer batches so the three-way split can be tested at the
same scale as the binary ones; check whether the `ViabilityLabel` result
is itself confounded).

This question is answered by two modules working together:
`reaction_icohp.py` (the original ICOHP-analog-of-formation-energy
metric, mission #5) and the schema-driven `reaction_analysis/` package
(meant to eventually supersede it) with its endobondic/exobondic
classifier. A second, independent line of evidence —
`cohp_extraction.py`'s antibonding-population-near-the-frontier metric,
mission #4 — asks a related but distinct question (how COHP is
distributed in energy, not the total reaction-level ICOHP balance) and
is documented alongside it below.

## Descriptors (ordered by relevance to the central question, not by mission number)

1. **[Reaction-ICOHP (Δ(ICOHP)) and endobondic/exobondic classification](#reaction-icohp-δicohp-and-endobondicexobondic-classification-mission-5--central-question)** — mission #5, `reaction_icohp.py` + `reaction_analysis/`. The project's central question as of 2026-08-16. Numerically the strongest correlation in the project (ρ=0.502 vs. `formation_energy_per_atom`) and the first descriptor whose signal survives an `is_metal` stratification; near-exact concordance with 7 reference validation reactions (CCC=0.999998); sign-based viability discrimination shows a real but not-yet-significant trend.
2. **[Antibonding population near the frontier (E_F/VBM)](#antibonding-population-near-the-frontier-ef-vbm-mission-4)** — mission #4, `cohp_extraction.py`. A related but distinct energy-resolved-COHP question. Its `bond_type=covalent` subgroup was long the most robust bond-type-stratified result in the project — retracted at the current 588-compound scale (see below); `bond_type=mixed` and `is_metal=False` now carry the surviving signal.

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
this test needed. **The full `ViabilityLabel`** (sign combined with
`delta_energy`, `classify_viability()`, run on real LOBSTER-derived
Δ(ICOHP) via `analysis/compute_case1_viability.py`) **discriminates far
more strongly still**: Fisher exact **p=2.8×10⁻¹⁵**, odds ratio **6.0**
(experimentally-realized compounds ~6× less likely to be predicted
`UNSTABLE_NONEXISTENT`). The finer three-way split
(exp_metastable/exp_stable/theo_metastable, still n=177) remains
not-yet-significant (Kruskal-Wallis p=0.278) because those batches
carry a `theoretical` flag but not yet a `family` label — the concrete
next step for this central question, alongside checking whether the
`ViabilityLabel` result is itself confounded by `bond_type`/`is_metal`
or by max-hull's deliberately-extreme selection (not yet checked, see
the report's §4).

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

Extended to every compound with a `COHPCAR.lobster` file in
`analysis/compute_antibonding_all_full.py` (a generic, full-scale
successor to the original 186/349-compound-scoped
`compute_antibonding_all.py`), currently 588 compounds (one, Cu4Au, excluded for a catastrophic
LOBSTER basis-projection failure — 100% k-point orthonormalization
failure, `mp_metadata.json`'s `quality_excluded_reason` has detail).
See
**[`analysis/REPORT_antibonding.md`](analysis/REPORT_antibonding.md)**.
Headline: the normalized metric reaches ρ=−0.155, p=1.8×10⁻⁴ (n=579) —
markedly weaker than the ρ=−0.282 (n=346) reported at the previous
scale, but still real and significant. The sign (more antibonding
population near the frontier associates with *more negative*, i.e. more
stable, formation energy) is the opposite of the naive Peierls/Jahn-Teller
reading, and diagnostics show the global number is substantially a
between-group effect. **The project's previously most robust
bond-type-stratified result is retracted**: `bond_type=covalent`, which
survived clearly at n=33 (ρ=−0.551, p=0.0009), does *not* survive at the
current n=54 (ρ=−0.046, p=0.74) — the same small-subgroup fragility this
project's own methodology has now caught repeatedly, including this
same descriptor's own earlier `bond_type=ionic` result (n=9→53, also
didn't replicate). **New at this scale**: `bond_type=mixed` (Zintl) is
now the best-surviving stratum (ρ=−0.538, p=6.6×10⁻⁶, n=62), and
`is_metal=False` remains the most consistently-surviving stratum across
this project's successive scales (ρ=−0.333, p=6.7×10⁻⁷, n=212,
strengthened rather than weakened by the dataset's growth). See the
report for the full within-group diagnostics, the ΔE-sensitivity check
(still robust across 0.5–2.0 eV), the `energy_above_hull` comparison
(a different, sign-inconsistent picture within `bond_type=metallic`),
and why the sign is not yet interpretable as confirming or refuting the
instability hypothesis.

`ICOBI`-based near-frontier windowing (as opposed to `ICOHP`/COHP) is a
natural extension not yet implemented; `cohp_extraction.py` is
COHP/ICOHP-only for now.
