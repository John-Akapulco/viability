# Does Δ(ICOHP) distinguish stable, metastable, and unstable compounds?

**This is the project's central research question, restated 2026-08-16**: not "which post-processed ICOHP/ICOBI descriptor correlates best with a stability target" (the framing every earlier mission used, `percolation_path.py` included — see `README.md`'s appendix section for that lineage), but specifically **does Δ(ICOHP) for the decomposition-into-elements reaction (`reaction_icohp.py`/`reaction_analysis`, the endobondic/exobondic axis) discriminate thermodynamically stable, metastable, and unstable compounds?** Two results, computed by `analysis/test_delta_icohp_viability.py`:

**1. Concordance against 7 reference validation reactions is essentially exact.** Lin's Concordance Correlation Coefficient (CCC — the correct test for *agreement*, not just correlation, between our computed values and the reference ones) is **0.999998** across the 7 worked reactions, sign agreement **7/7**, mean absolute residual **0.76 kJ/mol** (max 1.7 kJ/mol, fully attributable to the reference values' own rounding of their reported per-bond eV values). This is not a new finding — `tests/test_reitz_dronskowski_validation.py` already checked each case individually — but it had not previously been summarized as a single concordance statistic.

**2. Extended to every element and compound computed across the whole project to date (281 case-1 reactions, pooled without regard to which historical campaign or extension batch a compound came from — see §2), Δ(ICOHP) shows a real but not-yet-significant trend in the expected direction, and no significant sign-based prediction of experimental viability.** The endobondic fraction decreases monotonically exp_metastable (27.1%) > exp_stable (20.3%) > theo_metastable (11.9%) — physically the right direction (compounds known to exist, especially metastable ones, are more often "protected" by a bonding-derived kinetic barrier than theoretical-only ones) — but Kruskal-Wallis across these three groups is not significant (p=0.245, n=175), and a 2×2 test of endobondic/exobondic against experimental-vs-theoretical-only is also not significant (Fisher exact p=0.10, n=275). Against continuous stability targets, Δ(ICOHP)/atom correlates strongly with `formation_energy_per_atom` (ρ=−0.50, p=6.3×10⁻¹⁹, n=275 — this is the same relationship reported in `REPORT_reaction_icohp.md`, sign-flipped to this module's products−reactants convention) but only weakly and non-significantly with `energy_above_hull` (ρ=−0.11, p=0.068, n=279).

## 1. Reference concordance test

| | |
|---|---:|
| n | 7 |
| Lin's CCC | **0.999998** |
| Pearson r | 0.999999 (p=3.5×10⁻¹⁵) |
| Sign agreement | 7/7 |
| Mean \|residual\| | 0.76 kJ/mol |
| Max \|residual\| | 1.70 kJ/mol |
| RMSE | 0.98 kJ/mol |

| Reaction | Reference (kJ/mol) | Computed (kJ/mol) | Residual |
|---|---:|---:|---:|
| Pb(N₃)₂ → Pb + 3N₂ | +1345 | +1344.3 | −0.7 |
| S₄N₂ → ½S₈ + N₂ | +258 | +258.5 | +0.5 |
| S₄N₄ → ½S₈ + 2N₂ | +399 | +397.3 | −1.7 |
| ZnSn → Zn + Sn | −337 | −336.6 | +0.4 |
| CaO[sphalerite] → CaO[rocksalt] | −79 | −78.9 | +0.1 |
| CaN → ⅓Ca₃N₂ + ⅙N₂ | −205 | −204.8 | +0.2 |
| Mn₂O₇ → 2MnO₂ + 3/2 O₂ | −186 | −187.7 | −1.7 |

Why CCC and not just Pearson/Spearman: a correlation coefficient is satisfied by any monotonic relationship, including one with a large constant offset or a scale factor (e.g. if every computed value were systematically 50 kJ/mol too high, Pearson r would still read ≈1). CCC penalizes both — it specifically tests whether points fall on the *identity line* (computed = reference), which is the actual claim being tested here (do we reproduce the reference numbers, not just track them). A CCC this close to 1 means the two methods agree essentially within the reference values' own reporting precision, not merely that they move together.

## 2. Extended to every compound computed to date, pooled across campaigns

**This section deliberately does not organize by which historical campaign or extension batch (main 186-compound campaign, `extension`/`extension2`/`extension3`/`extension4`) a compound came from** — those boundaries reflect the order in which compute jobs happened to be submitted across sessions, not anything about the compounds' chemistry or stability, and grouping by them would be exactly the kind of between-group confound this project's own methodology exists to catch (see `REPORT_dimensionality_mincut.md` §5 for a worked example of what happens when that check isn't applied). Every compound with a computed case-1 (decomposition-into-elements) reaction is pooled into one population: **281 reactions**, spanning the original 186-compound campaign and all four extension batches together.

### 2.1 Continuous stability targets

| Target | n | ρ | p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | 275 | **−0.5016** | **6.3×10⁻¹⁹** |
| `energy_above_hull_eV_per_atom` | 279 | −0.1096 | 0.068 |

(Sign convention: `reaction_analysis`'s own products−reactants, i.e. the *opposite* of `reaction_icohp.py`'s reactant−products column reported elsewhere in this project — same magnitude, flipped sign; see `delta.py`'s docstring.) The formation-energy result is the same relationship already reported in `REPORT_reaction_icohp.md` (there ρ=+0.5016 in `reaction_icohp.py`'s own sign convention); repeated here to confirm it holds on the exact same pooled population this section uses, not a different subset. The hull-distance result remains weak and borderline, same conclusion as every earlier mission that tested `energy_above_hull` in this project.

### 2.2 Sign-based (endobondic/exobondic) discrimination

**Experimental-viability proxy (`theoretical` flag: `False` = the compound has been experimentally realized, `True` = theoretical-only, never synthesized as far as Materials Project's record indicates):**

| | endobondic | exobondic |
|---|---:|---:|
| experimental (n=177) | 47 | 130 |
| theoretical-only (n=98) | 17 | 81 |

Fisher's exact test: p=0.10, odds ratio=1.72 (endobondic is ~1.7× more common among experimentally-realized compounds than theoretical-only ones — the expected direction, since endobondic is meant as a kinetic-barrier signature that lets a compound persist despite favorable decomposition thermodynamics — but not significant at this sample size). Mann-Whitney on the continuous Δ(ICOHP)/atom itself, same split: experimental median −0.868 eV, theoretical-only median −1.442 eV (more exobondic), p=0.14 — same trend, same non-significance.

**Finer-grained (`family`, three-way, main-campaign compounds only, n=175):**

| family | n | median Δ(ICOHP)/atom (eV) | % endobondic |
|---|---:|---:|---:|
| exp_metastable | 59 | −1.077 | **27.1%** |
| exp_stable | 59 | −1.589 | 20.3% |
| theo_metastable | 57 | −1.659 | **11.9%** |

Kruskal-Wallis across the three groups: H=2.81, p=0.245 — **not significant**, but the ordering is exactly the physically-motivated one: `exp_metastable` compounds (known to exist despite being off the convex hull — the category the endobondic/exobondic framework is specifically about) have the highest endobondic fraction and the least-negative median Δ(ICOHP), `theo_metastable` compounds (never synthesized) have the lowest endobondic fraction and the most-negative median. This is a real, monotonic, mechanistically-sensible trend that does not (yet) clear the bar of statistical significance.

## 3. Reading

**The 7 reference validation cases are reproduced essentially exactly (§1) — the implementation is correct.** What §2 shows is that extending the same sign-based logic to a much larger, much more heterogeneous, and generally much more thermodynamically stable population (routine oxides/halides/nitrides/sulfides, not compounds specifically chosen because they sit near a decomposition edge) does not yet produce a statistically significant stable/metastable/unstable discriminator, though the trend points the right way. This is not a contradiction: the 7 reference cases were themselves chosen because they are informative edge cases (Pb(N₃)₂, S₄N₂/S₄N₄, ZnSn, and Mn₂O₇ are all borderline-stability compounds where a kinetic-barrier explanation is specifically needed) — most of this project's 281 compounds are not that; they are comfortably stable or comfortably not, where the endobondic/exobondic distinction may simply carry less discriminating information because the thermodynamic driving force (§2.1's formation-energy correlation) already dominates.

## 4. Limits and next steps

- **The family-based Kruskal-Wallis (§2.2) is the most promising unresolved lead in this report** — real monotonic trend, plausible mechanism, not yet significant at n=175. Worth revisiting once more `exp_metastable`-labeled compounds (the category most analogous to the 7 reference test cases) are added to the dataset, rather than assumed to be noise.
- **`classify.py`'s full `ViabilityLabel` machinery (not just `BondingLabel`) was not used here.** It requires a `delta_energy` for the *same* decomposition-into-elements reaction; for the vast majority of this project's compounds, decomposing into elements is strongly endothermic (i.e. `delta_energy` in `reaction_analysis`'s convention is strongly positive — formation from elements is exothermic, the normal case for stable inorganic compounds), which would trivially classify almost everything as `STABLE_ON_HULL` and tell us nothing beyond what `energy_above_hull` already says directly. `BondingLabel` alone (this report's approach) is the informative half of the classifier for this population; a proper `ViabilityLabel` test would need a population deliberately chosen to sit near a decomposition edge, closer to the 7 reference cases than to this project's general-purpose 186+extension dataset.
- **This section's `formation_energy_per_atom` result is not a new finding** — it is `REPORT_reaction_icohp.md`'s own headline, recomputed here in `reaction_analysis`'s sign convention to confirm consistency on the pooled (not campaign-split) population. Readers wanting the full stratified breakdown (`bond_type`, `is_metal`) should read that report; this one focuses specifically on the stable/metastable/unstable framing requested for this mission.
- **Case 2 (polymorph transition) and case 3 (decomposition into a compound + elements) are out of scope for this report** — both use different reaction topologies than the case-1 decomposition-into-elements reactions tested here; see `REPORT_reaction_icohp.md` §6 for case 2's own (also inconclusive) result.

---

*Code*: `analysis/test_delta_icohp_viability.py`. *Data*: `analysis/stats_summary_delta_icohp_viability.json`, reusing `analysis/reaction_icohp_case1.csv`, `analysis/reaction_analysis_case1_full.csv`, `mp_dataset/formation_energies.json` (no new DFT calculations for this report). *Validation fixtures*: `tests/test_reitz_dronskowski_validation.py`, `tests/fixtures/reitz_dronskowski_cases.yaml`.
