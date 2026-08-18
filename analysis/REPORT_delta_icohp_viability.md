# Does Δ(ICOHP) distinguish stable, metastable, and unstable compounds?

This is the project's central research question: does Δ(ICOHP) for the decomposition-into-elements reaction (`reaction_icohp.py`/`reaction_analysis`, the endobondic/exobondic axis) discriminate thermodynamically stable, metastable, and unstable compounds? Two results, computed by `analysis/test_delta_icohp_viability.py`:

**1. Concordance against 7 reference validation reactions is essentially exact.** Lin's Concordance Correlation Coefficient (CCC — the correct test for *agreement*, not just correlation, between our computed values and the reference ones) is **0.999998** across the 7 worked reactions, sign agreement **7/7**, mean absolute residual **0.76 kJ/mol** (max 1.7 kJ/mol, fully attributable to the reference values' own rounding of their reported per-bond eV values).

**2. Extended to every element and compound with a computed case-1 reaction (518 reactions, pooled without regard to which batch a compound came from — see §2), the coarse experimental-vs-theoretical-only split is now a significant discriminator.** Fisher's exact test on endobondic/exobondic vs. experimentally-realized/theoretical-only is **p=0.0053** (odds ratio 1.82), and Mann-Whitney on the continuous Δ(ICOHP)/atom, same split, is **p=0.0088** — endobondic character is markedly more common, and Δ(ICOHP)/atom markedly less negative, among compounds Materials Project records as experimentally realized. The finer three-way split (`exp_metastable`/`exp_stable`/`theo_metastable`) remains a real, monotonic, but not-yet-significant trend (Kruskal-Wallis p=0.278, n=177 — unchanged from before, since `family` labels are only assigned within the original campaign, not the newer batches; see §4). Against continuous stability targets, Δ(ICOHP)/atom correlates with `formation_energy_per_atom` (ρ=−0.297, p=7.3×10⁻¹²) — weaker than on the earlier, smaller population (ρ=−0.50) but still highly significant — and no longer even borderline with `energy_above_hull` (ρ=−0.058, p=0.19, up from p=0.068).

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

## 2. Extended to every compound with a computed case-1 reaction

Every compound with a computed case-1 (decomposition-into-elements) reaction is pooled into one population, regardless of which batch produced it — batch boundaries reflect submission order, not compound chemistry, and splitting by them would reintroduce exactly the kind of between-group confound this project's methodology is built to catch. Current population: **518 reactions**.

### 2.1 Continuous stability targets

| Target | n | ρ | p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | 511 | **−0.2969** | **7.3×10⁻¹²** |
| `energy_above_hull_eV_per_atom` | 515 | −0.0579 | 0.190 |

(Sign convention: `reaction_analysis`'s own products−reactants, i.e. the *opposite* of `reaction_icohp.py`'s reactant−products column reported elsewhere in this project — same magnitude, flipped sign; see `delta.py`'s docstring.) Both correlations are weaker than on the earlier, smaller population: the formation-energy result was ρ=−0.50 before the marginal-formation-energy and max-hull-distance batches were folded in, and the hull-distance result has moved from borderline (p=0.068) to clearly non-significant (p=0.19). This is consistent with those batches being deliberately selected for marginal or extreme thermodynamics rather than typical compounds — they add real chemical diversity but dilute a rank correlation that was partly carried by the original population's narrower stability range.

### 2.2 Sign-based (endobondic/exobondic) discrimination

**Experimental-viability proxy** (`theoretical` flag: `False` = experimentally realized, `True` = theoretical-only per Materials Project):

| | endobondic | exobondic |
|---|---:|---:|
| experimental (n=267) | 79 | 188 |
| theoretical-only (n=245) | 46 | 199 |

Fisher's exact test: **p=0.0053**, odds ratio=1.82 — endobondic character is significantly more common among experimentally-realized compounds than theoretical-only ones, the expected direction (endobondic is meant as a kinetic-barrier signature that lets a compound persist despite favorable decomposition thermodynamics). Mann-Whitney on the continuous Δ(ICOHP)/atom itself, same split: experimental median −0.868 eV, theoretical-only median −1.457 eV (more exobondic), **p=0.0088** — same direction, also significant.

**Finer-grained (`family`, three-way, n=177 — unchanged from the original campaign; the extension/marginal-formation-energy/max-hull batches carry a `theoretical` flag but not a `family` label, so this specific test has not grown):**

| family | n | median Δ(ICOHP)/atom (eV) | % endobondic |
|---|---:|---:|---:|
| exp_metastable | 59 | −1.038 | **28.8%** |
| exp_stable | 59 | −1.391 | 22.0% |
| theo_metastable | 59 | −1.601 | **13.6%** |

Kruskal-Wallis across the three groups: H=2.56, p=0.278 — still not significant, but the ordering is the same physically-motivated one as before: `exp_metastable` compounds have the highest endobondic fraction and the least-negative median Δ(ICOHP), `theo_metastable` compounds have the lowest endobondic fraction and the most-negative median.

### 2.3 Full `ViabilityLabel` (not just the sign) — now large enough to test

`classify_viability()` combines `BondingLabel`'s sign with `delta_energy` (`-formation_energy_per_atom`): a compound is only `UNSTABLE_NONEXISTENT` if decomposition into elements is *both* exobondic *and* thermodynamically favorable, a strictly stronger, more specific claim than the sign-only split in §2.2. Run via `analysis/compute_case1_viability.py` on the same pooled, no-batch-exclusion population (518 reactions, 511 classifiable — 7 excluded for missing `formation_energy_per_atom`, all COD-sourced or otherwise mp_id-less compounds):

| | experimental (n=267) | theoretical-only (n=245) |
|---|---:|---:|
| `STABLE_ON_HULL` | 233 | 131 |
| `UNSTABLE_NONEXISTENT` | 26 | 96 |
| `METASTABLE_VIABLE` | 8 | 17 |

Collapsing to viable (`STABLE_ON_HULL`+`METASTABLE_VIABLE`) vs. `UNSTABLE_NONEXISTENT`: Fisher's exact test **p=2.8×10⁻¹⁵**, odds ratio **6.01** (χ²=59.9, p=1.0×10⁻¹⁴) — a much larger effect than §2.2's sign-only Fisher test (p=0.0053, odds ratio 1.82). Experimentally-realized compounds are about 6× less likely to be predicted `UNSTABLE_NONEXISTENT` than theoretical-only ones. Of the 188 experimental compounds that are exobondic (§2.2), only 26 are actually predicted `UNSTABLE_NONEXISTENT` — the other 162 are `STABLE_ON_HULL`, since their decomposition into elements, despite being exobondic, is not thermodynamically favorable. Reading §2.2's Fisher table alone as "exobondic = not viable" would be a real misreading `classify.py` itself warns against; §2.3's test is the one that actually asks the viability question the `theoretical` flag is a proxy for, and it is by far the strongest single result in this report.

## 3. Reading

The 7 reference validation cases are still reproduced essentially exactly (§1). §2 now shows a substantially more informative picture than the earlier, smaller population, at every level tested except the finest one: the coarse binary split (has this compound ever been made, yes/no) is a real, significant sign-based discriminator (§2.2), and the full `ViabilityLabel` combining sign with thermodynamics is a *much* stronger one (§2.3, p=2.8×10⁻¹⁵) — the population growth this test was extended for (marginal-formation-energy and max-hull-distance batches, both deliberately chosen for borderline or extreme thermodynamics) supplied exactly the contrast needed for both. The continuous rank correlations against `formation_energy_per_atom` and `energy_above_hull` both weakened at the same time (§2.1), which is the expected price of adding chemically diverse, thermodynamically extreme compounds to what was a narrower population — Spearman correlation is sensitive to exactly this kind of range and composition shift, and is a different kind of test (continuous rank agreement) from the categorical §2.2/§2.3 results, so the two are not in tension. The three-way `family` split (§2.2) has not moved because it is not yet populated for the newer batches (§4) — it is untested at the larger scale, not disconfirmed.

## 4. Limits and next steps

- **Assign `family` labels to the extension/marginal-formation-energy/max-hull batches.** They already carry `theoretical` (used successfully in §2.2–2.3) but not `exp_metastable`/`exp_stable`/`theo_metastable`. Doing so would let the three-way Kruskal-Wallis test draw on the same population growth that made the binary tests significant, and is the most direct way to check whether that finer split would also turn significant.
- **§2.3's `ViabilityLabel` result is new and large — worth a closer look at what's actually driving it**, e.g. whether the odds ratio is itself confounded by `bond_type`/`is_metal` the way §2.1's continuous correlations are (not checked here), or by the fact that max-hull's theoretical polymorphs were deliberately selected to sit far above the hull (a `formation_energy_per_atom`-adjacent selection effect, not necessarily an independent confirmation of the bonding-kinetics mechanism).
- Case 2 (polymorph transition) and case 3 (decomposition into a compound + elements) remain out of scope for this report; see `REPORT_reaction_icohp.md` §6 for case 2's own result.

---

*Code*: `analysis/test_delta_icohp_viability.py`, `analysis/compute_case1_viability.py` (§2.3). *Data*: `analysis/stats_summary_delta_icohp_viability.json`, `analysis/case1_viability.csv`, reusing `analysis/reaction_icohp_case1.csv`, `analysis/reaction_analysis_case1_full.csv`, `mp_dataset/formation_energies.json` (no new DFT calculations for this report). *Validation fixtures*: `tests/test_reitz_dronskowski_validation.py`, `tests/fixtures/reitz_dronskowski_cases.yaml`.
