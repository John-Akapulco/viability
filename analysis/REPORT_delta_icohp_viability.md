# Does Δ(ICOHP) distinguish stable, metastable, and unstable compounds?

This is the project's central research question: does Δ(ICOHP) for the decomposition-into-elements reaction (`reaction_icohp.py`/`reaction_analysis`, the sign of Δ(ICOHP)/atom) discriminate thermodynamically stable, metastable, and unstable compounds? Computed by `analysis/test_delta_icohp_viability.py`:

**Extended to every element and compound with a computed case-1 reaction (518 reactions, pooled without regard to which batch a compound came from — see §1), the coarse experimental-vs-theoretical-only split is a significant discriminator.** Fisher's exact test on Δ(ICOHP)-sign vs. experimentally-realized/theoretical-only is **p=0.0053** (odds ratio 1.82), and Mann-Whitney on the continuous Δ(ICOHP)/atom, same split, is **p=0.0088** — a Δ(ICOHP) ≥ 0 sign (breaking the reactant's bonds costs more than the products recover) is markedly more common, and Δ(ICOHP)/atom markedly less negative, among compounds Materials Project records as experimentally realized. The finer three-way split (`exp_metastable`/`exp_stable`/`theo_metastable`) remains a real, monotonic, but not-yet-significant trend (Kruskal-Wallis p=0.278, n=177 — unchanged from before, since `family` labels are only assigned within the original campaign, not the newer batches; see §3). Against continuous stability targets, Δ(ICOHP)/atom correlates with `formation_energy_per_atom` (ρ=−0.297, p=7.3×10⁻¹²) — weaker than on the earlier, smaller population (ρ=−0.50) but still highly significant — and no longer even borderline with `energy_above_hull` (ρ=−0.058, p=0.19, up from p=0.068).

## 1. Extended to every compound with a computed case-1 reaction

Every compound with a computed case-1 (decomposition-into-elements) reaction is pooled into one population, regardless of which batch produced it — batch boundaries reflect submission order, not compound chemistry, and splitting by them would reintroduce exactly the kind of between-group confound this project's methodology is built to catch. Current population: **518 reactions**.

### 1.1 Continuous stability targets

| Target | n | ρ | p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | 511 | **−0.2969** | **7.3×10⁻¹²** |
| `energy_above_hull_eV_per_atom` | 515 | −0.0579 | 0.190 |

(Sign convention: `reaction_analysis`'s own products−reactants, i.e. the *opposite* of `reaction_icohp.py`'s reactant−products column reported elsewhere in this project — same magnitude, flipped sign; see `delta.py`'s docstring.) Both correlations are weaker than on the earlier, smaller population: the formation-energy result was ρ=−0.50 before the marginal-formation-energy and max-hull-distance batches were folded in, and the hull-distance result has moved from borderline (p=0.068) to clearly non-significant (p=0.19). This is consistent with those batches being deliberately selected for marginal or extreme thermodynamics rather than typical compounds — they add real chemical diversity but dilute a rank correlation that was partly carried by the original population's narrower stability range.

### 1.2 Sign-based discrimination

**Experimental-viability proxy** (`theoretical` flag: `False` = experimentally realized, `True` = theoretical-only per Materials Project):

| | Δ(ICOHP) ≥ 0 | Δ(ICOHP) < 0 |
|---|---:|---:|
| experimental (n=267) | 79 | 188 |
| theoretical-only (n=245) | 46 | 199 |

Fisher's exact test: **p=0.0053**, odds ratio=1.82 — a Δ(ICOHP) ≥ 0 sign is significantly more common among experimentally-realized compounds than theoretical-only ones, the expected direction (this sign is meant as a kinetic-barrier signature that lets a compound persist despite favorable decomposition thermodynamics). Mann-Whitney on the continuous Δ(ICOHP)/atom itself, same split: experimental median −0.868 eV, theoretical-only median −1.457 eV (more negative), **p=0.0088** — same direction, also significant.

**Finer-grained (`family`, three-way, n=177 — unchanged from the original campaign; the extension/marginal-formation-energy/max-hull batches carry a `theoretical` flag but not a `family` label, so this specific test has not grown):**

| family | n | median Δ(ICOHP)/atom (eV) | % with Δ(ICOHP) ≥ 0 |
|---|---:|---:|---:|
| exp_metastable | 59 | −1.038 | **28.8%** |
| exp_stable | 59 | −1.391 | 22.0% |
| theo_metastable | 59 | −1.601 | **13.6%** |

Kruskal-Wallis across the three groups: H=2.56, p=0.278 — still not significant, but the ordering is the same physically-motivated one as before: `exp_metastable` compounds have the highest Δ(ICOHP) ≥ 0 fraction and the least-negative median Δ(ICOHP), `theo_metastable` compounds have the lowest fraction and the most-negative median.

### 1.3 Full `ViabilityLabel` (not just the sign) — now large enough to test

`classify_viability()` combines `BondingLabel`'s sign with `delta_energy` (`-formation_energy_per_atom`): a compound is only `UNSTABLE_NONEXISTENT` if decomposition into elements is *both* Δ(ICOHP)-negative *and* thermodynamically favorable, a strictly stronger, more specific claim than the sign-only split in §1.2. Run via `analysis/compute_case1_viability.py` on the same pooled, no-batch-exclusion population (518 reactions, 511 classifiable — 7 excluded for missing `formation_energy_per_atom`, all COD-sourced or otherwise mp_id-less compounds):

| | experimental (n=267) | theoretical-only (n=245) |
|---|---:|---:|
| `STABLE_ON_HULL` | 233 | 131 |
| `UNSTABLE_NONEXISTENT` | 26 | 96 |
| `METASTABLE_VIABLE` | 8 | 17 |

Collapsing to viable (`STABLE_ON_HULL`+`METASTABLE_VIABLE`) vs. `UNSTABLE_NONEXISTENT`: Fisher's exact test **p=2.8×10⁻¹⁵**, odds ratio **6.01** (χ²=59.9, p=1.0×10⁻¹⁴) — a much larger effect than §1.2's sign-only Fisher test (p=0.0053, odds ratio 1.82). Experimentally-realized compounds are about 6× less likely to be predicted `UNSTABLE_NONEXISTENT` than theoretical-only ones. Of the 188 experimental compounds with Δ(ICOHP) < 0 (§1.2), only 26 are actually predicted `UNSTABLE_NONEXISTENT` — the other 162 are `STABLE_ON_HULL`, since their decomposition into elements, despite the negative sign, is not thermodynamically favorable. Reading §1.2's Fisher table alone as "Δ(ICOHP) < 0 = not viable" would be a real misreading `classify.py` itself warns against; §1.3's test is the one that actually asks the viability question the `theoretical` flag is a proxy for, and it is by far the strongest single result in this report.

## 2. Reading

§1 shows a substantially more informative picture than the earlier, smaller population, at every level tested except the finest one: the coarse binary split (has this compound ever been made, yes/no) is a real, significant sign-based discriminator (§1.2), and the full `ViabilityLabel` combining sign with thermodynamics is a *much* stronger one (§1.3, p=2.8×10⁻¹⁵) — the population growth this test was extended for (marginal-formation-energy and max-hull-distance batches, both deliberately chosen for borderline or extreme thermodynamics) supplied exactly the contrast needed for both. The continuous rank correlations against `formation_energy_per_atom` and `energy_above_hull` both weakened at the same time (§1.1), which is the expected price of adding chemically diverse, thermodynamically extreme compounds to what was a narrower population — Spearman correlation is sensitive to exactly this kind of range and composition shift, and is a different kind of test (continuous rank agreement) from the categorical §1.2/§1.3 results, so the two are not in tension. The three-way `family` split (§1.2) has not moved because it is not yet populated for the newer batches (§3) — it is untested at the larger scale, not disconfirmed.

## 3. Limits and next steps

- **Assign `family` labels to the extension/marginal-formation-energy/max-hull batches.** They already carry `theoretical` (used successfully in §1.2–1.3) but not `exp_metastable`/`exp_stable`/`theo_metastable`. Doing so would let the three-way Kruskal-Wallis test draw on the same population growth that made the binary tests significant, and is the most direct way to check whether that finer split would also turn significant.
- **§1.3's `ViabilityLabel` result is new and large — worth a closer look at what's actually driving it**, e.g. whether the odds ratio is itself confounded by `bond_type`/`is_metal` the way §1.1's continuous correlations are (not checked here), or by the fact that max-hull's theoretical polymorphs were deliberately selected to sit far above the hull (a `formation_energy_per_atom`-adjacent selection effect, not necessarily an independent confirmation of the bonding-kinetics mechanism).
- Case 2 (polymorph transition) and case 3 (decomposition into a compound + elements) remain out of scope for this report; see `REPORT_reaction_icohp.md` §6 for case 2's own result.
- A parallel test of the same viability question using Δ(ICOHP antibonding)/Δ(ICOBI antibonding) instead of total Δ(ICOHP) — mission #4b, see `REPORT_antibonding.md` — is now the project's strongest single correlation and its own viability discrimination (Mann-Whitney p=3.7×10⁻¹², Fisher p=2.0×10⁻⁷); the two descriptors' viability calls have not yet been cross-tabulated against each other to see how much they agree case-by-case.

---

*Code*: `analysis/test_delta_icohp_viability.py`, `analysis/compute_case1_viability.py` (§1.3). *Data*: `analysis/stats_summary_delta_icohp_viability.json`, `analysis/case1_viability.csv`, reusing `analysis/reaction_icohp_case1.csv`, `analysis/reaction_analysis_case1_full.csv`, `mp_dataset/formation_energies.json` (no new DFT calculations for this report).
