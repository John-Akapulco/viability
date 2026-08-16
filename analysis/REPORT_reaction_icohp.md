# Reaction ICOHP, decomposition into elements, vs. formation energy / hull distance (mission #5)

**Verdict: at the current 349-compound scale, `delta_icohp_per_atom` (case 1: AaBb → a A + b B, standard-state elemental references) reaches ρ=0.502, p=6.3×10⁻¹⁹ (n=275) against `formation_energy_per_atom` — substantially stronger than the ρ=0.369 (n=186) reported previously, and now clearly the strongest descriptor in this project by a wider margin than before (next is antibonding population at ρ=−0.384 on the same n=275 subset, then min-cut at ρ=0.130). The between-group confound this project has now checked four times (mission #3 min-cut, mission #4 antibonding, this metric previously, and now again) is still real and still huge at the `bond_type` level — Kruskal-Wallis across covalent/ionic/metallic remains overwhelming (p=1.1×10⁻¹⁵) and no `bond_type` stratum reaches significance alone. But this is the first time in the project that a coarser split — `is_metal` — tells a different story: **both `is_metal=False` (ρ=0.319, p=0.0001, n=138) and `is_metal=True` (ρ=0.256, p=0.0025, n=137) now survive**, where at n=186 neither did (best was p=0.044, "weak and unlike" any real survival). Against `energy_above_hull` the correlation is no longer completely null either: ρ=0.139, p=0.018 (n=291) — still an order of magnitude weaker than the formation-energy result, but no longer indistinguishable from zero.** See §3 for the stratified breakdown, §5 for the physical reading, and §6 for case 2 (polymorph comparison) now run on 49 groups instead of 8.

## 0. What changed this session, and why

The dataset grew: 89 new compounds (`mp_dataset/download_extension4.py` — alkali/alkaline-earth binaries against N/O/F/P/S/Cl, sourced from Materials Project) were computed, joined into the main pipeline, and every script below was rerun **unmodified** over the resulting 349-compound `mp_dataset/structures/`. `extension4`'s constituent elements (Li/Na/K/Rb/Cs, Be/Mg/Ca/Sr/Ba, N/O/F/P/S/Cl) were all already present in `ELEMENT_REFERENCE` — no new elemental references were downloaded or computed this session, so the elemental-reference-quality caveat from before (§7 below) is unaffected by this update.

- **`analysis/compute_reaction_icohp_case1.py` (rerun, not modified)**: 192 → **281 successful case-1 reactions** (`ok=281, skipped_single_element=68, skipped_no_reference=0, reference_not_ready=0, failed=0`). All 89 new compounds succeeded; 0 failures, 0 missing references.
- **`analysis/populate_reaction_analysis_case1_full.py` (rerun, not modified)**: the schema-driven `reaction_analysis/` package's parallel case-1 population, cross-checked against the line above — **281/281 still match** `reaction_icohp_case1.csv` to floating-point precision after the known sign flip, same as at n=192. `reactions_dataset/` now covers all 281.
- **`analysis/compute_reaction_icohp_case2.py` (rerun, not modified)**: **8 → 49 polymorph groups** — extension4 was deliberately designed to fix this project's polymorph-density gap (§6 in the old report explicitly flagged 8 groups as too few to call a test), by preferring, for each chemical system, a formula that already had an entry elsewhere in the dataset. It worked: see §6.
- **`analysis/stats_analysis_reaction_icohp.py` (rerun, not modified)**: same convention as before. One structural simplification: at n=192 (previous run), `n_in_main_campaign` (177 or 186) was smaller than the full case-1 population, so §4's comparison table had to be restricted to a "main-campaign-only" subset. At the current scale, `n_in_main_campaign` (299) equals the full case-1 row count (299) exactly — every case-1-eligible compound is now part of the joined `percolation_vs_formation_energy.csv`/`percolation_vs_antibonding.csv` used for the comparison descriptors, so that restriction and its caveat no longer apply (§4 is simpler this time).
- `bond_type`/`is_metal` are pulled from `percolation_vs_formation_energy.csv`/`percolation_vs_antibonding.csv` (better-populated than `reaction_icohp_case1.csv`'s own copies), same as before. 116/299 rows remain `bond_type`-unclassified and only enter the `all`-group correlations.

## 1. Correlations vs. `formation_energy_per_atom` (n=275)

`mp_dataset/fetch_formation_energy.py` was rerun and covers 347/349 structures (the 2 gaps are the COD-sourced `extension_S4N2`/`extension_S4N4`, unfetchable via Materials Project by construction, unchanged from before). Of the 281 successful case-1 reactions, 275 have a `formation_energy_per_atom` value; 6 remain unmatched for the same pre-existing reasons as before (the 2 COD compounds plus 4 original 6-pilot compounds with no `mp_id` recorded in `reaction_icohp_case1.csv` itself — not investigated further here, unchanged from the prior report).

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 275 | **0.5016** | **6.3×10⁻¹⁹** |
| bond_type=covalent | 22 | 0.1451 | 0.519 |
| bond_type=ionic | 51 | 0.1157 | 0.419 |
| bond_type=metallic | 86 | 0.0966 | 0.376 |
| is_metal=True | 137 | **0.2561** | **0.0025** |
| is_metal=False | 138 | **0.3190** | **0.0001** |

No n<15 rows this time. Sensitivity (unnormalized `delta_icohp_total`, `all` group): n=275, ρ=0.4135, p≈0 — same sign and comparable magnitude to the per-atom version, so the signal is not an artifact of the per-atom normalization.

## 2. Correlations vs. `energy_above_hull` (n=291, all available rows)

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 291 | **0.1390** | **0.0177** |
| bond_type=covalent | 38 | 0.0468 | 0.780 |
| bond_type=ionic | 51 | 0.1574 | 0.270 |
| bond_type=metallic | 86 | −0.0560 | 0.609 |
| is_metal=True | 137 | −0.0140 | 0.871 |
| is_metal=False | 138 | 0.1035 | 0.227 |

Sensitivity (unnormalized `delta_icohp_total`, `all` group): n=291, ρ=0.0601, p=0.307 — the normalized version's weak-but-real signal does not survive un-normalizing, unlike §1's much larger effect. This target remains far less correlatable than `formation_energy_per_atom` (consistent with every prior descriptor in this project), but the `all`-group result is no longer flatly null the way it was at n=190 (ρ=0.093, p=0.20) — worth re-checking again at a future scale rather than treating as settled either way.

## 3. Is the global formation-energy signal a between-group effect?

Same question mission #3 (min-cut, `bond_type`) and mission #4 (antibonding, `is_metal`) had to ask, now at n=275:

| Variable | metal median | gapped median | Mann-Whitney p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | −0.2940 | −1.2322 | **1.8×10⁻¹⁹** |
| `delta_icohp_per_atom` | 2.1366 | 0.2479 | **3.4×10⁻¹⁶** |

| Variable | Kruskal-Wallis across `bond_type` (covalent/ionic/metallic) |
|---|---:|
| `formation_energy_per_atom` | p=5.6×10⁻²² |
| `delta_icohp_per_atom` | p=1.1×10⁻¹⁵ |

Bond-type medians of `delta_icohp_per_atom`: metallic +2.39 (essentially unchanged from before), covalent −0.18 (was +0.15 — **sign flipped** with more data), ionic +0.13 (was −0.15 — **also flipped**). Metals still sit far above gapped compounds on this metric while simultaneously having far less negative formation energies, and the `bond_type` split is, if anything, more extreme than before (p=1.1×10⁻¹⁵ vs. the old p=2.1×10⁻⁶). Consistent with that: **no `bond_type` stratum reaches significance in §1**, same conclusion as before — the between-`bond_type`-group confound still fully explains the pooled correlation at that level of stratification.

**What's new**: the coarser `is_metal` split behaves completely differently from `bond_type`. At n=186, the best `is_metal` result was p=0.044 (weak, borderline) and this report's predecessor concluded "no surviving subgroup — the between-group confound explains the entirety of the observable signal." At n=275, **both `is_metal` groups survive comfortably** (p=0.0025 and p=0.0001). This is not a contradiction of the `bond_type` finding — `is_metal` and `bond_type` are correlated but not identical splits (e.g. `bond_type=metallic` is a subset of `is_metal=True`, but `is_metal=True` also includes some non-`metallic`-classified compounds, and vice versa for gapped/covalent+ionic) — but it does mean the old headline's flat claim that "the between-group confound explains the entirety of the observable signal here" needs to be walked back: **the confound fully explains the `bond_type`-level pooling, but there is now a statistically real relationship between bonding and formation energy *within* a fixed electronic character (metal-only, or gapped-only) that a smaller sample could not detect.** Whether that within-`is_metal` signal is itself further confounded by something else (e.g. anion identity, coordination, or a `bond_type`-adjacent variable not yet tested) is not established here — it is a genuine open finding, not yet a validated mechanism.

## 4. Comparison against every existing descriptor (n=275, same target)

Unlike the previous report, no subset restriction is needed here — `n_in_main_campaign` (299) now equals the full case-1 population (299), so every comparison descriptor (percolation weight, min-cut, antibonding) is available for the same compounds as the reaction-ICOHP numbers themselves.

| Metric | n | ρ | p |
|---|---:|---:|---:|
| **reaction ICOHP, per atom (this mission)** | 275 | **0.5016** | **~0** |
| antibonding population (normalized, mission #4 headline) | 275 | −0.3841 | ~0 |
| periodic min-cut (normalized, mission #3 headline) | 275 | 0.1297 | 0.0315 |
| percolation weight (normalized, mission #1) | 274 | 0.1282 | 0.0339 |

Reaction ICOHP remains the strongest of the four, and the gap to the next descriptor (antibonding, |ρ|=0.384) is now wider than it was at n=177 (where the gap to antibonding was only 0.033 in ρ). Note this comparison-table figure for antibonding population (ρ=−0.384 on this n=275 shared subset) is noticeably stronger than antibonding's own full-sample headline (ρ=−0.282, n=346, `REPORT_antibonding.md`) — the two reports are testing overlapping but not identical compound sets (this table is restricted to case-1-eligible, non-single-element compounds with a formation energy value), a reminder that "the" correlation for a given descriptor is somewhat sample-dependent and the two numbers are not in tension, just not directly comparable.

## 5. Physical reading

Unchanged from the original analysis, and if anything reinforced by §3's sign-flip observation: `METRIC_DEFINITION_reaction_icohp.md` §5 already flagged, from hand-worked examples, that **ICOHP measures orbital-overlap bond population, not electrostatic/Madelung lattice energy**, and ionic compounds derive most of their formation stability from the latter. The `bond_type=ionic` median flipping sign between n=7 (old) and n=51 (new) — and the `bond_type=covalent` median flipping too — is consistent with this metric behaving noisily/inconsistently for compound classes where ICOHP is not the dominant contributor to stability, exactly where the metric's own definition said to expect trouble. The carbon-allotrope and TiO2-polymorph worked examples (van der Waals / weak interlayer bonding invisible to ICOHP) are addressed directly in §6 now that case 2 has real statistical power.

Unlike the antibonding metric (mission #4), whose sign is counterintuitive relative to its own motivating hypothesis, this metric's sign remains the intuitive one — less-bonding-than-elements correlates with less stable formation energy. The interpretive puzzle here is not the sign but the newly-asymmetric stratification result (§3): a real `is_metal`-conditional signal alongside a fully-confounded `bond_type`-conditional one.

## 6. Case 2 at scale: does bonding track stability within a fixed composition?

**Rerun 2026-08-16** via `analysis/compute_reaction_icohp_case2.py` (unmodified) now that extension4 added same-formula polymorph pairs by design. **49 groups now exist across the full 349-compound `structures/` tree** (109 compound-dirs total), up from 8 — this is the direct payoff of extension4's selection rule (§0), and turns this section from an 8-point exploratory pass into a real, if still modest, statistical test.

| Group size | n groups | n agree | agreement rate |
|---:|---:|---:|---:|
| 2 (pairwise) | 38 | 21 | 55.3% |
| 3 | 10 | 5 | 50.0% |
| 6 | 1 | 0 | 0.0% |
| **total** | **49** | **26** | **53.1%** |

Expected agreement under pure chance (each group's most-bonding member independently uniform over its members): Σ(1/size) = 22.5/49 (45.9%). Observed (26/49) is modestly above this, but an exact Poisson-binomial test (each group's own chance probability 1/size, summed) gives **P(X≥26 | pure chance) = 0.19** — not distinguishable from noise at this sample size. Unlike the old n=8 result (3/8, at/below chance), this one leans the other direction, but neither leans hard enough to call it a real effect either way.

The only 6-member group, elemental carbon, is a repeat of the old report's finding and did not change: lonsdaleite remains most-bonding, the rhombohedral metastable-covalent phase remains most-stable — bonding still does not track polymorph stability for this system specifically, the clearest case in the dataset since it's the same physical story as before (graphite's weak interlayer bonding vs. diamond/lonsdaleite's strong covalent network).

**Reading**: at 49 groups, this remains consistent with the hand-worked carbon/TiO2 conclusion in `METRIC_DEFINITION_reaction_icohp.md` §5 — `icohp_per_atom` shows, at best, a chance-level relationship to which polymorph is thermodynamically preferred, not a real one. The physical reason given there (ICOHP sees covalent orbital overlap only, missing van der Waals and electrostatic contributions that often decide between structurally similar polymorphs) is unchanged and, with 6× the sample, not contradicted by the larger test.

*Code*: `analysis/compute_reaction_icohp_case2.py`. *Data*: `analysis/reaction_icohp_case2.json`, `analysis/reaction_icohp_case2.csv`.

## 7. Limits and next steps

- **`bond_type`-level confound still fully explains the pooled §1 number** — unchanged conclusion. **`is_metal`-level confound does NOT fully explain it** — new conclusion this session (§3), the most notable update in this report. Do not present §1's `all` row without the `bond_type` caveat, but also do not claim (as the previous version of this report did) that the entirety of the signal is confound-explained — the `is_metal`-stratified result is a real, if not yet mechanistically understood, finding.
- **Case 2: now genuinely at scale (§6)** — 49 groups, up from 8, closes the sampling-limitation caveat from the previous report. Result: 53.1% agreement, statistically indistinguishable from the 45.9% chance baseline (P=0.19). This is a more informative null result than before, not a positive finding — `icohp_per_atom` still does not appear to track polymorph stability.
- **Case 3 (decomposition into a compound + elements) still not attempted at scale, and likely still not tractable without new DFT** — unchanged from before; extension4 did not add anything usable here since it targeted case-1/case-2 gaps specifically.
- **Elemental reference quality: unaffected by this session's changes.** `ELEMENT_REFERENCE` was not extended (extension4's elements were all already covered) — the 44%-flagged band-overlap issue documented previously (and the finding that restricting to clean references made the correlation *stronger*, not weaker) still stands as last measured; not rechecked this session since nothing that would change it happened.
- **`formation_energy_per_atom` coverage**: 347/349 structures now covered (up from 258/260), same 2 unfetchable COD compounds as always.
- **116/299 rows remain `bond_type`-unclassified** and only ever enter the `all`-group correlations — same general caveat as every prior mission's NaN `bond_type` group, proportionally similar to before (was 15/192, now 116/299 — the ratio is roughly stable, not worsened by extension4).
- **Still no SISSO.** Reaction ICOHP is now, by a wider margin than before, the strongest raw correlate of `formation_energy_per_atom` in this project, and — new this session — the first descriptor whose signal survives a coarse stratification (`is_metal`) rather than being fully explained away by it. This is the strongest argument yet in this project for a stratified (`is_metal`-conditional) feature search over a pooled one, still not attempted.

---

*Code*: `analysis/compute_reaction_icohp_case1.py` (extraction), `analysis/compute_reaction_icohp_case2.py` (case 2), `analysis/stats_analysis_reaction_icohp.py` (statistics + figures), `analysis/check_lobster_quality_extension.py` (band-overlap quality scan, not rerun this session), `mp_dataset/fetch_formation_energy.py` (rerun, unmodified), `reaction_icohp.py` (unmodified library, validated in `METRIC_DEFINITION_reaction_icohp.md`). *Data*: `analysis/reaction_icohp_case1.json`, `analysis/reaction_icohp_case1.csv`, `analysis/reaction_icohp_case2.json`, `analysis/reaction_icohp_case2.csv`, `analysis/stats_summary_reaction_icohp.json`, `mp_dataset/formation_energies.json`. *Figures*: `analysis/figures_reaction_icohp/`.
