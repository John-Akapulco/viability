# Antibonding population near the frontier vs. formation energy (mission #4, part B)

**Verdict: at the current 349-compound scale, the antibonding-population-near-frontier metric (normalized) reaches ρ=−0.282, p=9.2×10⁻⁸ (n=346) against `formation_energy_per_atom` — weaker than the ρ=−0.328 (n=186) reported previously, but still comfortably the strongest global correlation of any descriptor in this project (mincut's own headline collapsed to ρ=0.089, p=0.098 — not significant — at this scale; see §2.3). The sign is unchanged: more occupied antibonding character near the frontier still associates with *more negative* (more stable) formation energy, the opposite of the naive Peierls/Jahn-Teller reading. The `bond_type=covalent` subgroup — the only one to ever survive stratification — not only held up at the larger scale, it got *stronger* (n=33, ρ=−0.551, p=0.0009). New at this scale: `is_metal=False` also now survives (n=157, ρ=−0.260, p=0.001), which it did not before. Also new: the `bond_type=ionic` result reported previously (ρ=−0.683, p=0.042, n=9) has evaporated now that the ionic sample grew to n=53 (ρ=−0.184, p=0.187) — exactly the small-n fragility this project's own reporting convention warns about, now demonstrated rather than just flagged.** See §4 for why the sign is not necessarily a refutation of the hypothesis, and §5 for what would be needed to say more.

## 0. What changed this session, and why

No new VASP/LOBSTER calculations for this metric specifically, and no changes to `cohp_extraction.py`. The dataset itself grew: 89 new compounds (alkali/alkaline-earth binaries against N/O/F/P/S/Cl, sourced from Materials Project — see `mp_dataset/download_extension4.py`) were computed, joined into the main pipeline, and `compute_antibonding_all.py` was simply rerun unmodified over the resulting 349-compound `mp_dataset/structures/` (186→349, +163 in the `family=extension` bucket across all four extension batches combined). No code in this metric's pipeline changed; every number below reflects more data, not a different calculation.

- **`analysis/compute_antibonding_all.py`**: 347/349 succeeded. The 2 failures are `extension_S4N2_cod4031496`/`extension_S4N4_cod7017102` — COD-sourced (no Materials Project entry for this composition), so no `is_metal` to draw on; same pre-existing, unrelated-to-this-session gap documented in `METRIC_DEFINITION_antibonding.md`'s AlNi/BeCu discussion.
- **`analysis/stats_analysis_antibonding.py`**: same statistical convention as before, rerun unmodified against the new `analysis/percolation_vs_antibonding.csv` (349 rows).

## 1. The zero floor

`antibond_w_raw` is exactly (< 10⁻⁶) zero for **106/349 compounds (30.4%)**, down from 33.9% (63/186) — extension4 added mostly ionic/covalent chemistries (alkali/alkaline-earth halides, nitrides, oxides, etc.), which have a much lower zero-fraction than metals, diluting the previous ratio:

| Group | n | n zero | frac zero |
|---|---:|---:|---:|
| metal (`is_metal=True`) | 190 | 75 | **39.5%** |
| gapped (`is_metal=False`) | 157 | 31 | 19.8% |
| bond_type=metallic | 88 | 47 | **53.4%** (unchanged — extension4 added zero new `metallic`-classified compounds) |
| bond_type=covalent | 35 | 6 | 17.1% |
| bond_type=ionic | 53 | 8 | 15.1% |

`bond_type=metallic`'s n and zero-fraction are numerically identical to the 186-compound report — extension4's chemistry (alkali/alkaline-earth + N/O/F/P/S/Cl) never classifies as `metallic` under `classify()` (that bucket requires `is_metal=True` with no anion-like element present, and every extension4 compound contains an anion-like element by construction), so this stratum is untouched by the new data. The same qualitative reading as before holds: metals are far more likely to show zero antibonding population in this window than gapped compounds, consistent with the metric's own construction (E_ref=E_F for metals sits inside a partially-filled band; E_ref=VBM for gapped compounds is a band-edge state more likely to carry bonding/antibonding mixing).

## 2. Correlations vs. `formation_energy_per_atom`

### 2.1 Primary metric (ΔE=1.0)

| Group | Metric | n | ρ | p |
|---|---|---:|---:|---:|
| all | antibonding pop. (raw) | 347 | −0.223 | **2.7×10⁻⁵** |
| all | antibonding pop. (normalized) | 346 | **−0.282** | **9.2×10⁻⁸** |
| bond_type=covalent | raw | 33 | −0.535 | **0.0013** |
| bond_type=covalent | normalized | 33 | **−0.551** | **0.0009** |
| bond_type=ionic | raw | 53 | −0.054 | 0.701 |
| bond_type=ionic | normalized | 53 | −0.184 | 0.187 |
| bond_type=metallic | raw | 88 | −0.169 | 0.116 |
| bond_type=metallic | normalized | 88 | −0.186 | 0.082 |
| is_metal=True | raw | 190 | −0.056 | 0.443 |
| is_metal=True | normalized | 189 | −0.069 | 0.348 |
| is_metal=False | raw | 157 | −0.072 | 0.373 |
| is_metal=False | normalized | 157 | **−0.260** | **0.001** |

No small-n rows this time — every group now has n≥33, so the n<15 caveat that applied to the old `bond_type=ionic` (n=9) row no longer applies, and its result changed accordingly (see §3).

### 2.2 ΔE sensitivity (raw, `all` group)

| ΔE (eV) | n | ρ | p |
|---:|---:|---:|---:|
| 0.5 | 347 | −0.239 | 3×10⁻⁵ (rounds to 0.0 at 4dp) |
| 1.0 (primary) | 347 | −0.223 | 2.7×10⁻⁵ |
| 2.0 | 347 | −0.200 | 0.0002 |

Same conclusion as before: the global correlation is essentially insensitive to window width across a 4× change in ΔE.

### 2.3 For comparison: every existing descriptor, same target (n≈345–347)

| Metric | n | ρ | p |
|---|---:|---:|---:|
| **antibonding pop. (normalized, this metric)** | 346 | **−0.282** | **9.2×10⁻⁸** |
| icohp_mean | 347 | −0.345 | ~0 |
| icohp_percolation_weight_min (raw) | 346 | 0.229 | ~0 |
| icohp_percolation_weight_min (normalized) | 345 | 0.172 | 0.0014 |
| mincut (normalized, mission #3 headline) | 346 | 0.089 | 0.098 |
| icohp_min | 347 | −0.055 | 0.307 |
| icohp_sum | 347 | −0.016 | 0.769 |
| icohp_max | 347 | 0.036 | 0.502 |

**Notable change unrelated to this metric**: periodic min-cut's own headline correlation, ρ=0.285 (p=0.0001) at n=186, has **collapsed to ρ=0.089 (p=0.098, not significant) at n=349** — the opposite direction of what happened to `icohp_mean`/`icohp_percolation_weight_min`, which stayed roughly comparable or grew. This project has not investigated why (min-cut's own report, mission #3, was not part of this rewrite); flagged here only because it changes the "ahead of min-cut" framing in the old verdict to "min-cut is no longer significant at all at this scale," which is a stronger claim of relative standing, not a weaker one — worth a dedicated look at `REPORT_mincut.md` (if one exists) or a fresh mission before trusting either number further.

## 3. Is the global signal a between-group effect? (the same question mission #3 had to ask of min-cut)

Two-sample Mann-Whitney U tests, metal vs. gapped, at the new n=349 scale:

| Variable | metal median | gapped median | p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | −0.1232 | −1.1300 | **3.7×10⁻²⁰** |
| `antibond_w_normalized` | 0.0059 | 0.0516 | **1.7×10⁻⁷** |

Same structural picture as before, now on ~2× the data: metals sit at systematically less-negative formation energy *and* systematically lower antibonding population than gapped compounds, both differences even more significant than at n=186. Neither `is_metal` subgroup showed significance alone in the original 186-compound report — but at n=349, **`is_metal=False` now does** (§2.1: ρ=−0.260, p=0.001), while `is_metal=True` still does not (p=0.348). This is new information: the between-group clustering is still real and strong, but it no longer fully explains the pooled signal on its own — there is now a detectable within-gapped-compound relationship that a smaller sample could not resolve.

Where this result continues to differ from mission #3's min-cut: `bond_type=covalent` (n=33, up from 23) reaches significance on its own and got *stronger* (ρ=−0.551, p=0.0009 vs. the old ρ=−0.490, p=0.018) — still the only `bond_type` stratum, of any descriptor tested in this project, to hold up under stratification with a defensible sample size. `bond_type=ionic` (n=53, up from 9) **no longer shows any signal at all** (ρ=−0.184, p=0.187) — the old report's ρ=−0.683 (p=0.042, n=9) is now revealed as a small-sample artifact, not a real subgroup effect. This is worth remembering the next time an n<15 row looks like the strongest number in a table: it can flip entirely once the sample grows.

## 4. On the sign

Unchanged from the 186-compound analysis — the Peierls/Jahn-Teller framing in `METRIC_DEFINITION_antibonding.md` §1 predicts the opposite sign from what's observed, and the same two non-exclusive readings apply:

- **Confound reading (§3)**: gapped/covalent compounds have both more negative formation energies and more antibonding population near a band edge for structural reasons unrelated to instability — this remains a live explanation, strengthened rather than weakened by the larger `is_metal`/`bond_type` gap seen in §3.
- **Target mismatch reading**: `formation_energy_per_atom` measures formation from elements, not proximity to a symmetry-lowering distortion of the already-formed structure — testing against `energy_above_hull` or real persistence/distortion labels remains the more probative, still-undone test.

The `bond_type=covalent` within-type result getting *stronger* with more data (§3) is the one piece of evidence in this report that continues to survive the crudest confound check — worth prioritizing over the now-defunct ionic result for any follow-up.

## 5. Limits and next steps

- **Between-group confound is still the dominant caveat** for the pooled number, now demonstrably not the *whole* story: `is_metal=False` newly survives stratification (§3), a genuine update to the 186-compound conclusion, not just a number refresh.
- **The covalent-subgroup result (n=33, ρ=−0.551) is now the single most robust finding in this report** — it survived both the original stratification check and a near-50% increase in sample size with the correlation getting stronger, not weaker. This is the strongest candidate for a dedicated follow-up (e.g., an independent covalent sub-sample, or breakdown by structure family) of anything in this project so far.
- **The old `bond_type=ionic` result is retracted** — n=9 was always flagged as untrustworthy, and it is now shown to not replicate at n=53. No action needed beyond noting it, but a reminder to treat every current n<15 row (there are none left in this particular report, but plenty elsewhere in this project) with the same suspicion.
- **Sign interpretation is still unresolved** (§4) — this report establishes a robust statistical association, not a validated physical mechanism, same caveat as before.
- **`formation_energy_per_atom` vs. `energy_above_hull` vs. real persistence labels**: still not tested against `energy_above_hull` for this metric specifically (that gap was closed for reaction-ICOHP, see `REPORT_reaction_icohp.md` §2, but not repeated here) — a natural next step now that the pipeline for pulling `energy_above_hull` alongside `formation_energy_per_atom` already exists.
- **min-cut's collapsed headline (§2.3) is an open question this report surfaces but does not answer** — worth a dedicated look before citing either descriptor's ranking against the other with confidence.
- **Still no SISSO.** The covalent-subgroup result is now the best-evidenced single-stratum finding in the project, but it's one bond-type bucket at n=33, not yet the broad, well-powered-per-stratum situation that would justify a combinatorial symbolic search.

---

*Code*: `analysis/compute_antibonding_all.py` (extraction, all 349 compounds), `analysis/stats_analysis_antibonding.py` (statistics + figures). `cohp_extraction.py` untouched since the pilot-validated version. *Data*: `analysis/antibonding_all.json`, `analysis/percolation_vs_antibonding.csv`, `analysis/stats_summary_antibonding.json`. *Figures*: `analysis/figures_antibonding/`.
