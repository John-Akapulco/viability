# Percolation weight vs. thermodynamic stability — overnight campaign report

**Generated**: 2026-08-14, following the overnight campaign started 2026-08-13 23:56 CEST.

## TL;DR (morning summary)

- **Targeted**: 180 new compounds (60 exp_stable + 60 exp_metastable + 60 theo_metastable,
  unary/binary, non-magnetic, no f-elements, ≤20 atoms/cell) + 6 pilot compounds already
  on disk = **186 total**.
- **Obtained**: **186/186** (100%) have a complete VASP+LOBSTER run and a
  `percolation_path.py` result. Nothing was skipped or left disconnected
  (`n_disconnected_or_missing_percolation = 0`).
- **5 jobs failed on the first pass, all fixed and recovered overnight** (see
  [Failures found and fixed](#failures-found-and-fixed-not-just-logged)) — final
  count above already reflects the recovery.
- **Main quantitative result**: no statistically significant correlation between
  percolation weight (raw or normalized) and `energy_above_hull` **overall**
  (Spearman ρ ≈ 0.06–0.07, p ≈ 0.34–0.38, n=186). A reference logistic regression
  (stable vs. metastable) gets **cross-validated AUC = 0.496 ± 0.141** — chance
  level. Within the metallic subset (n=88) there's a weak, nominally significant
  positive correlation (ρ ≈ 0.19–0.22, p ≈ 0.04–0.08) — but the classic ICOHP
  aggregates show a **comparable-magnitude** correlation in the same subset, so
  percolation does not demonstrate a clear predictive edge over them at this
  sample size. Full numbers in [Correlation results](#correlation-results).
- **Priority actions for next steps** (see [full list](#prioritized-next-actions)):
  1. Extend `metallic` bond_type sample specifically (currently n=88, the only
     subgroup showing a borderline signal) rather than growing the dataset
     uniformly.
  2. Investigate `classify()`'s 66/186 (35%) unclassified compounds (transition
     metal pnictides/chalcogenides, H-containing compounds) — a real coverage
     gap, not analyzed here because the mission brief said reuse `classify()`
     unchanged.
  3. Re-run the conventional-cell question flagged in the original pilot report
     (`report/rapport_fr.tex` §6) — still unaddressed, still the most likely
     confound for the near-zero overall correlation.
  4. `lobster_wall_time_s` is missing for 13/186 compounds (pre-existing pilot
     + individually-resubmitted jobs whose `submit.sh` predates the timing
     instrumentation) — cosmetic gap, not a data-quality issue, but worth
     backfilling if the timing table matters downstream.
  5. Do **not** move to SISSO/symbolic regression yet — see
     [Limits](#limits-and-next-steps) for why.

---

## 1. Dataset

| Family | n | Definition | `energy_above_hull` range (eV/at) |
|---|---:|---|---|
| `exp_stable` | 63 | ICSD-referenced, on the Materials Project convex hull | 0.0000 (by construction) |
| `exp_metastable` | 63 | ICSD-referenced, 0 < E_hull ≤ 0.100 eV/atom | 0.0005 – 0.0907 |
| `theo_metastable` | 60 | Not ICSD-referenced (theoretical-only), 0 < E_hull ≤ 0.200 eV/atom | 0.0013 – 0.1991 |

All unary/binary chemical systems, `|total_magnetization| ≤ 0.01 μ_B`/formula unit
(non-magnetic), no lanthanide/actinide elements, ≤20 sites/cell. Selection code:
`mp_dataset/select_campaign.py`. This is the dataset the user explicitly asked to
keep ("garde les 180 que ta campagne") in place of the mission brief's originally
proposed bond_type × hull-window scheme — the bond_type stratification below is
still produced by reusing `classify()`, just as a secondary axis rather than the
primary sampling axis.

**bond_type** (via `fetch_candidates.py::classify()`, reused unmodified, applied
post hoc from `is_metal` + `elements`):

| bond_type | n |
|---|---:|
| metallic | 88 |
| covalent | 23 |
| ionic | 9 |
| *(unclassified)* | 66 |

35% of the dataset falls outside `classify()`'s three categories — mostly
transition-metal pnictides/chalcogenides (e.g. AsPd₂, FeTe₂, Mo₃Se₄) that are
metallic but contain an element on the "anion-like" exclusion list, and
H-containing compounds (H₂O, CsH, MgH₂ — hydrogen isn't in any of `classify()`'s
element sets). This is a real coverage gap in the heuristic, not a bug in this
analysis; flagged here rather than fixed, per the instruction to reuse `classify()`
as-is.

## 2. Correlation results

Spearman rank correlation with `energy_above_hull_eV_at` (not Pearson: no reason
to expect a linear relationship between a bond-strength descriptor and a
thermodynamic quantity).

### 2.1 Percolation weight (the descriptor under test)

| Group | n | metric | ρ | p | note |
|---|---:|---|---:|---:|---|
| all | 186 | raw | 0.064 | 0.383 | — |
| all | 186 | normalized | 0.071 | 0.339 | — |
| bond_type=covalent | 23 | raw | 0.165 | 0.453 | — |
| bond_type=covalent | 23 | normalized | 0.110 | 0.619 | — |
| bond_type=ionic | 9 | raw | -0.402 | 0.284 | **n<15, do not over-interpret** |
| bond_type=ionic | 9 | normalized | -0.493 | 0.178 | **n<15, do not over-interpret** |
| bond_type=metallic | 88 | raw | 0.191 | 0.075 | borderline |
| bond_type=metallic | 88 | normalized | **0.222** | **0.038** | only nominally-significant result at n≥15 |

Normalization (`weight_min / |icohp_min|`) does not change the qualitative
picture — it nudges the metallic-subset correlation from borderline (p=0.075) to
nominally significant (p=0.038), a small effect, not a step change.

### 2.2 Classic aggregates (for comparison — does percolation add anything?)

| Group | n | metric | ρ | p |
|---|---:|---|---:|---:|
| all | 186 | icohp_sum | 0.032 | 0.670 |
| all | 186 | icohp_mean | -0.058 | 0.436 |
| all | 186 | icohp_min | 0.093 | 0.209 |
| all | 186 | icohp_max | 0.079 | 0.281 |
| bond_type=metallic | 88 | icohp_sum | **0.223** | **0.037** |
| bond_type=metallic | 88 | icohp_min | 0.198 | 0.065 |
| bond_type=ionic | 9 | icohp_sum | -0.694 | 0.038 | n<15 |

**Answer to the question the previous (6-compound) report left open**: the
percolation descriptor is not just *different* from the classic aggregates, it
needs to also be shown *useful* (point 2 of the mission brief). At n=186, in the
one subgroup where anything reaches nominal significance (metallic, n=88),
`icohp_percolation_weight_min_normalized` (ρ=0.222, p=0.038) and `icohp_sum`
(ρ=0.223, p=0.037) are statistically indistinguishable in strength. **No
additional predictive power over the classic aggregates is demonstrated by this
dataset** — a fair, if less exciting, conclusion than the previous report's
"qualitatively different" framing, and the more decision-relevant one.

### 2.3 Reference logistic regression (stable vs. metastable)

Features: `icohp_percolation_weight_min_normalized`, `icohp_mean`, `bond_ionic`,
`bond_metallic` (dummy-coded, `bond_covalent` as reference level). Rows with
unclassified `bond_type` dropped (n=120 of 186).

| | |
|---|---|
| n (stable / metastable) | 120 (45 / 75) |
| CV folds | 5 (StratifiedKFold) |
| **CV AUC (mean ± std)** | **0.496 ± 0.141** |
| per-fold AUC | 0.556, 0.430, 0.252, 0.637, 0.607 |

AUC ≈ 0.5 = coin flip. This 4-feature reference model has no discriminative power
for stable-vs-metastable classification on this dataset. The wide per-fold spread
(0.25–0.64) at n=120 is itself informative: this is not a stable, learnable
signal at this sample size, not just a "weak but real" one.

Figures: `figures/percolation_vs_ehull.png` (scatter, log-scale, colored by
bond_type), `figures/roc_stable_vs_metastable.png` (illustrative full-fit ROC;
the number that matters is the cross-validated AUC above, not this curve).

## 3. Compute cost

| | VASP | LOBSTER |
|---|---:|---:|
| n with timing | 186 | 173 |
| mean (s) | 124 | 397 |
| median (s) | 72 | 145 |
| max (s) | 1440 (24 min) | 6675 (111 min) |
| **total** | **6.4 h** | **19.1 h** |

(Summed *sequential* time; actual wall-clock was ~3.5h thanks to the 8-way
concurrent array.) LOBSTER dominates total compute cost by a factor of ~3 over
VASP — expected, given the orbitalwise COHP/COBI/COOP generation for every
bond up to 6 Å. VASP time scales roughly with cell size (41s mean for ≤4-site
cells vs. 281s mean for 12–20-site cells).

`lobster_wall_time_s` is `NaN` for 13/186 compounds: the 6 pilot compounds
(pre-dating any timing instrumentation) plus 7 individually-resubmitted jobs
(Zn, Ta5Ge3, and the 5 W-containing recoveries below) whose `submit.sh` was
generated by `prepare_vasp_lobster.py`'s per-compound template, which never had
a `time` wrapper around the `lobster-5.1.1` call — only `mp_dataset/submit_array.sh`
did. `vasp_wall_time_s` has no such gap (186/186) because it's parsed from
`vasp.log`, which every job produces regardless of template.

## 4. Failures found and fixed (not just logged)

The mission brief was explicit: *"chaque échec de job doit être loggé avec la
raison probable... regarde OSZICAR/stderr avant de conclure."* Both root causes
below were found this way and are now fixed in the codebase, not just
worked around for tonight's run.

**5/178 jobs failed on the first array pass** (all binary compounds containing
**tungsten**): `Al12W`, `BW`, `WO2`, `TcW`, `TiW`.

- **Root cause**: `mp_dataset/prepare_vasp_lobster.py` overrode W's POTCAR to
  `W_sv` (a "richer" semicore variant) because `W_pv` — the MP/pymatgen
  recommended choice — wasn't present in the local PSP mirror. `POTCAR.W_sv.gz`
  turned out to be **missing its `LEXCH` field entirely** (confirmed by direct
  inspection: `zcat ... | grep LEXCH` returns nothing for `W_sv`, but does for
  `W` and every other element checked). VASP read garbage past the missing
  field and refused to run ("I REFUSE TO CONTINUE WITH THIS SICK JOB").
- **Fix**: `_POTCAR_OVERRIDES = {"W": "W"}` (plain `W`, verified valid PBE
  PAW). 4/5 recovered immediately on retry.

**1/5 (`Al12W`) failed again after the POTCAR fix, with a different error**:
VASP SIGSEGV (stack trace through `__kmpc_fork_call`/OpenMP).

- **Root cause**: exactly the known issue the mission brief flagged
  ("`ulimit -s unlimited` avant `vasp_std` — fix SIGSEGV connu") — but that
  `ulimit` line was present in `mp_dataset/submit_array.sh` (used for the bulk
  178-job array) and had been **forgotten in `prepare_vasp_lobster.py`'s
  per-compound `SLURM_TEMPLATE`**, used for individually-submitted jobs
  (the 6 pilots, smoke tests, and this recovery batch).
- **Fix**: added `ulimit -s unlimited` to `SLURM_TEMPLATE` in
  `prepare_vasp_lobster.py`. Recovered on the next retry.

Net effect: **186/186 final success rate**, both root causes fixed at the
source so they won't recur in any future extension of this dataset.

## 5. Limits and next steps

- **Sample size**: 186 total, but only 88/23/9 per bond_type once split — the
  `ionic` correlations (n=9) are not meaningful and are reported only for
  completeness, flagged inline every time. Even the `metallic` n=88 result
  should be read as "worth extending," not "established."
- **Primitive vs. conventional cell** (carried over from the 6-compound pilot,
  `report/rapport_fr.tex` §6): structures are used in the primitive cell as
  supplied by Materials Project. For high-symmetry structures (rock-salt,
  zinc-blende, etc.) the reported `a/b/c` percolation directions don't
  coincide with conventional crystallographic axes, and the strongest bond
  can sit entirely inside the primitive cell (contributing to no
  non-contractile cycle at all) — this was the headline finding of the pilot
  report and is not re-derived here, but it's the most likely reason the
  correlation with `energy_above_hull` is this weak: the "weakest link" being
  measured may not be the physically meaningful one for a given structure's
  actual failure mode. **This is priority #3** in the actions above precisely
  because it could change the correlation numbers in this report, not just
  add noise reduction.
- **No true "unstable" class**: `theo_metastable` (theoretical-only, up to
  200 meV/atom above hull) is used as a *proxy* for "not persistent," per the
  original framing of this project (a metastable compound never seen
  experimentally may simply not have been made yet, or may be kinetically
  inaccessible — `theoretical=True` in Materials Project conflates both).
  There is no dynamical-stability (phonon) filter here at all — a compound
  could be dynamically unstable and still appear in any of the three families.
- **bond_type coverage**: 66/186 (35%) unclassified by `classify()` (see
  §1) — a real limitation of the reused heuristic, not fixed here per
  instructions.
- **Why not SISSO / symbolic regression yet**: SISSO-class tools are built to
  combine *many* candidate descriptors into interpretable analytic
  expressions. Here there is exactly **one** candidate descriptor (percolation
  weight, raw or normalized) being compared against classic aggregates that
  are themselves algebraically related to it (both are functions of the same
  underlying ICOHP values) — there is no meaningful combinatorial search space
  yet. Justifying a jump to SISSO would require: (a) several more candidate
  descriptors of comparable conceptual novelty (not just more transformations
  of ICOHP), (b) a dataset large enough per stratum (n≥50-100 per bond_type
  at minimum) that a symbolic search wouldn't just overfit noise, and (c) a
  clearer signal than "AUC≈0.5" to suggest there's a learnable relationship
  worth a more expressive search for. None of the three hold yet.

## Prioritized next actions

1. **Extend the `metallic` bond_type sample** (currently the only subgroup
   with a nominally significant, if weak, signal) rather than growing the
   dataset uniformly across all three families.
2. **Investigate the 66 `classify()`-unclassified compounds** — decide whether
   a broader/different bonding-character heuristic is worth building (out of
   scope for tonight per the mission brief, but worth scoping next).
3. **Re-test with conventionalized cells** (`SpacegroupAnalyzer.
   get_conventional_standard_structure`) before graph construction — flagged
   in the pilot report, still open, most likely to move the correlation
   numbers in this report.
4. **Backfill `lobster_wall_time_s`** for the 13 compounds missing it (cosmetic,
   low priority).
5. **Do not start SISSO/symbolic regression** until points 1–3 are addressed
   (see rationale above).

---

*Data*: `analysis/percolation_vs_hull.csv` (186 rows). *Numbers*:
`analysis/stats_summary.json`. *Figures*: `analysis/figures/`. *Pipeline*:
`mp_dataset/select_campaign.py` → `download_campaign.py` →
`prepare_vasp_lobster.py` → `mp_dataset/submit_array.sh` (SLURM, 8-way
capped) → `percolation_path.py` (unmodified) → `analysis/build_dataset.py` →
`analysis/stats_analysis.py`.
