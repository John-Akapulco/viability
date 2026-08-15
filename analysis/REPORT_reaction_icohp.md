# Reaction ICOHP, decomposition into elements, vs. formation energy / hull distance (mission #5)

**Verdict: `delta_icohp_per_atom` (case 1: AaBb → a A + b B, standard-state elemental references) reaches ρ=0.369, p=2.2×10⁻⁷ (n=186, updated 2026-08-15 after extending formation-energy coverage to extension-only compounds — see §1 note) against `formation_energy_per_atom` — the strongest global correlation found in this project so far, just ahead of the antibonding-population metric (|ρ|=0.337 on the shared n=177 main-campaign subset) and clearly ahead of min-cut (ρ=0.313) and percolation weight (ρ=0.104). The sign is intuitive for once: a compound that is *less* bonding (in ICOHP terms) than its elemental references tends to have a *less negative* (less stable) formation energy. But the same diagnostic this project has now run three times in a row (mission #3 min-cut, mission #4 antibonding, this one) says the same thing: it is substantially a between-group effect. Metal vs. gapped compounds differ enormously in both `delta_icohp_per_atom` (Mann-Whitney p=4.5×10⁻¹⁰) and `formation_energy_per_atom` (p=6.6×10⁻¹¹) — a Kruskal-Wallis across the three `bond_type` strata is even stronger (p=2.1×10⁻⁶ and p=5.0×10⁻⁸ respectively) — and no single `bond_type` stratum reaches significance on its own. Against `energy_above_hull` the correlation vanishes entirely (ρ=0.093, p=0.20, n=190).** See §3 for the stratified breakdown and §5 for the physical reading.

## 0. What changed this session, and why

Continuing directly from `analysis/METRIC_DEFINITION_reaction_icohp.md` §6, which explicitly left two things undone: extending case 1 past its n=1–5 worked examples, and running any statistical test at all.

- **`analysis/compute_reaction_icohp_case1.py` (rerun, not modified)**: at the point it was first run, 24 of the 62 elemental reference calculations (`extension_Sn`, `extension_Se`, `extension_Zr`, `extension_B`, ... `extension_Ru`) were still queued or running on SLURM, so 107/236 compounds were skipped with `reference(s) not yet computed`. All 24 references finished since; rerunning with no code change now succeeds for 192/260 compounds (`ok=192, skipped_single_element=68, skipped_no_reference=0, reference_not_ready=0, failed=0`). The 68 skips are correct exclusions (pure elements decompose into themselves — degenerate reaction, belongs to case 2), not a data gap.
- **`analysis/stats_analysis_reaction_icohp.py` (new)**: same convention as `stats_analysis.py`/`stats_analysis_antibonding.py` — Spearman only, `bond_type`/`is_metal` stratification, n<15 flagged, head-to-head comparison against every prior descriptor on the same target. Two targets tested separately and not assumed to transfer: `energy_above_hull_eV_per_atom` (available for all 192 case-1 rows) and `formation_energy_per_atom` (originally only the 177/192 rows in the 186-compound main campaign, extended to 186/192 the same day — see §1's update note). Writes `analysis/stats_summary_reaction_icohp.json` and figures under `analysis/figures_reaction_icohp/`.

`bond_type`/`is_metal` as stored on `reaction_icohp_case1.csv` itself are sparse (184/192 and 183/192 NaN — most of these compounds' `mp_metadata.json` was never run through the main campaign's classification step). Where a compound overlaps the main 186-campaign, the better-populated `bond_type`/`is_metal` from `percolation_vs_formation_energy.csv` / `percolation_vs_antibonding.csv` is used instead. All 15 extension-only rows remain unclassified (`bond_type=NaN`) — they fall into the `all`-group correlations but no stratified one.

## 1. Correlations vs. `formation_energy_per_atom` (n=186)

**Updated 2026-08-15**: `mp_dataset/fetch_formation_energy.py` was rerun (its `collect_mp_ids()` already scans every `structures/*/mp_metadata.json`, no code change needed) and now covers 258/260 compounds instead of only the original 186 — the 2 gaps are `extension_S4N2`/`extension_S4N4` (COD-sourced, no MP `mp_id`, genuinely unfetchable this way). `stats_analysis_reaction_icohp.py` was updated to fall back to this file when `percolation_vs_formation_energy.csv` (main-campaign only) has no value, bringing 9 more extension-only compounds (Ca3N2, Mn2O7, TiO2 polymorphs, ...) into this table. Result barely moved (ρ 0.3702→0.3690) — the 9 new points are consistent with the existing pattern, not a correction to it. 6 rows remain unmatched for an unrelated, pre-existing reason: 4 of the original 6-pilot compounds (NaCl, AlNi, LiBr, BeCu) and the 2 COD compounds have no `mp_id` in `reaction_icohp_case1.csv` itself — not investigated further here.

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 186 | **0.3690** | **2.2×10⁻⁷** |
| bond_type=covalent | 21 | 0.0766 | 0.741 |
| bond_type=ionic$^*$ | 7 | −0.6429 | 0.119 |
| bond_type=metallic | 86 | 0.0966 | 0.376 |
| is_metal=True | 120 | 0.1843 | **0.044** |
| is_metal=False | 57 | 0.2274 | 0.089 |

$^*$ n<15: do not over-interpret in isolation. Stratified rows are unchanged from the original n=177 run — all 9 newly-added points are `bond_type`-unclassified extension compounds, so they only affect the `all` row.

Sensitivity (unnormalized `delta_icohp_total`, `all` group): n=186, ρ=0.3173, p=1.0×10⁻⁵ — same sign and similar magnitude to the per-atom version, so the signal is not an artifact of the per-atom normalization.

## 2. Correlations vs. `energy_above_hull` (n=190, all available rows)

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 190 | 0.0925 | 0.204 |
| bond_type=covalent | 21 | 0.1136 | 0.624 |
| bond_type=ionic$^*$ | 7 | −0.0591 | 0.900 |
| bond_type=metallic | 86 | −0.0560 | 0.609 |
| is_metal=True | 120 | −0.0432 | 0.639 |
| is_metal=False | 57 | 0.0677 | 0.617 |

$^*$ n<15: do not over-interpret in isolation.

Nothing here, anywhere — same pattern as every prior descriptor in this project against this target. `formation_energy_per_atom` and `energy_above_hull` are confirmed once again to be non-interchangeable targets: the headline §1 result does not transfer.

## 3. Is the global formation-energy signal a between-group effect?

Same question mission #3 (min-cut, `bond_type`) and mission #4 (antibonding, `is_metal`) had to ask, and the same answer:

| Variable | metal median | gapped median | Mann-Whitney p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | −0.291 | −0.930 | **6.6×10⁻¹¹** |
| `delta_icohp_per_atom` | 2.188 | 0.139 | **4.5×10⁻¹⁰** |

| Variable | Kruskal-Wallis across `bond_type` (covalent/ionic/metallic) |
|---|---:|
| `formation_energy_per_atom` | p=5.0×10⁻⁸ |
| `delta_icohp_per_atom` | p=2.1×10⁻⁶ |

Bond-type medians of `delta_icohp_per_atom`: metallic +2.39, covalent +0.15, ionic −0.15 — metals sit far above gapped compounds on this metric, and simultaneously have far less negative formation energies (general chemistry: metals form with less exothermic formation energy than strongly-bonded covalent/ionic compounds). This mechanically produces a positive pooled correlation even absent any true within-group relationship, exactly the same structural situation as mission #3/#4. Consistent with that reading: **no `bond_type` stratum reaches significance in §1** (best is `is_metal=True` at p=0.044, weak and unlike the antibonding metric's covalent-subgroup result which held at p<0.03). Unlike mission #4, this metric has **no surviving subgroup** — the between-group confound explains the entirety of the observable signal here, not just most of it.

## 4. Comparison against every existing descriptor (n=177 main-campaign-only subset, same target)

Restricted to `in_main_campaign` rows specifically (unlike §1's n=186) because the comparison descriptors (percolation weight, min-cut, antibonding) were only ever computed for the original 186-compound campaign.

| Metric | n | ρ | p |
|---|---:|---:|---:|
| **reaction ICOHP, per atom (this mission)** | 177 | **0.3702** | **3.9×10⁻⁷** |
| antibonding population (normalized, mission #4 headline) | 177 | −0.3369 | 4.5×10⁻⁶ |
| periodic min-cut (normalized, mission #3 headline) | 177 | 0.3126 | 2.3×10⁻⁵ |
| percolation weight (normalized, mission #1) | 177 | 0.1037 | 0.170 |

Numerically the strongest of the four on this common subset — but it is also the one with the *weakest* within-group survival (§3), which is the opposite of what would make a headline number trustworthy. Rank by raw |ρ| is not the same as rank by "how much of this is real."

## 5. Physical reading

`METRIC_DEFINITION_reaction_icohp.md` §5 already flagged this from three hand-worked examples before any statistics were run: Ca3N2 → 3 Ca + N2 and Mn2O7 → MnO2 + O2 both showed strongly *positive* `delta_icohp_per_atom` (less bonding than elemental references) for compounds that are in fact strongly ionic and thermodynamically very stable — because **ICOHP measures orbital-overlap bond population, not electrostatic/Madelung lattice energy**, and ionic compounds derive most of their formation stability from the latter. The `bond_type=ionic` row in §1 (ρ=−0.643, n=7, not significant) is numerically consistent with that warning — small sample, but the sign flips relative to the pooled result, exactly where the metric's own definition said to expect trouble. The carbon-allotrope and TiO2-polymorph worked examples in the same document showed a second failure mode (van der Waals / weak interlayer bonding invisible to ICOHP) that this batch run cannot speak to directly, since case 2 (polymorph comparison) was not run at scale here.

Unlike the antibonding metric (mission #4), whose sign was *counterintuitive* relative to its own motivating hypothesis, this metric's sign is the intuitive one — less-bonding-than-elements correlates with less stable formation energy, no unresolved sign puzzle. The problem here is not interpretation but that the signal doesn't survive stratification at all.

## 6. Case 2 at scale: does bonding track stability within a fixed composition?

**Added 2026-08-15**, via `analysis/compute_reaction_icohp_case2.py` (new, uses `reaction_icohp.compare_polymorphs` unmodified). Case 2 needs no elemental references or reaction balancing — any group of ≥2 compounds sharing a reduced formula is directly comparable via `icohp_per_atom`. Only **8 such groups exist across the full 260-compound `structures/` tree** (21 compound-dirs total) — this dataset was built for chemical-system breadth (`select_campaign.py` caps 2 entries per chemsys), not polymorph density, so this is a small, exploratory pass, not a powered test.

| Formula | n | Most-bonding | Most-stable (min EAH) | Agree? |
|---|---:|---|---|:---:|
| C | 6 | lonsdaleite | rhombohedral (metastable_covalent) | No |
| MgSn | 2 | mp-1094801 | mp-1094801 | Yes |
| PbSe | 2 | exp_stable (mp-2201) | exp_stable (mp-2201) | Yes |
| TaP | 2 | theo_metastable | exp_stable | No |
| TiNi | 2 | exp_stable | exp_stable | Yes |
| TiO2 | 3 | baddeleyite | TiO2-II | No |
| TiSi2 | 2 | theo_metastable | exp_stable | No |
| Zn | 2 | theo_metastable (mp-2647117) | extension (mp-79) | No |

**3/8 agree** — at or slightly below the ≈3.5/8 expected by chance alone for this mix of group sizes (six 2-way, one 3-way, one 6-way: $\sum 1/\text{size} = 6\times\tfrac12+\tfrac13+\tfrac16=3.5$). This confirms, at the only scale this dataset currently supports, what the hand-worked carbon-allotrope and TiO2-polymorph examples in `METRIC_DEFINITION_reaction_icohp.md` §5 already showed individually: **`icohp_per_atom` does not track which polymorph is thermodynamically preferred.** Consistent with the same physical reading as §5 — ICOHP sees covalent orbital overlap only, missing the van der Waals and electrostatic contributions that often decide between structurally similar polymorphs (graphite's weak interlayer bonding vs. diamond/lonsdaleite's strong covalent network is the clearest case: diamond and lonsdaleite dominate `icohp_per_atom` but graphite, not either of them, is the true ground state within the C group here too).

One footnote: the `Zn` group's theoretical member (`theo_metastable_Zn_mp-2647117`) is the same compound `compute_reaction_icohp_case1.py`'s own comments already flagged as "wrong structure type, not a valid elemental reference" for case 1 — that concern doesn't apply here (case 2 just compares two real structures head-to-head), but it's worth knowing the two caveats share a compound.

*Code*: `analysis/compute_reaction_icohp_case2.py`. *Data*: `analysis/reaction_icohp_case2.json`, `analysis/reaction_icohp_case2.csv`.

## 7. Limits and next steps

- **Between-group confound explains the entire headline number** (§3) — worse than mission #3/#4 in this specific respect, since neither prior mission's confound was total; this one's plausibly is. Do not present §1's `all` row without this caveat attached.
- **Case 2: done at scale (§6), but severely sample-limited.** Only 8 polymorph groups exist in the current 260-compound dataset — 3/8 agreement (at/below chance) confirms the worked-example finding but cannot be called a statistical test. Growing this would require deliberately sampling multiple polymorphs per chemical system, which cuts directly against `select_campaign.py`'s chemsys-breadth criterion — a real tension between this project's two sampling goals, not something fixable by just computing more compounds from the existing selection logic.
- **Case 3 (decomposition into a compound + elements) not attempted at scale, and likely not tractable without new DFT.** Unlike case 1 (products are always pure elements, and 62 elemental references now cover most of the periodic table present in this dataset) or case 2 (any same-formula pair already in the dataset works), case 3 needs the *specific* lower-energy compound(s) a given compound would decompose into — generally not itself in this project's structure set unless hand-picked (as Mn2O7 → MnO2 + O2 was). Materials Project's own phase-diagram machinery (`pymatgen.analysis.phase_diagram`, via `PDEntry`/`PhaseDiagram`) can compute the correct equilibrium decomposition products for any composition programmatically — but those products would then need their own VASP+LOBSTER calculations if not already present, i.e. this is a new calculation campaign, not a batch script over existing data like cases 1 and 2 were.
- **`formation_energy_per_atom` coverage gap: resolved 2026-08-15** (§1 note) — 9 more extension-only compounds now included via `fetch_formation_energy.py`, n=177→186. 6 rows remain unmatched (2 COD compounds with no `mp_id`, 4 pre-existing pilot-compound gaps unrelated to this fix) and §4's comparison table is still main-campaign-only (n=177) since the comparison descriptors themselves were never computed for extension compounds — that part of the gap is unchanged.
- **Elemental reference quality: checked 2026-08-15, and it is worse than expected.** `analysis/check_lobster_quality_extension.py` parsed `bandOverlaps.lobster` `maxDeviation` for all 74 extension compounds: **27/62 elemental references (44%) exceed 0.1** (threshold chosen well above the well-behaved-calculation range, <0.03, seen throughout the 6-pilot dataset). 26 of the 28 flagged extension compounds overall are elemental references actually used in `ELEMENT_REFERENCE` — Ir 4.88, Hf 3.78, Os 3.56, Ag 3.14, Ti 2.39, Au 1.75, Cd 1.62, Co 1.14, Ni 1.02, Pd 0.96, Tc 0.94, Li 0.90, Zn 0.76, Sc 0.73, Pt 0.70, Mo 0.67, Re 0.67, Ru 0.66, Cr 0.62, C-graphite 0.55, V 0.47, Rh 0.34, B 0.28, Fe 0.23, Ta 0.16, and W at 15.98 (comparable in magnitude to `extension_CaN_mp-1058549`/`extension_CaO_mp-2605`, 17.25/15.44, the two compounds-under-test already known to have this issue). The pattern — almost entirely transition metals — is consistent with under-converged k-meshes for metals with complex Fermi surfaces, not isolated data errors. **This means most reactions in the §1 result use at least one flagged reference, which was not accounted for in that correlation.** The obvious next check — restricting §1's correlation to reactions built entirely from unflagged references — was run immediately after this quality scan: 104/192 rows (54%) use at least one flagged reference. Restricted to the 78 clean-reference rows with `formation_energy_per_atom` available, ρ=0.392, p=3.8×10⁻⁴ — *stronger*, not weaker, than the pooled §1 result; the 99-row flagged-reference subset alone is weaker (ρ=0.233, p=0.020). **The LOBSTER band-overlap problem does not appear to be inflating the headline correlation — if anything it dilutes it.** This does not retroactively validate the flagged references' own ICOHP values (still provisional), but it rules out "the correlation is an artifact of bad references" as an explanation for §1's result specifically. Separately, 7 of the 62 references were manually overridden against an automated rule-based pick (Ag, In, Rb, Cs, Se, Sn, Ta — see `METRIC_DEFINITION_reaction_icohp.md` §6); Ag is also band-overlap-flagged, compounding both caveats on the same reference. Full results: `analysis/lobster_quality_extension.json`.
- **15/192 rows (all extension-only) remain `bond_type`-unclassified** and only ever enter the `all`-group correlations, never a stratified one — same general caveat as every prior mission's NaN `bond_type` group.
- **Still no SISSO.** Four descriptors now show global signal against `formation_energy_per_atom` (percolation weight the weak exception); this is the strongest by raw ρ but also the one whose signal is most cleanly explained away by between-group clustering — if anything, an argument for stratified/within-group feature search over a pooled one, not yet attempted.

---

*Code*: `analysis/compute_reaction_icohp_case1.py` (extraction), `analysis/stats_analysis_reaction_icohp.py` (statistics + figures), `analysis/check_lobster_quality_extension.py` (band-overlap quality scan), `mp_dataset/fetch_formation_energy.py` (rerun 2026-08-15, unmodified), `reaction_icohp.py` (unmodified library, validated in `METRIC_DEFINITION_reaction_icohp.md`). *Data*: `analysis/reaction_icohp_case1.json`, `analysis/reaction_icohp_case1.csv`, `analysis/stats_summary_reaction_icohp.json`, `analysis/lobster_quality_extension.json`, `mp_dataset/formation_energies.json`. *Figures*: `analysis/figures_reaction_icohp/`.
