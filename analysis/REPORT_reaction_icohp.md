# Reaction ICOHP, decomposition into elements, vs. formation energy / hull distance (mission #5)

**Verdict: `delta_icohp_per_atom` (case 1: AaBb → a A + b B, standard-state elemental references) reaches ρ=0.370, p=3.9×10⁻⁷ (n=177) against `formation_energy_per_atom` — the strongest global correlation found in this project so far, just ahead of the antibonding-population metric (|ρ|=0.337 on the same n=177 subset) and clearly ahead of min-cut (ρ=0.313) and percolation weight (ρ=0.104). The sign is intuitive for once: a compound that is *less* bonding (in ICOHP terms) than its elemental references tends to have a *less negative* (less stable) formation energy. But the same diagnostic this project has now run three times in a row (mission #3 min-cut, mission #4 antibonding, this one) says the same thing: it is substantially a between-group effect. Metal vs. gapped compounds differ enormously in both `delta_icohp_per_atom` (Mann-Whitney p=4.5×10⁻¹⁰) and `formation_energy_per_atom` (p=6.6×10⁻¹¹) — a Kruskal-Wallis across the three `bond_type` strata is even stronger (p=2.1×10⁻⁶ and p=5.0×10⁻⁸ respectively) — and no single `bond_type` stratum reaches significance on its own. Against `energy_above_hull` the correlation vanishes entirely (ρ=0.093, p=0.20, n=190).** See §3 for the stratified breakdown and §5 for the physical reading.

## 0. What changed this session, and why

Continuing directly from `analysis/METRIC_DEFINITION_reaction_icohp.md` §6, which explicitly left two things undone: extending case 1 past its n=1–5 worked examples, and running any statistical test at all.

- **`analysis/compute_reaction_icohp_case1.py` (rerun, not modified)**: at the point it was first run, 24 of the 62 elemental reference calculations (`extension_Sn`, `extension_Se`, `extension_Zr`, `extension_B`, ... `extension_Ru`) were still queued or running on SLURM, so 107/236 compounds were skipped with `reference(s) not yet computed`. All 24 references finished since; rerunning with no code change now succeeds for 192/260 compounds (`ok=192, skipped_single_element=68, skipped_no_reference=0, reference_not_ready=0, failed=0`). The 68 skips are correct exclusions (pure elements decompose into themselves — degenerate reaction, belongs to case 2), not a data gap.
- **`analysis/stats_analysis_reaction_icohp.py` (new)**: same convention as `stats_analysis.py`/`stats_analysis_antibonding.py` — Spearman only, `bond_type`/`is_metal` stratification, n<15 flagged, head-to-head comparison against every prior descriptor on the same target. Two targets tested separately and not assumed to transfer: `energy_above_hull_eV_per_atom` (available for all 192 case-1 rows) and `formation_energy_per_atom` (only for the 177/192 rows whose `mp_id` is in the original 186-compound main campaign — the 15 extension-only compounds, e.g. `Ca3N2`, `Mn2O7`, the TiO2 high-pressure polymorphs, were never sent through the mission #3 formation-energy fetch, and are excluded here rather than imputed). Writes `analysis/stats_summary_reaction_icohp.json` and figures under `analysis/figures_reaction_icohp/`.

`bond_type`/`is_metal` as stored on `reaction_icohp_case1.csv` itself are sparse (184/192 and 183/192 NaN — most of these compounds' `mp_metadata.json` was never run through the main campaign's classification step). Where a compound overlaps the main 186-campaign, the better-populated `bond_type`/`is_metal` from `percolation_vs_formation_energy.csv` / `percolation_vs_antibonding.csv` is used instead. All 15 extension-only rows remain unclassified (`bond_type=NaN`) — they fall into the `all`-group correlations but no stratified one.

## 1. Correlations vs. `formation_energy_per_atom` (n=177, main-campaign overlap)

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 177 | **0.3702** | **3.9×10⁻⁷** |
| bond_type=covalent | 21 | 0.0766 | 0.741 |
| bond_type=ionic$^*$ | 7 | −0.6429 | 0.119 |
| bond_type=metallic | 86 | 0.0966 | 0.376 |
| is_metal=True | 120 | 0.1843 | **0.044** |
| is_metal=False | 57 | 0.2274 | 0.089 |

$^*$ n<15: do not over-interpret in isolation.

Sensitivity (unnormalized `delta_icohp_total`, `all` group): n=177, ρ=0.3345, p=5.3×10⁻⁶ — same sign and similar magnitude to the per-atom version, so the signal is not an artifact of the per-atom normalization.

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

## 4. Comparison against every existing descriptor (same n=177 subset, same target)

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

## 6. Limits and next steps

- **Between-group confound explains the entire headline number** (§3) — worse than mission #3/#4 in this specific respect, since neither prior mission's confound was total; this one's plausibly is. Do not present §1's `all` row without this caveat attached.
- **Case 1 only.** Case 2 (polymorph comparison, no stoichiometric balancing needed) and case 3 (decomposition into a compound + elements) are defined in `reaction_icohp.py` and validated on n=1–5 examples each in the metric definition, but neither has been run at scale or tested statistically. Case 2 in particular (TiO2 polymorphs, carbon allotropes) already showed a real anti-correlation with `energy_above_hull` in the worked examples (§5) that a batch run could test directly, without the elemental-reference machinery this report depends on.
- **`formation_energy_per_atom` coverage gap.** The 15 extension-only compounds (Ca3N2, Mn2O7, the TiO2 high-pressure polymorphs, S4N2/S4N4, ...) have no `formation_energy_per_atom` on file and were excluded from §1/§3/§4 rather than imputed. Fetching it for them (same MP API call pattern as mission #3's `fetch_formation_energy.py`) would let the `all` and `bond_type` correlation tables use the full n=192/n=~207 rather than the 177-row main-campaign-only subset.
- **Elemental reference quality is provisional, not verified.** Per `METRIC_DEFINITION_reaction_icohp.md` §6: 7 of the 62 references were manually overridden against an automated rule-based pick (Ag, In, Rb, Cs, Se, Sn, Ta); the rest carry a "near-degenerate, not individually verified" flag. No LOBSTER `bandOverlaps.lobster` `maxDeviation` quality check has been run on any of the 50 bulk-downloaded element calculations. Two compounds already known to have anomalous band-overlap quality from an earlier mission, `extension_CaN_mp-1058549` and `extension_CaO_mp-2605`, are both present in this 192-row batch (both are compounds under test here, not references, but the same projection-quality concern applies to their own ICOHP data) — their individual `delta_icohp_per_atom` values should be treated as lower-confidence, consistent with the existing project-wide caveat on those two.
- **15/192 rows (all extension-only) remain `bond_type`-unclassified** and only ever enter the `all`-group correlations, never a stratified one — same general caveat as every prior mission's NaN `bond_type` group.
- **Still no SISSO.** Four descriptors now show global signal against `formation_energy_per_atom` (percolation weight the weak exception); this is the strongest by raw ρ but also the one whose signal is most cleanly explained away by between-group clustering — if anything, an argument for stratified/within-group feature search over a pooled one, not yet attempted.

---

*Code*: `analysis/compute_reaction_icohp_case1.py` (extraction), `analysis/stats_analysis_reaction_icohp.py` (statistics + figures), `reaction_icohp.py` (unmodified library, validated in `METRIC_DEFINITION_reaction_icohp.md`). *Data*: `analysis/reaction_icohp_case1.json`, `analysis/reaction_icohp_case1.csv`, `analysis/stats_summary_reaction_icohp.json`. *Figures*: `analysis/figures_reaction_icohp/`.
