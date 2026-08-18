# Reaction ICOHP, decomposition into elements, vs. formation energy / hull distance (mission #5)

**Verdict: at the current 518-reaction scale, `delta_icohp_per_atom` (case 1: AaBb → a A + b B, standard-state elemental references) reaches ρ=0.297, p≈0 (n=511) against `formation_energy_per_atom`.** This is still a real, highly significant correlation given the sample size, but weaker than it looked at smaller scale, and the confound picture has changed shape: the `bond_type`-level split is **no longer uniformly confound-explained** — `bond_type=ionic` (ρ=0.360, p=0.0004) and `bond_type=mixed` (ρ=−0.525, p≈0, opposite sign) both now survive on their own, `bond_type=covalent` also survives but flips sign (ρ=−0.354, p=0.023), and only `bond_type=metallic` stays non-significant (ρ=0.096, p=0.089). The coarser `is_metal` split — which used to be the one that survived stratification — is now essentially redundant with `bond_type=metallic` (same n, same ρ, same p in every row below) and does **not** clear significance on its own (`is_metal=True` p=0.089, `is_metal=False` p=0.122). Against `energy_above_hull` the correlation stays weak and, at the `all`-group level, is no longer even borderline (ρ=0.058, p=0.190). See §3 for the stratified breakdown, §4 for how this metric now compares to the other three descriptors on a shared subset (antibonding edges it out there), §5 for the physical reading, and §6 for case 2 (polymorph comparison) at 74 groups.

## 1. Correlations vs. `formation_energy_per_atom` (n=511)

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 511 | **0.2969** | **≈0** |
| bond_type=covalent | 41 | **−0.3539** | **0.023** |
| bond_type=ionic | 93 | **0.3604** | **0.0004** |
| bond_type=metallic | 316 | 0.0959 | 0.089 |
| bond_type=mixed | 61 | **−0.5249** | **≈0** |
| is_metal=True | 316 | 0.0959 | 0.089 |
| is_metal=False | 195 | −0.1111 | 0.122 |

Sensitivity (unnormalized `delta_icohp_total`, `all` group): n=511, ρ=0.2095, p≈0 — same sign, smaller magnitude than the per-atom version, so most (not all) of the signal survives un-normalizing.

## 2. Correlations vs. `energy_above_hull` (n=515, all available rows)

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 515 | 0.0579 | 0.190 |
| bond_type=covalent | 41 | −0.2933 | 0.063 |
| bond_type=ionic | 95 | 0.0556 | 0.592 |
| bond_type=metallic | 318 | −0.0698 | 0.214 |
| bond_type=mixed | 61 | **−0.3073** | **0.016** |
| is_metal=True | 318 | −0.0698 | 0.214 |
| is_metal=False | 197 | −0.1265 | 0.077 |

Sensitivity (`delta_icohp_total`, `all` group): n=515, ρ=−0.0562, p=0.203. This target remains far less correlatable than `formation_energy_per_atom`, consistent with every prior descriptor in this project — the one exception is `bond_type=mixed`, the only stratum that reaches significance against either target here.

## 3. Is the global formation-energy signal a between-group effect?

Same question mission #3 (min-cut, `bond_type`) and mission #4 (antibonding, `is_metal`) had to ask, now with four `bond_type` categories (a "mixed"/Zintl category was added alongside covalent/ionic/metallic):

| Variable | metal median | gapped median | Mann-Whitney p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | −0.0689 (n=316) | −0.7609 (n=195) | **1.1×10⁻²⁵** |
| `delta_icohp_per_atom` | 2.1881 (n=319) | 0.0199 (n=197) | **1.7×10⁻³⁹** |

| Variable | Kruskal-Wallis across `bond_type` (covalent/ionic/metallic/mixed) |
|---|---:|
| `formation_energy_per_atom` | H=153.5, p=4.7×10⁻³³ |
| `delta_icohp_per_atom` | H=179.0, p=1.5×10⁻³⁸ |

Bond-type medians of `delta_icohp_per_atom`: metallic +2.19 (least-bonding-relative-to-elements, by far the largest group), ionic +0.15, covalent −0.36, mixed −0.67 (the most exobondic group on average). Metals still sit far above gapped compounds on this metric while simultaneously having far less negative formation energies, and both group-level tests remain overwhelming — pooling across `bond_type` is still not a valid way to read the `all` row in §1.

**What §1 shows within that caveat**: unlike the `bond_type`-vs-`is_metal` asymmetry found at a smaller scale, three of the four `bond_type` strata now carry their own significant, non-trivial signal — ionic and mixed strongly, covalent more weakly and with the opposite sign from ionic. `bond_type=metallic` (316 of 511 rows, the majority) does not. Because `is_metal=True` now maps almost exactly onto `bond_type=metallic` (identical n, ρ, and p in every row of §1–§2 above — the two splits have converged to the same population at this dataset composition), `is_metal` no longer tells a different story from `bond_type`: it simply inherits `bond_type=metallic`'s non-significant result. The most informative single fact in this section is not "does the confound explain everything" (yes, at the pooled level) but that **the sign of the surviving relationship differs by bond type** (positive for ionic, negative for covalent and mixed) — a genuine open finding about what "more/less bonding than the elements" means differently across bonding character, not yet a validated mechanism.

## 4. Comparison against every existing descriptor (n=288, shared subset)

Restricted to the subset with all four descriptors available (`n_in_main_campaign`=288 of 518 — the comparison CSVs this table depends on, `percolation_vs_formation_energy.csv`/`percolation_vs_antibonding.csv`, do not yet cover every batch in the full case-1 population).

| Metric | n | ρ (vs `formation_energy_per_atom`) | p | ρ (vs `energy_above_hull`) | p |
|---|---:|---:|---:|---:|---:|
| antibonding population (mission #4) | 288 | **−0.3940** | ≈0 | 0.1407 | 0.017 |
| **reaction ICOHP, per atom (this mission)** | 288 | 0.3704 | ≈0 | 0.0970 | 0.100 |
| periodic min-cut (mission #3) | 288 | 0.1176 | 0.046 | 0.0543 | 0.359 |
| percolation weight (mission #1) | 287 | 0.1055 | 0.074 | 0.0835 | 0.158 |

On this shared subset, antibonding population's `|ρ|` (0.394) edges out reaction ICOHP's (0.370) against `formation_energy_per_atom` — the two are close, and reaction ICOHP is the only one of the four that stays inside the same order of magnitude against `energy_above_hull` too, but the plain "reaction ICOHP is the strongest descriptor in the project" statement no longer holds cleanly on this particular subset. It remains the strongest against `energy_above_hull` here (0.097 vs. antibonding's 0.141 — antibonding is actually stronger on hull-distance too). Min-cut and percolation weight trail both by a wide margin on every target, consistent with every prior report.

## 5. Physical reading

`METRIC_DEFINITION_reaction_icohp.md` §5 flagged, from hand-worked examples, that **ICOHP measures orbital-overlap bond population, not electrostatic/Madelung lattice energy** — ionic compounds derive most of their formation stability from the latter. §3's bond-type medians (ionic positive, covalent and mixed negative) and §1's now-significant, opposite-signed correlations within `bond_type=ionic` vs. `bond_type=covalent`/`mixed` are consistent with this: the metric is not a uniform proxy for stability across bonding character, it behaves differently — and in the ionic case, apparently still usefully — depending on which physical mechanism actually carries the formation energy for that class of compound.

This metric's sign remains the intuitive one throughout (less-bonding-than-elements correlates with less stable formation energy, within groups where the correlation exists at all) — unlike the antibonding metric (mission #4), whose sign is counterintuitive relative to its own motivating hypothesis. The interpretive puzzle here is the opposite-signed strata (§3), not the sign itself.

## 6. Case 2: does bonding track stability within a fixed composition?

`analysis/compute_reaction_icohp_case2.py` finds every reduced formula with ≥2 computed polymorphs and asks whether the most-bonding member (least `icohp_per_atom`... i.e. most negative) is also the most-stable one (lowest `energy_above_hull`). **74 groups** with ≥2 hull-labeled members exist in the current dataset (163 compound-dirs total).

| Group size | n groups | n agree | agreement rate |
|---:|---:|---:|---:|
| 2 (pairwise) | 48 | 28 | 58.3% |
| 3 | 14 | 5 | 35.7% |
| 4 | 8 | 2 | 25.0% |
| 5 | 1 | 0 | 0.0% |
| 6 | 2 | 0 | 0.0% |
| 7 | 1 | 0 | 0.0% |
| **total** | **74** | **35** | **47.3%** |

Expected agreement under pure chance (each group's most-bonding member independently uniform over its members): Σ(1/size) = 31.3/74 (42.4%). Observed (35/74) sits slightly above this, but an exact Poisson-binomial test (each group's own chance probability 1/size, summed) gives **P(X≥35 | pure chance) = 0.22** — not distinguishable from noise.

The two 6-member groups are elemental carbon and NaN₃; the one 7-member group is C₃N₄. For carbon specifically — the clearest physical story in the dataset, same finding as at smaller scale — lonsdaleite remains most-bonding while the metastable rhombohedral covalent phase remains most-stable: bonding does not track polymorph stability for this system, the same graphite/diamond-family story (weak interlayer bonding vs. strong covalent network) `METRIC_DEFINITION_reaction_icohp.md` §5 already worked out by hand.

**Reading**: at 74 groups, this remains consistent with that hand-worked conclusion — `icohp_per_atom` shows, at best, a chance-level relationship to which polymorph is thermodynamically preferred, not a real one. ICOHP sees covalent orbital overlap only, missing van der Waals and electrostatic contributions that often decide between structurally similar polymorphs.

*Code*: `analysis/compute_reaction_icohp_case2.py`. *Data*: `analysis/reaction_icohp_case2.json`, `analysis/reaction_icohp_case2.csv`.

## 7. Limits and next steps

- **The `bond_type`/`is_metal` confound picture has genuinely changed shape (§3) and needs re-reading, not just re-numbering.** At smaller scale, `bond_type` was fully confound-explained and `is_metal` was the split that survived; now `is_metal` has converged onto `bond_type=metallic` (identical population) and inherited its null result, while three of four `bond_type` strata (ionic, covalent, mixed) carry their own significant, opposite-signed signal. Any future write-up of this metric should stratify by `bond_type`, not `is_metal` — `is_metal` no longer adds information beyond `bond_type=metallic` at this dataset composition.
- **Reaction ICOHP is no longer unambiguously the strongest descriptor in the project (§4)** — antibonding population edges it out on the shared 288-compound subset against both targets. Worth re-running §4's comparison once the comparison CSVs (`percolation_vs_formation_energy.csv`/`percolation_vs_antibonding.csv`) cover the full 518-row population rather than the 288-row overlap subset, since the ranking could plausibly change again at full scale.
- **Case 2: 74 groups, still a chance-level null result (§6)** — same conclusion as at 8 and 49 groups, now with much more statistical power behind it. Treat `icohp_per_atom` as not informative for polymorph-stability ranking; this looks settled rather than under-sampled at this point.
- **Case 3 (decomposition into a compound + elements) remains not attempted at scale** — unchanged, still likely not tractable without new DFT targeted at that reaction topology specifically.
- **`bond_type` classification now has 0 unclassified rows** (was 116/299) — the ICOBI-first-shell classifier (`analysis/compute_icohp_icobi_bondtype.py`) now covers every case-1 compound with LOBSTER data, closing the coverage caveat every earlier version of this report carried.
- **Still no SISSO.** Given §3's finding that the useful signal is bond-type-conditional and opposite-signed across strata, a pooled feature search would be actively misleading here — any future search should be run separately per `bond_type`, not on the pooled population.

---

*Code*: `analysis/compute_reaction_icohp_case1.py` (extraction), `analysis/compute_reaction_icohp_case2.py` (case 2), `analysis/stats_analysis_reaction_icohp.py` (statistics + figures), `reaction_icohp.py` (library, validated in `METRIC_DEFINITION_reaction_icohp.md`), `mp_dataset/fetch_formation_energy.py`. *Data*: `analysis/reaction_icohp_case1.json`, `analysis/reaction_icohp_case1.csv`, `analysis/reaction_icohp_case2.json`, `analysis/reaction_icohp_case2.csv`, `analysis/stats_summary_reaction_icohp.json`, `mp_dataset/formation_energies.json`. *Figures*: `analysis/figures_reaction_icohp/`.
