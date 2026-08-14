# Antibonding population near the frontier vs. formation energy (mission #4, part B)

**Verdict: the antibonding-population-near-frontier metric (normalized) reaches ρ=−0.328, p=5.0×10⁻⁶ (n=186) against `formation_energy_per_atom` — the strongest global correlation of any descriptor in this project so far, ahead of periodic min-cut's ρ=0.285, p=0.0001. The sign says compounds with *more* occupied antibonding character right at the frontier tend to have *more negative* (more stable) formation energy — the opposite direction the Peierls/Jahn-Teller motivation naively predicts. Diagnostics show this global number is substantially, but not entirely, a between-group (metal vs. gapped) effect: within the `is_metal` split neither group reaches significance alone, but within `bond_type=covalent` (n=23) the correlation *does* survive (p=0.018–0.029) — a genuinely different pattern from min-cut, where no bond_type subgroup held up at all.** See §4 for why the sign is not necessarily a refutation of the hypothesis, and §5 for what would be needed to say more.

## 0. What changed this session, and why

No new VASP/LOBSTER calculations, no changes to `cohp_extraction.py` or `percolation_path.py`. Two things, continuing directly from `analysis/METRIC_DEFINITION_antibonding.md` (validated on 6 pilots only, extension explicitly left open):

- **`analysis/compute_antibonding_all.py`**: runs `cohp_extraction.antibonding_population_near_frontier` on all 186 compounds already in `mp_dataset/structures/`, at ΔE ∈ {0.5, 1.0, 2.0} eV, with `is_metal` fetched from Materials Project (single batched query, same convention as the fetch for `formation_energy_per_atom` in mission #3 — not derived locally, per the AlNi/BeCu pitfall documented in the metric definition). All 186/186 succeeded (no `antibond_error` rows). Writes `analysis/antibonding_all.json` (full per-compound detail) and `analysis/percolation_vs_antibonding.csv` (186×44, primary ΔE=1.0 columns plus `antibond_w_raw_dE{0.5,1.0,2.0}` for sensitivity, merged onto every prior descriptor already on file: percolation weight, ICOHP/ICOBI aggregates, dimensionality, min-cut, `formation_energy_per_atom`).
- **`analysis/stats_analysis_antibonding.py`**: same statistical convention as `stats_analysis.py` (Spearman, no assumed linearity, n<15 groups flagged) and the same report structure as mission #3 — correlation table (`all` / `bond_type` / **`is_metal`**, the last one new here because the metric's own definition treats metals and gapped compounds differently, E_ref=E_F vs E_ref=VBM), a ΔE-sensitivity check, and a head-to-head comparison against every existing descriptor on the same target. Writes `analysis/stats_summary_antibonding.json` and figures under `analysis/figures_antibonding/`.

## 1. The zero floor

Before reading correlation numbers: `antibond_w_raw` is exactly (< 10⁻⁶) zero for **63/186 compounds (33.9%)** — i.e. no occupied antibonding COHP at all within 1 eV of the frontier for a third of the dataset. This is not evenly spread:

| Group | n | n zero | frac zero |
|---|---:|---:|---:|
| metal (`is_metal=True`) | 125 | 56 | **44.8%** |
| gapped (`is_metal=False`) | 61 | 7 | 11.5% |
| bond_type=metallic | 88 | 47 | 53.4% |
| bond_type=covalent | 23 | 2 | 8.7% |
| bond_type=ionic | 9 | 1 | 11.1% |

Metals are far more likely to show zero antibonding population in this window than gapped compounds. This matches the metric's own construction — for a metal E_ref=E_F sits inside a partially-filled band where whether the immediately-sub-Fermi states carry net antibonding character is compound-specific and often small, whereas a gapped compound's VBM is a band-edge state that more often carries some bonding/antibonding mixing. It also means correlations computed on the metallic subgroup are testing a variable that is zero more often than not — worth keeping in mind reading §2.

## 2. Correlations vs. `formation_energy_per_atom`

### 2.1 Primary metric (ΔE=1.0)

| Group | Metric | n | ρ | p |
|---|---|---:|---:|---:|
| all | antibonding pop. (raw) | 186 | −0.312 | **1.5×10⁻⁵** |
| all | antibonding pop. (normalized) | 186 | **−0.328** | **5.0×10⁻⁶** |
| bond_type=covalent | raw | 23 | −0.455 | **0.029** |
| bond_type=covalent | normalized | 23 | −0.490 | **0.018** |
| bond_type=ionic$^*$ | raw | 9 | −0.433 | 0.244 |
| bond_type=ionic$^*$ | normalized | 9 | −0.683 | **0.042** |
| bond_type=metallic | raw | 88 | −0.169 | 0.116 |
| bond_type=metallic | normalized | 88 | −0.186 | 0.082 |
| is_metal=True | raw | 125 | −0.117 | 0.192 |
| is_metal=True | normalized | 125 | −0.126 | 0.163 |
| is_metal=False | raw | 61 | −0.100 | 0.443 |
| is_metal=False | normalized | 61 | −0.237 | 0.065 |

$^*$ n<15: do not over-interpret in isolation.

### 2.2 ΔE sensitivity (raw, `all` group)

| ΔE (eV) | n | ρ | p |
|---:|---:|---:|---:|
| 0.5 | 186 | −0.321 | 8.1×10⁻⁶ |
| 1.0 (primary) | 186 | −0.312 | 1.5×10⁻⁵ |
| 2.0 | 186 | −0.309 | 1.7×10⁻⁵ |

The global correlation is essentially insensitive to the window width — sign, magnitude, and significance are stable across a 4× change in ΔE. This rules out the signal being an artifact of the specific 1.0 eV cutoff.

### 2.3 For comparison: every existing descriptor, same target

| Metric | n | ρ | p |
|---|---:|---:|---:|
| **antibonding pop. (normalized, this metric)** | 186 | **−0.328** | **5.0×10⁻⁶** |
| mincut (normalized, mission #3 headline) | 186 | 0.285 | 0.0001 |
| icohp_mean | 186 | −0.226 | 0.002 |
| icohp_max | 186 | 0.205 | 0.005 |
| percolation weight (raw) | 186 | 0.111 | 0.132 |
| percolation weight (normalized) | 186 | 0.090 | 0.223 |
| icohp_sum | 186 | −0.087 | 0.239 |
| icohp_min | 186 | 0.050 | 0.498 |

The antibonding-population metric is now the single strongest correlate of `formation_energy_per_atom` found in this project, ahead of min-cut.

## 3. Is the global signal a between-group effect? (the same question mission #3 had to ask of min-cut)

Two-sample Mann-Whitney U tests, metal vs. gapped:

| Variable | metal median | gapped median | p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | −0.272 | −0.930 | **9.4×10⁻¹¹** |
| `antibond_w_normalized` | 0.0015 | 0.0516 | **1.2×10⁻⁷** |

Metals sit at systematically less-negative formation energy *and* systematically lower antibonding population than gapped compounds — both differences highly significant, and in the direction that would mechanically produce a negative pooled correlation even with zero true within-group relationship. This is the same structural situation as min-cut in mission #3 (there: between-`bond_type` clustering; here: between-`is_metal` clustering), and neither `is_metal` subgroup reaches significance alone (metal p=0.16–0.19, gapped p=0.065–0.44), so **the global number should not be read as "this metric predicts formation energy within a given electronic character."**

Where this result differs from mission #3's: `bond_type=covalent` (n=23, entirely gapped compounds by construction of the bond-type labeling) reaches significance on its own (p=0.018–0.029) — the only subgroup, of any descriptor tested so far in this project including min-cut, to hold up under stratification at a defensible sample size. `bond_type=ionic` (n=9) shows the numerically strongest normalized correlation (ρ=−0.683, p=0.042) but is too small to trust in isolation, same caveat as every n<15 group in every prior report.

## 4. On the sign

The Peierls/Jahn-Teller framing in `METRIC_DEFINITION_antibonding.md` §1 predicts that *more* occupied antibonding character near the frontier signals *closer to an electronic instability* — naively, less stable, i.e. a less negative (higher) `formation_energy_per_atom`. The observed correlation is negative: more antibonding population associates with *more* negative (more stable) formation energy. Two non-exclusive readings, neither tested further here:

- **Confound reading (§3)**: gapped/covalent compounds happen to have both more negative formation energies (general chemistry: strong covalent/ionic bonds form very exothermically) and more antibonding population near a band edge (VBM states are the top of a bonding-to-antibonding-mixed manifold by construction) — the correlation could be picking up "this compound is covalent/gapped" rather than anything about electronic instability specifically.
- **Target mismatch reading**: `formation_energy_per_atom` measures energy released forming the compound from *elements*, not proximity to a symmetry-lowering distortion of the *already-formed* structure. A compound can be a very exothermic formation product (very negative formation energy) and still sit close to a Peierls-type instability in its adopted structure — these are different physical questions, same limitation flagged for min-cut in mission #3 and for the original percolation weight in mission #1. The metric was motivated by structural instability, not thermodynamic formation stability; testing it against `energy_above_hull` (relative stability among the compound's own polymorphs/decomposition products) or, better, real persistence/distortion labels, is the more probative test and was not done here.

The `bond_type=covalent` within-type result (§3) at least survives the crudest version of the confound check (it holds within a single bond-type stratum), so it is not purely an artifact of pooling covalent and metallic compounds together — but it is one subgroup at n=23, not yet a validated finding.

## 5. Limits and next steps

- **Between-group confound is the dominant caveat**, same structural issue as mission #3's min-cut result — now for `is_metal` rather than `bond_type`. The covalent-subgroup result (§3) is the one piece of this report that survives stratification and is worth following up specifically, e.g. checking whether it holds on an independent covalent sub-sample or breaks down further by structure family.
- **Sign interpretation is unresolved** (§4) — this report establishes a robust statistical association, not a validated physical mechanism. No compound in this 186-set is independently documented as a Peierls/Jahn-Teller case, so nothing here confirms or refutes the instability hypothesis specifically; it tests correlation with thermodynamic formation stability, a related but distinct question.
- **`formation_energy_per_atom` vs. `energy_above_hull` vs. real persistence labels**: same limitation noted in every report since mission #1. `energy_above_hull` was the project's original target; this report (like mission #3) used `formation_energy_per_atom` because it is more correlatable in general, but that also means the two targets are not interchangeable and this result should not be assumed to hold against hull distance without checking.
- **Zero floor (§1)**: 34% of compounds have exactly-zero raw antibonding population in this window; for the metallic subgroup that's 53%. A metric that is a point mass at zero for over half a stratum has limited discriminating power within that stratum regardless of correlation testing — consistent with the metallic subgroup being the weakest (p=0.08–0.12) of the three bond types tested.
- **Still no SISSO.** Three descriptors now show *some* signal against `formation_energy_per_atom` (min-cut, and now antibonding population more strongly), but the strongest one here is also the one with the clearest known confound and an unresolved sign — not yet the clean, well-powered-per-stratum situation that would justify a combinatorial symbolic search.

---

*Code*: `analysis/compute_antibonding_all.py` (extraction, all 186 compounds), `analysis/stats_analysis_antibonding.py` (statistics + figures). `cohp_extraction.py` untouched from the pilot-validated version. *Data*: `analysis/antibonding_all.json`, `analysis/percolation_vs_antibonding.csv`, `analysis/stats_summary_antibonding.json`. *Figures*: `analysis/figures_antibonding/`.
