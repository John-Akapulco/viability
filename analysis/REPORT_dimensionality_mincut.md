# Network dimensionality + periodic min-cut vs. formation energy (mission #3)

**Verdict: the periodic min-cut descriptor (normalized) is the first metric in this project to show a statistically significant *global* correlation against a stability target — ρ=0.285, p=0.0001 (n=186) against `formation_energy_per_atom` — clearly stronger than the original percolation weight ever achieved against either target (ρ≈0.06-0.11, p≥0.13 throughout). But the signal is likely driven substantially by between-bond-type clustering rather than a uniform within-class relationship (no individual `bond_type` subgroup reaches reliable significance on its own), and network dimensionality alone does *not* separate `formation_energy_per_atom` (Kruskal-Wallis p=0.19).** This is a real, positive, quantified development — but not yet strong enough or clean enough to justify anything beyond continuing to refine these two descriptors (see limits).

## 0. What changed this session, and why

No new VASP/LOBSTER calculations — pure post-processing on the 186 compounds already computed. Two new modules, `percolation_path.py` untouched:

- **`network_dimensionality.py`**: `build_graph()` in `percolation_path.py` uses every bond up to the LOBSTER geometric cutoff (~6 Å), at which radius nearly any compound — including genuinely layered or molecular ones — ends up connected in all 3 directions, so the "disconnected" status essentially never fires and carries no physical dimensionality meaning. This module instead filters bonds to a per-compound relative-strength threshold θ (LobsterPy/George et al., *ChemPlusChem* 2022 convention: θ=10% of the strongest bond) before a plain BFS connectivity check per direction — giving an integer 0D–3D classification.
- **`periodic_mincut.py`**: the existing percolation weight measures *traversability* (cheapest closed walk back to a periodic image) — a compound can have a cheap cycle while still being hard to cleave, if many independent weak paths exist. Min-cut instead measures *separability*: the minimum total bond strength whose removal splits the crystal into two semi-infinite halves along a direction, via a finite "ribbon" graph (nodes = (atom, layer)) and `networkx.minimum_cut` (no hand-rolled flow algorithm).

Both were validated on synthetic cases with known answers **before** running on any real compound (mandatory per the mission brief) — see `tests/test_network_dimensionality.py` (3 cases: isotropic 3D, layered 2D with θ-sensitivity, molecular 0D) and `tests/test_periodic_mincut.py` (4 cases, including one specifically designed to catch an implementation that confuses min-cut with shortest-path: a layered material with 3 independent weak inter-layer bonds must give 3× a single bond's weight, not 1×). All pass; full suite (20 tests total, including the pre-existing 16) is green.

`formation_energy_per_atom` was fetched via a single batched MP query for the 186 `mp_id` already on file (`mp_dataset/fetch_formation_energy.py` → `mp_dataset/formation_energies.json`), 186/186 resolved.

## 1. Network dimensionality distribution

| Dimensionality | n (θ=5%) | n (θ=10%, main) | n (θ=20%) |
|---|---:|---:|---:|
| 0D | 12 | 13 | 18 |
| 1D | 1 | 2 | 5 |
| 2D | 10 | 12 | 16 |
| 3D | 163 | 159 | 147 |

**θ-sensitivity**: 4/186 compounds (2%) change classification between θ=5% and θ=10%; 14/186 (8%) change between θ=10% and θ=20%. The classification is reasonably stable in the 5–10% range and more sensitive going to 20% (expected: a looser threshold keeps more marginal long-range bonds, nudging some 2D/1D/0D compounds up toward 3D). At θ=10%, `bond_type=metallic` is almost entirely 3D (87/88) — physically expected, metals rarely have a genuinely low-dimensional bonding network — while `covalent` and `ionic` show real spread (covalent: 5×0D, 6×2D, 12×3D; ionic: 1×0D, 1×2D, 7×3D), consistent with layered/molecular covalent structures being common and mostly-3D ionic salts having a handful of low-dimensional outliers.

`metastable_covalent_C_rhombohedral_mp-169` (rhombohedral graphite, flagged in the very first pilot report for its strong intra-layer/weak inter-layer anisotropy) comes out **2D at all three θ**, a clean and physically confident result — exactly the behavior the original percolation-weight-only approach could not surface, since at the full 6 Å cutoff graphite's interlayer coupling is nonzero and the compound reads as "connected" in all 3 directions regardless.

## 2. Correlations vs. `formation_energy_per_atom`

Spearman rank correlation (same convention as prior reports: no assumed linearity, small-n groups flagged explicitly).

### 2.1 Min-cut (the new descriptor)

| Group | Metric | n | ρ | p |
|---|---|---:|---:|---:|
| all | mincut (raw) | 186 | 0.261 | **0.0003** |
| all | mincut (normalized) | 186 | **0.285** | **0.0001** |
| bond_type=metallic | mincut (raw) | 88 | 0.148 | 0.168 |
| bond_type=metallic | mincut (normalized) | 88 | 0.143 | 0.183 |
| bond_type=covalent | mincut (raw) | 23 | −0.176 | 0.422 |
| bond_type=ionic$^*$ | mincut (raw) | 9 | 0.500 | 0.171 |
| bond_type=ionic$^*$ | mincut (normalized) | 9 | 0.733 | 0.025 |

$^*$ n<15: do not over-interpret in isolation.

### 2.2 For comparison: the original percolation weight and classic aggregates, same target

| Group | Metric | n | ρ | p |
|---|---|---:|---:|---:|
| all | percolation weight (raw) | 186 | 0.111 | 0.132 |
| all | percolation weight (normalized) | 186 | 0.090 | 0.223 |
| all | icohp_sum | 186 | −0.087 | 0.239 |
| all | icohp_mean | 186 | −0.226 | **0.002** |
| all | icohp_min | 186 | 0.050 | 0.498 |
| all | icohp_max | 186 | 0.205 | **0.005** |

**Reading these together**: `formation_energy_per_atom` is, across the board, a more correlatable target than `energy_above_hull` was — even `icohp_mean` and `icohp_max` reach significance here, which none of the classic aggregates ever did against the hull distance. That's expected: formation energy is a direct cohesion measure against elemental references, mechanically closer to what ICOHP-based descriptors capture, while hull distance is a *relative*-stability measure against competing phases, a different physical question. Against this more favorable target, **min-cut (normalized) is still the single strongest correlation of everything tested** — a genuine, if modest, win for the new descriptor, not just an artifact of an easier target (the percolation weight and icohp_sum, tested against the exact same target, stay non-significant).

**The important caveat**: the figure (`figures_v2/mincut_vs_formation_energy.png`) shows the global trend is visually dominated by *where bond-type clusters sit*, not a uniform slope — metallic compounds (n=88, the majority) sit near `formation_energy≈0` with high, scattered min-cut values (high coordination number → mechanically more bonds per cross-section, independent of true bond strength); ionic/covalent compounds sit further left (more negative formation energy) with uniformly low min-cut. **No individual bond_type subgroup reaches reliable significance on its own** (metallic p=0.17-0.18 at n=88, the only subgroup large enough to trust). This means the strong global number should be read as "min-cut separates bond-type clusters that also happen to differ in formation energy," not yet as "within a given bonding chemistry, min-cut predicts formation energy" — those are different claims, and only the first one is supported here.

## 3. Does dimensionality alone predict formation energy?

| Dimensionality | n | mean (eV/at) | median (eV/at) |
|---|---:|---:|---:|
| 0D | 13 | −0.760 | −0.345 |
| 1D | 2 | −0.548 | −0.548 |
| 2D | 12 | −0.805 | −0.805 |
| 3D | 159 | −0.609 | −0.368 |

Kruskal-Wallis H=4.71, p=0.195 (4 groups). One-way ANOVA F=0.346, p=0.792. **Neither test finds a significant difference in formation energy across dimensionality classes.** Dimensionality on its own is not a useful univariate predictor of `formation_energy_per_atom` in this dataset — group sizes are also very unbalanced (159 vs. 2 for 1D), which limits the power of either test regardless. See `figures_v2/formation_energy_by_dimensionality.png`.

## 4. Limits and next steps

- **Min-cut's global signal is likely a between-bond-type effect, not (yet) a within-type one.** The single most important next step is testing whether this holds up *within* the metallic subgroup specifically once it's the only one with enough points (n=88) to say anything reliable — it currently sits at ρ=0.14-0.15, p≈0.17, i.e. the same "not quite there" territory the original percolation weight was in against `energy_above_hull`.
- **Coordination-number confound for raw min-cut**: absolute min-cut weight mechanically scales with how many bonds cross a given cross-section, which correlates with coordination number/cell composition independent of bond *strength*. The normalized version (÷ strongest bond) partially controls for this but doesn't normalize by bond *count* — a per-bond or per-atom-in-cross-section normalized min-cut is a natural next variant to test, not done here.
- **Sample size per stratum**: 88/23/9 by bond_type once split, same limitation as every prior report in this project. The ionic normalized correlation (ρ=0.733, p=0.025, n=9) is the numerically strongest single result in this whole report and the least trustworthy one.
- **`formation_energy_per_atom` vs. real persistence labels**: this report tests against a second Materials Project thermodynamic quantity, not against actual experimental persistence/decomposition behavior (e.g. the NaBeH₃ vs. CsBeH₃ quenchable/non-quenchable distinction mentioned in the mission brief). That remains the more probative target and the logical next step once (or if) these descriptors show a within-bond-type signal against something DFT-derived first.
- **Still no SISSO.** Two descriptors now show *some* signal (min-cut clearly, dimensionality not at all as a standalone variable) against one target; that is still short of the multi-descriptor, well-powered-per-stratum situation that would justify a combinatorial symbolic search.

## Prioritized next actions

1. Re-test mincut (normalized) correlation restricted to `bond_type=metallic` only, with a larger metallic-only sample if feasible, since that's the one subgroup with both real numbers and enough n to matter.
2. Build a per-bond-count-normalized min-cut variant to check whether the current global signal survives controlling for coordination number.
3. Re-run the dimensionality/min-cut analysis on the conventional-cell pilot outputs from mission #2 (`mp_dataset/structures_conventional_pilot/`, 6 compounds) as a quick cross-check — dimensionality and min-cut are hypothesized to be much less sensitive to the primitive/conventional cell choice than the percolation weight was (min-cut in particular is validated by construction to require cutting *every* parallel path, which is exactly the property that made the percolation weight fragile to cell size in the mission-#2 finding).
4. Scope what a real persistence-label dataset (quenchable/non-quenchable, or a phonon-stability filter as mentioned in the mission-#1 alternatives list) would take to assemble, since that remains the most likely path to a genuinely decision-relevant result.

---

*Data*: `analysis/percolation_vs_formation_energy.csv` (186 rows, all mission-#3 columns joined onto the existing `percolation_vs_hull.csv` columns). *Figures*: `analysis/figures_v2/`. *New code*: `network_dimensionality.py`, `periodic_mincut.py`, `analysis/build_dataset_v2.py`, `mp_dataset/fetch_formation_energy.py` (all at repo root / existing subfolders, `percolation_path.py` and its tests untouched). *Tests*: `tests/test_network_dimensionality.py`, `tests/test_periodic_mincut.py` (20/20 passing project-wide).
