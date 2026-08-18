# Antibonding population near the frontier vs. formation energy (mission #4, part B)

**Verdict: at the current 588-compound scale (Cu4Au excluded, see below), the antibonding-population-near-frontier metric (normalized) reaches ρ=−0.155, p=1.8×10⁻⁴ (n=579) against `formation_energy_per_atom` — markedly weaker than the ρ=−0.282 (n=346) reported at the previous scale, though still real and significant.** The sign is unchanged: more occupied antibonding character near the frontier still associates with *more negative* (more stable) formation energy, the opposite of the naive Peierls/Jahn-Teller reading. **The `bond_type=covalent` subgroup — previously the single most robust stratified result in the entire project (ρ=−0.551, p=0.0009 at n=33) — does NOT survive at the current scale** (ρ=−0.046, p=0.74, n=54): the earlier result does not replicate once the covalent sample nearly doubled, the same small-n fragility this project's own methodology has now caught four times (mincut/`bond_type`, this metric's own old `bond_type=ionic` result, min-cut/anion, and now this). **New at this scale: `bond_type=mixed` (Zintl-type) is the strongest surviving stratum** (ρ=−0.538, p=6.6×10⁻⁶, n=62) — mirroring the same pattern found for reaction-ICOHP (`REPORT_reaction_icohp.md`) at the same scale. `is_metal=False` still survives (ρ=−0.333, p=6.7×10⁻⁷, n=212), and `is_metal=True` has converged onto an identical population to `bond_type=metallic` (same n, ρ, p — same convergence pattern documented for reaction-ICOHP).

## 0. Data-quality note: Cu4Au excluded

One compound, `extension_Cu4Au_mp-1225761`, was excluded from this and every other analysis in the project after a systematic quality audit (`analysis/audit_lobster_quality.py`) found its LOBSTER run catastrophically corrupted: 100% of k-points failed orthonormalization and `bandOverlaps.lobster` shows a maximum deviation of 11.4 (vs. this project's 0.1 red-flag threshold). Confirmed a unique case in the whole 597-structure dataset, not a wider problem. See `mp_metadata.json`'s `quality_excluded_reason` field and the activity report's Limits section for detail.

## 1. The zero floor

`antibond_w_raw` is exactly (< 10⁻⁶) zero for **173/588 compounds (29.4%)**:

| Group | n | n zero | frac zero |
|---|---:|---:|---:|
| metal (`is_metal=True`) | 372 | 126 | **33.9%** |
| gapped (`is_metal=False`) | 216 | 47 | 21.8% |
| bond_type=metallic | 372 | 126 | **33.9%** |
| bond_type=covalent | 56 | 11 | 19.6% |
| bond_type=ionic | 98 | 21 | 21.4% |
| bond_type=mixed | 62 | 15 | 24.2% |

Same qualitative reading as at every prior scale: metals are far more likely to show zero antibonding population in this window than gapped compounds, consistent with the metric's own construction (E_ref=E_F for metals sits inside a partially-filled band; E_ref=VBM for gapped compounds is a band-edge state more likely to carry bonding/antibonding mixing).

## 2. Correlations vs. `formation_energy_per_atom`

### 2.1 Primary metric (ΔE=1.0)

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 579 | **−0.1550** | **1.8×10⁻⁴** |
| bond_type=covalent | 54 | −0.0464 | 0.739 |
| bond_type=ionic | 96 | −0.1958 | 0.056 |
| bond_type=metallic | 367 | 0.0563 | 0.282 |
| bond_type=mixed | 62 | **−0.5376** | **6.6×10⁻⁶** |
| is_metal=False | 212 | **−0.3334** | **6.7×10⁻⁷** |
| is_metal=True | 367 | 0.0563 | 0.282 |

Sensitivity (raw, `all` group): n=581, ρ=−0.1266, p=2.2×10⁻³ — same sign, smaller magnitude, so most of the normalized signal survives un-normalizing.

### 2.2 ΔE sensitivity (raw, `all` group)

| ΔE (eV) | n | ρ | p |
|---:|---:|---:|---:|
| 0.5 | 581 | −0.1300 | 1.7×10⁻³ |
| 1.0 (primary) | 581 | −0.1266 | 2.2×10⁻³ |
| 2.0 | 581 | −0.1175 | 4.6×10⁻³ |

Same conclusion as at every prior scale: the global correlation is essentially insensitive to window width.

## 3. Is the global signal a between-group effect?

Two-sample Mann-Whitney U tests, metal vs. gapped:

| Variable | metal median | gapped median | p |
|---|---:|---:|---:|
| `formation_energy_per_atom` | −0.0435 (n=367) | −0.5831 (n=212) | **1.9×10⁻²⁴** |
| `antibond_w_normalized` | 0.0069 | 0.0224 | **6.4×10⁻⁵** |

Kruskal-Wallis across `bond_type` (covalent/ionic/metallic/mixed — a "mixed"/Zintl category has been added since the previous version of this report):

| Variable | H | p |
|---|---:|---:|
| `formation_energy_per_atom` | 155.1 | 2.1×10⁻³³ |
| `antibond_w_normalized` | 32.9 | 3.3×10⁻⁷ |

Both group-level tests remain overwhelming — pooling across `bond_type` is still not a valid way to read the `all` row in §2.1.

**What has genuinely changed since the last scale**: previously, `bond_type=covalent` was the standout stratum (the only one to survive, and the strongest result of any kind in the project) while `is_metal` split into two survivors. At the current scale, **`bond_type=covalent` no longer survives at all**, `bond_type=ionic` is borderline (p=0.056), `bond_type=metallic` stays non-significant with a *flipped* sign (now positive, was negative in the pooled/covalent readings), and the stratum that now carries the strongest signal is the newly-populated `bond_type=mixed` (Zintl) category — a category that essentially did not exist in the dataset at the scale the covalent result was originally reported. `is_metal=True` has converged onto an identical population to `bond_type=metallic` (same n/ρ/p in every row), the same convergence already documented for reaction-ICOHP at this scale (`REPORT_reaction_icohp.md`) — it is `bond_type`, not `is_metal`, that carries whatever stratified signal remains.

## 4. Vs. `energy_above_hull`

| Group | n | ρ | p |
|---|---:|---:|---:|
| all | 585 | 0.0818 | 0.048 |
| bond_type=covalent | 56 | −0.2505 | 0.063 |
| bond_type=ionic | 98 | −0.0169 | 0.869 |
| bond_type=metallic | 369 | **0.2323** | **6.5×10⁻⁶** |
| bond_type=mixed | 62 | −0.1412 | 0.274 |

Notable: against this target, it is `bond_type=metallic` that survives (positive sign — more antibonding character near the frontier associates with *higher* hull distance, the opposite direction from the `formation_energy_per_atom` finding in that same stratum). The pooled `all` row is only borderline (p=0.048) and should not be over-read given the stratified picture is this inconsistent target-to-target.

## 5. On the sign

Unchanged from every prior scale — the Peierls/Jahn-Teller framing in `METRIC_DEFINITION_antibonding.md` §1 predicts the opposite sign from what's observed against `formation_energy_per_atom`, and the same two non-exclusive readings apply:

- **Confound reading**: gapped/covalent compounds have both more negative formation energies and more antibonding population near a band edge for structural reasons unrelated to instability.
- **Target mismatch reading**: `formation_energy_per_atom` measures formation from elements, not proximity to a symmetry-lowering distortion of the already-formed structure — §4 above (a genuinely different, sign-inconsistent picture against `energy_above_hull`) is itself evidence for this reading: if the metric tracked a single real physical effect, both targets should agree, and they don't.

## 6. Limits and next steps

- **The `bond_type=covalent` result is retracted.** It does not survive at n=54 (was reported at n=33). Treat it as a small-sample artifact of the earlier scale, not a real finding — the fourth time this project's own methodology has caught exactly this failure mode.
- **`bond_type=mixed` (Zintl) is now the strongest surviving stratum** (ρ=−0.538, p=6.6×10⁻⁶) — new, not previously testable (the category did not exist at the prior scale). Worth prioritizing as the next dedicated follow-up, echoing the same finding already flagged for reaction-ICOHP.
- **`is_metal=False` remains the most consistently-surviving stratum across scales** (p=6.7×10⁻⁷ now, was p=0.001) — the one result in this report that has strengthened, not weakened, with more data.
- **Sign interpretation is still unresolved** (§5) — now further complicated by the `energy_above_hull` result (§4) actively disagreeing in direction with the `formation_energy_per_atom` result within `bond_type=metallic`. This report establishes a set of statistical associations, not a validated physical mechanism.
- **Still no SISSO.** The signal is bond-type-conditional, inconsistent in sign across targets, and its previously-headline stratum has since retracted — a pooled feature search would be premature here.

---

*Code*: `analysis/compute_antibonding_all_full.py` (full-scale, generic re-implementation superseding the old `compute_antibonding_all.py`/`compute_antibonding_extension.py`/`compute_icohp_antibonding_maxhull.py` batch-specific scripts), `analysis/audit_lobster_quality.py` (Cu4Au exclusion). `cohp_extraction.py` untouched since the pilot-validated version. *Data*: `analysis/icohp_antibonding_full.csv` (588 rows), `analysis/icohp_icobi_bondtype.csv` (bond-type/is_metal source). *Figures*: not regenerated this pass — `analysis/figures_antibonding/` still reflects the previous scale.
