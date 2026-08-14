# Primitive vs. conventional cell — pilot verdict (mission #2)

**Verdict: the primitive/conventional cell bias explains essentially NONE of the near-zero correlation from the 186-compound campaign — and it does so for a different reason than hypothesized.** Switching to the conventional cell changes the percolation weight substantially for 4/6 pilot compounds (NaCl -55x, LiBr -20x, Si -3x, C-rhomb -3x — all well beyond the ~40% run-to-run DFT noise floor established by the two negative-control compounds), but in every single case the weight gets **smaller**, not larger/closer to the strongest bond as the pilot report's NaCl hypothesis predicted. Direct path inspection confirms the strongest Na–Cl bond (−0.593 eV) still does **not** participate in NaCl's minimum-weight cycle in the conventional cell either — the algorithm instead finds an even weaker long-range path, because the conventional cell simply has 4x more atoms and therefore more opportunities for a cheap indirect route. **Do not extend to the 186-compound scale under this hypothesis** (see recommendation at the end).

## Method

For each of the 6 pilot compounds, `SpacegroupAnalyzer(CONTCAR).get_conventional_standard_structure()` from the *relaxed* primitive-cell structure, then a fresh static VASP+LOBSTER run and `percolation_path.py` (unmodified) on the conventional cell.

**Methodology confound found and fixed before drawing conclusions**: the original 6 primitive-cell results (used throughout the 186-compound campaign) were computed with the pilot campaign's old fixed `NBANDS=100`, while `prepare_vasp_lobster.py` was since updated (180-compound campaign) to derive `NBANDS` dynamically from the LOBSTER basis size. Comparing conventional-cell (new NBANDS) against original primitive-cell (old NBANDS) would have confounded the cell-choice question with a DFT-methodology change — confirmed by re-running AlNi's *primitive* cell (identical structure, byte-for-byte same lattice vectors and positions) and finding a 7x difference purely from `NBANDS=100` vs `NBANDS=13`. **All numbers below use a matched primitive-v2 (dynamic NBANDS) vs. conventional (dynamic NBANDS) comparison**, both re-run with current code.

## Results

| Compound | Expansion | SG (conventional) | Weight primitive-v2 (eV) | Weight conventional (eV) | Ratio conv/prim |
|---|---:|---|---:|---:|---:|
| NaCl (hull, ionic) | ×4.0 | Fm-3m | 0.03088 | 0.00055 | **0.018** |
| Si (hull, covalent) | ×4.0 | Fd-3m | 0.00318 | 0.00101 | **0.318** |
| AlNi (hull, metallic) | ×1.0 | Pm-3m | 0.02069 | 0.02932 | 1.417 *(negative control)* |
| LiBr (metastable, ionic) | ×4.0 | Fm-3m | 0.09880 | 0.00490 | **0.050** |
| C rhomb. (metastable, covalent) | ×2.0 | C2/m | 0.00027 | 0.00009 | **0.333** |
| BeCu (metastable, metallic) | ×1.0 | Pm-3m | 0.00054 | 0.00034 | 0.630 *(near noise floor)* |

AlNi and BeCu have primitive = conventional cell (already-simple-cubic Pm-3m); re-running them anyway gives the **DFT noise floor** for "no structural change at all": AlNi drifts 42%, BeCu 37%. Every compound whose cell actually changed drifted 3x–55x — unambiguously above that floor.

### Direct path check for NaCl (the core diagnostic)

```
Primitive:     weight=0.03088, direction=b
  path = [{'from_atom': 1, 'to_atom': 1, 'translation': [0, 1, 0], 'weight': 0.03088, 'bond_label': '39'}]
  (2nd-neighbor Cl-Cl, single hop)

Conventional:  weight=0.00055, direction=a
  path = [{'from_atom': 0, 'to_atom': 0, 'translation': [1, 0, 0], 'weight': 0.00055, 'bond_label': '3'}]
  (a much longer-range Na-Na or Cl-Cl shell, single hop)
```

`icohp_min` (strongest bond, Na–Cl nearest-neighbor) = **−0.593 eV** in both cells. Neither the primitive nor the conventional minimum-weight path uses it. The conventional cell's winning path is a *single* long-range bond of weight 0.00055 eV — over 1000x weaker than the strongest bond, and weaker even than the primitive cell's already-weak 2nd-shell path. Going conventional did not surface the strong bond; it surfaced a *weaker* one, because 8 atoms/cell instead of 2 gives the minimum-weight search strictly more candidate routes to be cheap with.

## Why this refutes the pilot report's mechanism (not just the numbers)

The pilot report's NaCl discussion hypothesized: *the strongest bond sits entirely inside the primitive cell (translation (0,0,0)) and therefore cannot contribute to a non-contractile cycle; a conventional cell would let it participate.* This is only half right: the strongest Na–Cl bond genuinely doesn't cross the primitive cell boundary. But the conventional cell doesn't fix that — the reference cell is a different shape, but the search still finds whatever is globally cheapest, and a bigger cell (more atoms, more bond-length shells represented) mechanically gives the minimum-weight search **more long-range near-zero ICOHP tail values** to exploit. The problem isn't "primitive vs. conventional" per se; it's that **the minimum-weight non-contractile cycle is structurally biased toward whatever the weakest bond in the entire ICOHPLIST cutoff sphere happens to be**, and larger cells have more such bonds to choose from. This is a distinct, more general problem than the one originally diagnosed.

## Recommendation

**Do not proceed to the conditional 186-compound (or even 88-compound metallic-only) conventional-cell extension.** The observed change, while real and substantial, does not support the hypothesized fix — re-running everything in conventional cells would very plausibly make the metric *more* noise-dominated (more atoms → more weak long-range escape routes), not less, likely leaving the correlation with `energy_above_hull` just as close to zero or worse. That compute budget is better spent elsewhere.

This result is exactly why mission #3 (started the same night, in parallel) pivots to **network dimensionality** (a relative-strength-thresholded connectivity question, immune to this failure mode by construction — it asks "does *any* sufficiently strong path exist", not "what is the single cheapest one") and **periodic min-cut** (which requires cutting *every* parallel path, not just finding one cheap one, and is validated in `tests/test_periodic_mincut.py` specifically against a synthetic case designed to catch this exact "collapses to the single cheapest path" failure mode). Both are better-motivated responses to this finding than a conventional-cell campaign would have been.

---

*Data*: `mp_dataset/results_primitive_pilot_v2.csv`, `mp_dataset/results_conventional_pilot.csv`, `mp_dataset/structures_primitive_pilot_v2/`, `mp_dataset/structures_conventional_pilot/`. *Code*: `mp_dataset/generate_conventional_pilot.py`.
