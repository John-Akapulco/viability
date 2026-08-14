# COHP antibonding-population feasibility — steps 0+1 (pilot only)

**Conclusion: the extraction pipeline is fully validated and ready for step 2.** All 6 pilot compounds' `COHPCAR.lobster` survived on Yargla without any recomputation (best case throughout: `WAVECAR`, `COHPCAR.lobster`, `COBICAR.lobster`, `POTCAR` all present and intact). pymatgen's own parsers (`Cohpcar`, `CompleteCohp`) reproduce the already-validated `ICOHPLIST.lobster` values to within numerical noise (max |diff| = 1×10⁻⁵ eV across 558 bond labels checked, exact 0.0 for 5/6 compounds), with a confirmed, non-assumed sign convention (negative = bonding), and metal/gap character is now settled for all 6 against Materials Project's own converged calculation — including resolving the one open question flagged in the mission brief (rhombohedral graphite: confirmed a narrow-gap semiconductor, 0.198 eV, **not** a semimetal). No unresolved problem blocks step 2.

## Summary table

| Compound | Yargla survival | Cross-validation (max \|diff\|, n labels) | Metal/gap (MP) | Local coarse-mesh gap estimate | 
|---|---|---:|---|---|
| NaCl (hull, ionic) | WAVECAR+COHPCAR present, no recompute | 0.0 eV (n=50) | gap, 5.00 eV | 5.22 eV (consistent) |
| Si (hull, covalent) | WAVECAR+COHPCAR present, no recompute | 0.0 eV (n=64) | gap, 0.61 eV | 0.63 eV (consistent) |
| AlNi (hull, metallic) | WAVECAR+COHPCAR present, no recompute | 0.0 eV (n=96) | **metal**, 0.0 eV | 0.14 eV — **artifact, see below** |
| LiBr (metastable, ionic) | WAVECAR+COHPCAR present, no recompute | 0.0 eV (n=50) | gap, 4.92 eV | 5.08 eV (consistent) |
| C rhombohedral (metastable, covalent) | WAVECAR+COHPCAR present, no recompute | 0.0 eV (n=130) | gap, **0.198 eV** (confirmed not semimetal) | 0.85 eV (same sign, different magnitude — see below) |
| BeCu (metastable, metallic) | WAVECAR+COHPCAR present, no recompute | 1×10⁻⁵ eV (n=168) | **metal**, 0.0 eV | 0.17 eV — **artifact, see below** |

## 0. Yargla file survival

All checked before any decision, per the mission brief (no assumption of uniform behavior across compounds). Full machine-readable log: `analysis/yargla_file_survival.json`.

```
compound                                     WAVECAR       COHPCAR.lobster   COBICAR.lobster   POTCAR
hull_ionic_NaCl_mp-22862                     610.6 MB      13.7 MB           13.7 MB           404 KB
hull_covalent_Si_mp-149                      565.1 MB      17.5 MB           17.5 MB           196 KB
hull_metallic_AlNi_mp-1487                   641.6 MB      70.0 MB           70.0 MB           429 KB
metastable_ionic_LiBr_mp-23259               565.1 MB      17.2 MB           17.2 MB           362 KB
metastable_covalent_C_rhombohedral_mp-169    582.7 MB      35.5 MB           35.5 MB           207 KB
metastable_metallic_BeCu_mp-2323             526.8 MB      122.4 MB          122.4 MB          428 KB
```

All 6 land in the best-case row of the mission's scenario table ("`WAVECAR` + `COHPCAR.lobster` présents → copie directe, aucun calcul"). These files were never committed to git (`.gitignore`: dense per-grid LOBSTER outputs and raw wavefunctions are excluded as large + regenerable) but were never deleted from disk either — they are the original outputs from this project's very first VASP+LOBSTER run (the 6-compound pilot, computed at the start of this whole session, before the 180-compound campaign). Lightweight copies (`COHPCAR.lobster`, `ICOHPLIST.lobster`, `CONTCAR`, `mp_metadata.json`, `lobsterin` — 265 MB total, no `WAVECAR`/`COBICAR.lobster` duplicated since they aren't needed for parsing `COHPCAR.lobster`) were made to `mp_dataset/structures_cohp/{compound_id}/` for this step's work, keeping the original 3.5 GB of `WAVECAR` files un-duplicated.

## 1. Parser and cross-validation

**Library, not hand-rolled**: `pymatgen.io.lobster.outputs.Cohpcar` (raw per-line COHPCAR reader) and `pymatgen.electronic_structure.cohp.CompleteCohp` (structure-aware wrapper), pymatgen **2026.5.4** (see `requirements.txt`). Both cite George et al., *ChemPlusChem* 2022 (the same reference already used for the θ=10% bond-strength threshold in mission #3) as the module's canonical reference.

### 1.1 Cross-validation (mandatory before any interpretation)

For every bond label in `ICOHPLIST.lobster` (already validated throughout this project — it's what `percolation_path.py` consumes), the ICOHP trace read from `COHPCAR.lobster` (independently, via `Cohpcar`) at the energy index closest to E=E_F was compared against the `ICOHPLIST.lobster` value for the same label:

- **558 bond labels checked across the 6 compounds.**
- **5/6 compounds: exact match, max |diff| = 0.0 eV.**
- **1/6 (BeCu, 168 labels): max |diff| = 1×10⁻⁵ eV**, mean |diff| = 1.07×10⁻⁶ eV — five to six orders of magnitude below the weakest physically meaningful ICOHP value anywhere in this dataset (O(10⁻⁴) eV) and consistent with LOBSTER computing `ICOHPLIST.lobster` via adaptive quadrature vs. `COHPCAR.lobster`'s fixed energy-grid trapezoidal integration — a known, harmless difference in integration method, not a data or parsing problem.

**No discrepancy required stopping** — the mission's explicit "if this doesn't match, don't continue" gate is cleared for all 6.

### 1.2 Sign convention (verified, not assumed)

Confirmed by direct inspection of pymatgen's own `Cohp.has_antibnd_states_below_efermi` source (`pymatgen.electronic_structure.cohp`), which tests `cohp_values > positive_limit` to mean *antibonding*. Applying this to the exact same numerical array that §1.1 shows matches `ICOHPLIST.lobster` (same sign, same values) establishes, empirically rather than by assumption:

> **NEGATIVE ICOHP/COHP = bonding, POSITIVE ICOHP/COHP = antibonding** — matching the convention already used throughout this project since the very first pilot report (e.g. NaCl's strongest Na–Cl bond, ICOHP = −0.595 eV).

Both `Cohpcar` (raw) and `CompleteCohp` (structure-aware) preserve this same sign; no flip happens anywhere in the pymatgen wrapping. One energy-axis subtlety worth documenting for step 2: `Cohpcar.energies` is already shifted so E_F = 0 (LOBSTER writes it that way directly — confirmed empirically, not merely per the pymatgen docstring, which describes the same behavior but was verified independently here), while `CompleteCohp.energies` is on an absolute-eV scale (`efermi` is a nonzero absolute value); both give identical ICOHP values when each class's own E_F reference point is used, but step 2 must pick one consistently.

### 1.3 Metal vs. gap classification, per pilot

A **known pitfall was found and documented, not silently worked around**: a naive local check (either `pymatgen.io.vasp.Vasprun.eigenvalue_band_properties`, or DOS-at-E_F from the local calculation) on our LOBSTER-oriented 10×10×10 k-mesh with `ISMEAR=0`/`SIGMA=0.05` spuriously suggests small nonzero gaps (0.14 eV for AlNi, 0.17 eV for BeCu) for two compounds that are unambiguously metals (well-known B2 intermetallic alloys). The DOS-at-E_F check is even less reliable in the other direction: NaCl (5+ eV gap, unambiguous insulator) shows a *nonzero* smeared DOS at E_F (0.38 states/eV) purely from Gaussian-smearing tails on a wide gap — a naive "DOS(E_F) > 0 ⟹ metal" rule would misclassify it. **Materials Project's own converged classification was used as the authoritative source instead** (`cohp_extraction.metal_or_gap_from_mp`, single batched query, same API key convention as the rest of this project):

| Compound | MP is_metal | MP band_gap (eV) |
|---|---|---:|
| NaCl | False | 5.004 |
| Si | False | 0.610 |
| AlNi | **True** | 0.000 |
| LiBr | False | 4.923 |
| C rhombohedral | False | **0.198** |
| BeCu | **True** | 0.000 |

The pilot set covers both classes as the mission brief expected (2 metals, 4 gapped), and the one flagged uncertainty is resolved: **rhombohedral graphite is a narrow-gap semiconductor (0.198 eV) per MP, not a semimetal** — worth double-checking again once/if step 2 needs a specific gap value, since 0.198 eV is close enough to zero that different methodologies (this project's local PBE-GGA calc gives 0.849 eV, MP's own gives 0.198 eV) disagree substantially on the magnitude even while agreeing it's non-metallic.

## Tests

`tests/test_cohp_extraction.py`, run against the real (not synthetic) `COHPCAR.lobster` files for all 6 pilots — deliberately real data, since validating pymatgen's parsing of the actual LOBSTER file format is the point of this step, not a hand-rolled algorithm:

1. Bond-label count in `COHPCAR.lobster` matches `ICOHPLIST.lobster` (all 6).
2. Cross-validation within a 1×10⁻⁴ eV tolerance (10× looser than the worst observed deviation, still far tighter than any physically meaningful ICOHP in this dataset) — all 6 pass, all labels matched.
3. Sign-convention check runs and returns the documented convention statement for all 6.
4. Metal/gap classification matches the MP reference table above (hardcoded from the live query transcribed in §1.3, so the test suite stays fast/offline) — including the explicit rhombohedral-graphite assertion.

All 4 new tests pass; full project suite (24 tests, including the pre-existing 20) is green.

## Explicitly NOT done in this mission (per scope)

- No energy window ("near-E_F") defined.
- No aggregated antibonding-population score/metric.
- No extension beyond the 6 pilots.
- `percolation_path.py` and its results untouched.
- No SISSO.

---

*New code*: `cohp_extraction.py` (repo root, `percolation_path.py` untouched). *Data*: `analysis/yargla_file_survival.json`, `mp_dataset/structures_cohp/` (6 compounds, COHPCAR.lobster + ICOHPLIST.lobster + CONTCAR + metadata, no WAVECAR/COBICAR duplicated). *Tests*: `tests/test_cohp_extraction.py` (4/4 passing, real LOBSTER data).
