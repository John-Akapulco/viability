# Reaction-ICOHP: definition and validation

## 1. Physical motivation

Every descriptor built so far in this project (percolation weight,
network dimensionality, min-cut, antibonding population) looks at a
compound's *own* bonding topology in isolation. None of them ask the more
direct thermochemical question: is a compound's total bond population
"worth more" in ICOHP terms than the same atoms would have in some
reference configuration (its constituent elements, a polymorph, another
compound)? That's exactly the question `formation_energy_per_atom`
answers for total DFT energy -- this descriptor is the ICOHP analog.

Three reaction types, user-requested in this order:
1. **decomposition into elements**: AaBb -> a A + b B (standard-state references)
2. **polymorphs**: same composition, different structure, direct comparison
3. **decomposition into another compound + elements**: general balanced reaction

## 2. Core intensive quantity: `icohp_per_atom`

Sum every bond in `ICOHPLIST.lobster` (LOBSTER's own accounting of every
symmetry-inequivalent periodic bond in the cell, unfiltered -- the same
raw source `percolation_path.py`'s `icohp_sum` column already uses for the
whole-cell case, `bond_pair=None`) and divide by the number of atoms in
the cell (from `CONTCAR`, the relaxed structure, same file
`cohp_extraction.py`'s `CompleteCohp` loader uses).

This is assumed to be a bulk/converged, transferable quantity -- exactly
the same assumption `formation_energy_per_atom` makes about total DFT
energy per atom. Sign convention (established project-wide, see
`cohp_extraction.py`'s `sign_convention_check`): more negative ICOHP =
more bonding.

## 3. Case 1 and case 3 share one primitive: `reaction_delta_icohp`

Given a reactant directory and a list of product directories,
`reaction_icohp.py` balances the reaction by composition alone
(`pymatgen.analysis.reaction_calculator.Reaction`, on reduced formulas),
scales the balanced reaction up so the reactant side matches the
reactant's *actual* cell (`n_sites` atoms, not the reduced formula), then:

```
delta_icohp_total = icohp_total(reactant's own cell)
                     - sum_p [ coeff_formula_units(p) * atoms_per_fu(p) * icohp_per_atom(p) ]
delta_icohp_per_atom = delta_icohp_total / n_sites(reactant)
```

Case 1 is not special-cased -- it is simply the case where every supplied
product happens to be a pure element. This is why case 3 (Mn2O7 -> MnO2 +
O2) uses the exact same function as case 1 (Mn2O7 -> Mn + O2).

**Mass-balance self-check** (built into validation, not the library code):
for every tested reaction, `sum(coeff_formula_units * atoms_per_fu)` over
all products must equal the reactant's `n_sites` exactly (atoms in = atoms
out). Confirmed to machine precision for Ca3N2 -> Ca + N2 (24 Ca atoms +
16 N atoms = 40, matching Ca3N2's own 40-site cell) and Mn2O7 -> MnO2 + O2.

## 4. Case 2: polymorphs need no reaction balancing

`compare_polymorphs()` takes a list of same-composition compound
directories and directly compares `icohp_per_atom` -- no stoichiometry
involved, since the atom count is already 1:1 comparable.

## 5. Validation on already-computed extension compounds

Three worked examples (real LOBSTER output, not synthetic):

**Ca3N2 -> 3 Ca + N2** (case 1): `delta_icohp_per_atom = +4.44`. Positive
-- Ca3N2 has *less* net ICOHP-bonding per atom than its elemental
references, dominated by N2's triple bond (icohp_per_atom = -16.4, one of
the strongest bonds physically possible) heavily outweighing Ca3N2's
mixed ionic/covalent bonding (icohp_total/40 sites = -2.3). This is
expected and important to flag: **ICOHP measures orbital-overlap bond
population, not electrostatic/Madelung lattice energy** -- a real ionic
nitride's true (very negative) formation energy is dominated by
electrostatics that ICOHP does not see at all. This metric should NOT be
expected to reproduce `formation_energy_per_atom` for strongly ionic
compounds; the more probable target here is `bond_type=covalent` or
`bond_type=metallic` subgroups, mirroring the one persistent finding of
the whole project so far (the antibonding metric's only surviving
subgroup was `bond_type=covalent`, [[project-viability-antibonding-status]]-linked
finding).

**Mn2O7 -> MnO2 + O2** (case 3): `delta_icohp_per_atom = +1.12`. Same
sign/direction as Ca3N2, consistent story (O2's double bond is a strong
per-atom bonding reference, same caveat applies).

**Carbon allotrope polymorphs** (case 2, 5-way: graphite, diamond,
lonsdaleite, M-carbon, W-carbon): lonsdaleite and diamond (sp3,
4-coordinate) have the most negative `icohp_per_atom` (most bonding);
graphite (sp2, weak interlayer coupling) has the least, by a wide margin
(+1.83 vs. lonsdaleite). This is a real, informative **anti-correlation**
with `energy_above_hull` -- graphite is the most thermodynamically stable
of the three real phases (EAH=0.003) yet the *least* bonding by this
metric, because ICOHP only sees covalent orbital overlap and does not
capture the weak interlayer van der Waals attraction that (along with
entropy) makes graphite the true ground state. Documented here as a
concrete demonstration that `icohp_per_atom` is not a drop-in stability
predictor -- exactly the target-mismatch caveat already flagged for the
antibonding metric (`REPORT_antibonding.md`, "Target mismatch reading").

**TiO2 high-pressure polymorphs** (case 2, 3-way): baddeleyite-type has
the most bonding, but this does NOT match the EAH ranking either
(TiO2-II has the lowest EAH but is not the most-bonding structure) --
same lesson, different compound family.

## 6. What this does NOT establish yet

- No statistical test against `formation_energy_per_atom` or
  `energy_above_hull` across a real sample size has been run -- these are
  n=1 to n=5 worked examples validating the *arithmetic*, not a
  correlation claim. That is the explicitly separate next mission once
  case 1 is extended to compounds with all-element references available
  (in progress: `download_elements_reference.py`, 50 elemental references
  submitted to SLURM 2026-08-14, `mp_dataset/structures/extension_Ti_mp-46`
  submitted separately for the TiO2 polymorphs).
- The strong ionic-compound sign/magnitude effect seen in Ca3N2 and Mn2O7
  strongly suggests this metric, like every other one in this project,
  needs `bond_type` stratification before any global correlation is
  trusted -- expect it to behave very differently for ionic vs. covalent
  vs. metallic compounds.
- Elemental reference structure selection for the 50-element bulk batch
  (`download_elements_reference.py`) was RULE-BASED (lowest-energy
  non-theoretical Materials Project entry), not individually
  literature-verified like the original hand-picked Ca/Mn/S/Ti references
  -- 7 elements got a manual override where the automated pick
  contradicted a well-known real standard state by a DFT-noise-level
  margin (Ag, In, Rb, Cs, Se, Sn, Ta); the rest carry an explicit
  "near-degenerate, not individually verified" flag in their
  `mp_metadata.json` note where applicable. Treat any single element's
  reference structure as provisional until spot-checked, same spirit as
  the CaN/CaO LOBSTER-quality caveat elsewhere in this project.
- No LOBSTER band-overlap `maxDeviation` quality check has been run yet
  on any of the 50 new element calculations (the same check that flagged
  CaN/CaO) -- do before trusting any individual reaction result at scale.
