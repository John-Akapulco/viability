# reactions_dataset/

Expected format for `reaction_analysis.cli` (see `reaction_analysis/schema.py`
for the exact field definitions):

```
reactions_dataset/
  entries/
    <compound_id>/
      entry.json          # a CompoundEntry, serialized (e.g. entry.model_dump_json())
      ICOHPLIST.lobster    # kept for provenance / re-parsing if needed
      ICOBILIST.lobster    # optional
      POSCAR (or CONTCAR/.cif)
  reactions/
    <reaction_id>.json    # a Reaction, serialized
```

`entries/<id>/entry.json` is expected to already exist by the time the
CLI runs -- produced ahead of time via `reaction_analysis.parse_lobster.parse_compound_entry()`,
not derived on the fly by the CLI itself.

**Populated at full case-1 scale.** First validated on a 6-compound batch
(`analysis/validate_reaction_analysis_case1.py`: Ca, O2, N2, As, Pd, plus
targets CaO / Ca3N2 / AsPd2, 3 `decomposition_to_elements` reactions) --
cross-checked delta_per_atom_eV against the existing ad hoc
`reaction_icohp.py`'s (mission #5) `delta_icohp_per_atom` for the same
reactions to floating-point precision (opposite sign convention, see
`delta.py` docstring). Then extended to every case-1-eligible compound in
`mp_dataset/structures/` (`analysis/populate_reaction_analysis_case1_full.py`,
same scope as `analysis/compute_reaction_icohp_case1.py`): 192 targets,
254 distinct `entries/` (targets + shared elemental references), all
192/192 cross-checked against `reaction_icohp_case1.csv` and matching.
Raw `ICOHPLIST.lobster`/`CONTCAR` are NOT duplicated into `entries/<id>/`
at this scale (unlike the 6-compound validation batch) -- they already
live in `mp_dataset/structures/<id>/`, which each entry's
`CompoundEntry.source_path` points back to.

See `tests/test_balance.py`/`tests/test_delta.py` for the synthetic
hand-verified fixtures (Na/K/NaK/Na2K toy compounds) and
`reaction_analysis/cli.py`'s docstring for how the CLI runs over a
dataset in this format.
