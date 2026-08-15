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

**Empty as of this commit.** No real LOBSTER production data has been run
through `reaction_analysis` yet -- this directory only documents the
expected layout. See `tests/test_balance.py`/`tests/test_delta.py` for a
fully worked, hand-verified example (Na/K/NaK/Na2K toy compounds) and
`reaction_analysis/cli.py`'s docstring for how a real dataset would be
run once populated.
