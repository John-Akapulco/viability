# reaction_analysis test fixtures

Only two new toy compounds live here: `compound_K` (pure K, single atom)
and `compound_Na2K` (2 Na + 1 K), both in the exact minimal style of
`examples/dataset/` (round ICOHP numbers, tiny cubic cells, hand-checkable
sums). They exist because `examples/dataset/` alone (`compound_A`=Na,
`compound_B`=NaK, `compound_C_disconnected`=Na) cannot cover all three
reaction types the schema needs to test:

- **decomposition_to_elements**: `Na2K -> 2 Na + K` needs a pure-K
  reference, which didn't exist yet (`compound_K`, added here).
- **decomposition_to_compound_and_elements**: `Na2K -> NaK + Na` reuses
  `compound_B` (NaK) and `compound_A` (Na) directly.
- **polymorph_transition**: `compound_A` vs `compound_C_disconnected`
  (both pure Na, different bond topology) already exist and are reused
  as-is, no new fixture needed.

`tests/test_parse_lobster.py`, `tests/test_balance.py`, and
`tests/test_delta.py` reference `examples/dataset/` for the reused
compounds and this directory only for `compound_K`/`compound_Na2K`.
