"""Merge select_marginal_formation_energy.py's periodic-table-wide pool
(50 entries, mostly metallic near-zero-FE alloys) and
select_marginal_ionic.py's chemistry-targeted pool (30 entries, alkali/
alkaline-earth nitrides and polyphosphides, mostly non-metallic) into
one candidate list, re-applying chemsys diversity across the COMBINED
pool (not just within each source separately) -- some chemsys (e.g.
N-Na with 3 NaN3 polymorphs, P-Li with 2 LiP5 polymorphs) exceed the cap
once merged even though neither source list did on its own.

Selection-only, no download. Writes
mp_dataset/marginal_candidates_merged.json.
"""

import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
MAX_PER_CHEMSYS = 2

SOURCES = [
    (HERE / "marginal_formation_energy_candidates.json", "periodic_wide"),
    (HERE / "marginal_ionic_candidates.json", "ionic_targeted"),
]


def main() -> None:
    merged = []
    for path, source in SOURCES:
        records = json.loads(path.read_text())
        for r in records:
            r["source"] = source
            merged.append(r)

    n_by_source = defaultdict(int)
    for r in merged:
        n_by_source[r["source"]] += 1
    print(f"By source: {dict(n_by_source)}")

    # Duplicate mp_id safety check (should not happen, sources are disjoint
    # element universes, but don't silently merge if it ever does).
    seen_ids = set()
    deduped = []
    for r in merged:
        if r["mp_id"] in seen_ids:
            print(f"  DROPPED duplicate mp_id {r['mp_id']} ({r['formula']}, already kept from another source)")
            continue
        seen_ids.add(r["mp_id"])
        deduped.append(r)

    # Chemsys diversity across the COMBINED pool: prefer experimental
    # over theoretical, then smallest |FE| (closest to the marginal
    # boundary either source was built around).
    deduped.sort(key=lambda r: (bool(r["theoretical"]), abs(r["formation_energy_per_atom_eV"])))
    by_chemsys = defaultdict(list)
    kept = []
    dropped_for_diversity = []
    for r in deduped:
        if len(by_chemsys[r["chemsys"]]) >= MAX_PER_CHEMSYS:
            dropped_for_diversity.append(r)
            continue
        by_chemsys[r["chemsys"]].append(r)
        kept.append(r)

    kept.sort(key=lambda r: (r["source"], r["chemsys"], r["formula"]))

    out_path = HERE / "marginal_candidates_merged.json"
    out_path.write_text(json.dumps(kept, indent=2))

    print(f"\n{len(deduped)} unique entries -> {len(kept)} after chemsys diversity "
          f"(cap {MAX_PER_CHEMSYS}/chemsys), {len(dropped_for_diversity)} dropped for diversity")
    n_by_source_final = defaultdict(int)
    n_metal = n_nonmetal = 0
    for r in kept:
        n_by_source_final[r["source"]] += 1
        if r["is_metal"]:
            n_metal += 1
        else:
            n_nonmetal += 1
    print(f"Final by source: {dict(n_by_source_final)}")
    print(f"Final metallic={n_metal}, non-metallic={n_nonmetal}")
    print(f"Wrote {len(kept)} entries to {out_path}\n")

    print(f"{'source':<16}{'formula':<12}{'mp_id':<12}{'chemsys':<10}{'FE(eV/at)':>11}  "
          f"sg{'':<10}nsites  is_metal  theo")
    for r in kept:
        sg = r["spacegroup"] or "?"
        print(f"{r['source']:<16}{r['formula']:<12}{r['mp_id']:<12}{r['chemsys']:<10}"
              f"{r['formation_energy_per_atom_eV']:>11.4f}  {sg:<12}{r['nsites']:<8}"
              f"{str(r['is_metal']):<10}{r['theoretical']}")


if __name__ == "__main__":
    main()
