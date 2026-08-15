"""CLI entry point for reaction_analysis, same spirit as
percolation_path.py: batch-process every reaction under
--dataset/reactions/*.json against the compound entries under
--dataset/entries/*/entry.json, writing one row per (reaction, metric)
to CSV or JSON. A reaction that fails (unbalanced, missing entry, bad
bond_pair) never stops the batch -- its error is recorded in that row's
`error` field and the batch continues.

entries/<id>/entry.json is expected to already exist (produced ahead of
time by parse_lobster.py + CompoundEntry.model_dump_json()) -- this CLI
does not parse raw LOBSTER output itself. As of this commit, no real
production dataset has been built in this format yet; see the package
docstring.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Optional, Sequence

from reaction_analysis.delta import compute_delta
from reaction_analysis.schema import CompoundEntry, Reaction, ReactionResult

RESULT_FIELDS = [
    "reaction_id", "type", "metric", "bond_pair",
    "delta_per_formula_unit_eV", "delta_per_atom_eV", "delta_per_bond_eV",
    "delta_per_bond_conservative", "warnings", "error",
]


def load_entries(dataset_root: Path) -> dict[str, CompoundEntry]:
    entries: dict[str, CompoundEntry] = {}
    entries_dir = dataset_root / "entries"
    if not entries_dir.exists():
        return entries
    for entry_dir in sorted(entries_dir.iterdir()):
        entry_path = entry_dir / "entry.json"
        if not entry_path.exists():
            continue
        entries[entry_dir.name] = CompoundEntry.model_validate_json(entry_path.read_text())
    return entries


def load_reactions(dataset_root: Path, reaction_type: Optional[str] = None) -> list[Reaction]:
    reactions = []
    reactions_dir = dataset_root / "reactions"
    if not reactions_dir.exists():
        return reactions
    for path in sorted(reactions_dir.glob("*.json")):
        reaction = Reaction.model_validate_json(path.read_text())
        if reaction_type is not None and reaction.type != reaction_type:
            continue
        reactions.append(reaction)
    return reactions


def run(dataset_root: Path, reaction_type: Optional[str] = None) -> list[ReactionResult]:
    entries = load_entries(dataset_root)
    results: list[ReactionResult] = []
    for reaction in load_reactions(dataset_root, reaction_type):
        missing = sorted({
            m.compound_id for m in (*reaction.reactants, *reaction.products)
            if m.compound_id not in entries
        })
        if missing:
            results.append(ReactionResult(
                reaction_id=reaction.reaction_id, type=reaction.type, metric=reaction.metric,
                bond_pair=reaction.bond_pair, delta_per_formula_unit_eV=math.nan,
                delta_per_atom_eV=math.nan, delta_per_bond_eV=math.nan,
                delta_per_bond_conservative=False, warnings=[],
                error=f"missing CompoundEntry(s): {missing}",
            ))
            continue
        results.extend(compute_delta(reaction, entries))
    return results


def _write_csv(results: list[ReactionResult], path: Path) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for r in results:
            row = r.model_dump()
            row["warnings"] = "; ".join(row["warnings"])
            writer.writerow(row)


def _write_json(results: list[ReactionResult], path: Path) -> None:
    path.write_text(json.dumps([r.model_dump() for r in results], indent=2, default=str))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--dataset", type=Path, required=True,
        help="Root of a reactions_dataset/-style tree (must contain entries/ and reactions/).",
    )
    p.add_argument(
        "--reaction-type", dest="reaction_type", default=None,
        choices=["decomposition_to_elements", "decomposition_to_compound_and_elements", "polymorph_transition"],
        help="Restrict to reactions of this type (default: all types).",
    )
    p.add_argument("--output", type=Path, required=True, help="Output path (.csv or .json).")
    p.add_argument("--also-json", type=Path, default=None, help="Additionally write full JSON detail here.")
    p.add_argument("-v", "--verbose", action="count", default=0)
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    results = run(args.dataset, args.reaction_type)

    if args.output.suffix == ".json":
        _write_json(results, args.output)
    else:
        _write_csv(results, args.output)
    if args.also_json is not None:
        _write_json(results, args.also_json)

    if args.verbose:
        n_err = sum(1 for r in results if r.error is not None)
        print(f"Processed {len(results)} result row(s) ({n_err} with an error) -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
