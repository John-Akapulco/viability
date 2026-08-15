"""Build a CompoundEntry from a compound directory (ICOHPLIST.lobster +
optional ICOBILIST.lobster + a structure file) already computed elsewhere
in this project.

Empirical finding this module's summation depends on (verified by
tests/test_parse_lobster.py on tests/fixtures/, reusing the exact
examples/dataset/ toy compounds percolation_path.py's own tests already
use rather than inventing a separate fixture set): **LOBSTER lists each
symmetry-inequivalent periodic bond exactly once in ICOHPLIST.lobster,
not once per direction.** A bond between atom i and atom j with
translation (nx,ny,nz) never has a separate reverse entry (atom j, atom
i, (-nx,-ny,-nz)) in the same file. This is why IcohpSummary.sum_total_eV
is a plain unfiltered sum over every label with no /2 factor --
`percolation_path.py`'s own graph builder adds the reverse-direction edge
itself (see its module docstring) precisely because the file does not
already contain it; `reaction_icohp.py`'s icohp_per_atom makes the same
unfiltered-sum assumption. This module follows the same convention for
consistency across the whole project, now with an explicit regression
test backing it instead of an inherited assumption.

Atom-species-pair extraction (for `by_bond_type`) reuses the exact
pymatgen access pattern already established in percolation_path.py's
`_ingest()` (see that module's comment): `IcohpCollection` has no public
accessor for atom1/atom2 alongside the value, so `._icohplist` (label ->
IcohpValue) plus `IcohpValue.as_dict()` (public) is used instead of
poking at private attributes directly.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

from pymatgen.core import Structure
from pymatgen.io.lobster.outputs import Icohplist
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from reaction_analysis.schema import BondTypeSummary, CompoundEntry, IcohpSummary

_ELEMENT_FROM_LABEL = re.compile(r"^([A-Za-z]+)")


def _element_from_atom_label(label: str) -> str:
    """"Na1" -> "Na", "Cl2" -> "Cl" -- LOBSTER atom labels are always an
    element symbol immediately followed by a 1-based site index."""
    m = _ELEMENT_FROM_LABEL.match(label)
    if not m:
        raise ValueError(f"Cannot extract element symbol from LOBSTER atom label {label!r}")
    return m.group(1)


def _load_structure(compound_dir: Path) -> Structure:
    for name in ("CONTCAR", "POSCAR"):
        candidate = compound_dir / name
        if candidate.exists():
            return Structure.from_file(candidate)
    cifs = sorted(compound_dir.glob("*.cif"))
    if cifs:
        return Structure.from_file(cifs[0])
    raise FileNotFoundError(f"No CONTCAR, POSCAR, or *.cif structure file found in {compound_dir}")


def raw_bond_records(compound_dir: Path, filename: str, are_cobis: bool) -> list[dict]:
    """One dict per symmetry-inequivalent bond label: atom1, atom2 (LOBSTER
    site labels, e.g. "Na1"), translation (tuple[int,int,int]), value
    (spin-summed ICOHP/ICOBI, eV). Same pymatgen access pattern as
    percolation_path.py's _ingest() -- see this module's docstring."""
    parsed = Icohplist(are_cobis=are_cobis, filename=str(compound_dir / filename))
    records = []
    for value in parsed.icohpcollection._icohplist.values():
        d = value.as_dict()
        records.append({
            "atom1": d["atom1"],
            "atom2": d["atom2"],
            "translation": tuple(int(x) for x in d["translation"]),
            "value": value.summed_icohp,
        })
    return records


def check_no_reverse_duplicates(records: list[dict]) -> list[str]:
    """Returns a list of human-readable problems (empty = OK) if any bond
    (atom_i, atom_j, translation) also has a separate reverse-direction
    entry (atom_j, atom_i, -translation) in the same record list -- the
    double-counting failure mode this module's summation assumes does not
    happen (see module docstring)."""
    seen = {(r["atom1"], r["atom2"], r["translation"]) for r in records}
    problems = []
    for r in records:
        reverse_translation = tuple(-x for x in r["translation"])
        reverse_key = (r["atom2"], r["atom1"], reverse_translation)
        if reverse_key in seen:
            problems.append(
                f"possible reverse-direction duplicate: ({r['atom1']}, {r['atom2']}, {r['translation']}) "
                f"and ({r['atom2']}, {r['atom1']}, {reverse_translation}) both present"
            )
    return problems


def _build_summary(records: list[dict], n_sites: int, z: int) -> IcohpSummary:
    n_bonds = len(records)
    sum_total = sum(r["value"] for r in records)

    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        el1 = _element_from_atom_label(r["atom1"])
        el2 = _element_from_atom_label(r["atom2"])
        key = "-".join(sorted((el1, el2)))
        by_pair[key].append(r)

    by_bond_type = {
        key: BondTypeSummary(
            mean_eV=sum(r["value"] for r in recs) / len(recs),
            sum_eV=sum(r["value"] for r in recs),
            n_bonds=len(recs),
        )
        for key, recs in by_pair.items()
    }

    return IcohpSummary(
        sum_total_eV=sum_total,
        sum_per_atom_eV=sum_total / n_sites,
        sum_per_formula_unit_eV=sum_total / z,
        mean_per_bond_eV=sum_total / n_bonds,
        n_bonds=n_bonds,
        by_bond_type=by_bond_type,
    )


def parse_compound_entry(
    compound_dir: Path,
    role: str,
    *,
    compound_id: Optional[str] = None,
    energy_total_eV: Optional[float] = None,
) -> CompoundEntry:
    """Build a CompoundEntry from `compound_dir`. `role` is required and
    never inferred (no chemical-role logic lives in this package -- see
    package docstring); the caller decides whether a given directory is
    playing "target", "element", "hull_neighbor", or "polymorph" in
    whatever Reaction it will be used for.

    `compound_id` defaults to `compound_dir.name` -- this project has no
    `build_candidate_name()`-style dedicated naming function (checked:
    not present anywhere in the current codebase, despite being an
    assumed convention in the original request for this module); every
    other module in this project (compute_reaction_icohp_case1.py,
    compute_antibonding_all.py, ...) already uses the directory name
    directly as compound_id, so this module follows the same actual
    convention rather than a function that does not exist.

    `energy_total_eV` is passthrough-only: nothing in this module parses
    OUTCAR/vasprun.xml (both gitignored project-wide, and most already
    deleted from disk to save space) to obtain it.
    """
    structure = _load_structure(compound_dir)
    n_sites = len(structure)
    composition = structure.composition
    reduced, z = composition.get_reduced_composition_and_factor()
    comp_dict = {str(el): int(round(amt)) for el, amt in composition.get_el_amt_dict().items()}

    sga = SpacegroupAnalyzer(structure)

    icohp_records = raw_bond_records(compound_dir, "ICOHPLIST.lobster", are_cobis=False)
    icohp_summary = _build_summary(icohp_records, n_sites, round(z))

    icobi_summary = None
    icobi_path = compound_dir / "ICOBILIST.lobster"
    if icobi_path.exists():
        icobi_records = raw_bond_records(compound_dir, "ICOBILIST.lobster", are_cobis=True)
        icobi_summary = _build_summary(icobi_records, n_sites, round(z))

    return CompoundEntry(
        compound_id=compound_id or compound_dir.name,
        formula=composition.reduced_formula,
        composition=comp_dict,
        Z=round(z),
        space_group_symbol=sga.get_space_group_symbol(),
        space_group_number=sga.get_space_group_number(),
        role=role,
        energy_total_eV=energy_total_eV,
        energy_per_atom_eV=(energy_total_eV / n_sites) if energy_total_eV is not None else None,
        icohp=icohp_summary,
        icobi=icobi_summary,
        source_path=str(compound_dir),
    )
