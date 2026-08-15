"""Reaction-ICOHP: the total-bonding-energy (ICOHP) analog of a formation
energy, for three kinds of reactions between compounds already computed in
this project (see analysis/METRIC_DEFINITION_reaction_icohp.md for the full
rationale and validation):

  1. decomposition into elements: AaBb -> a A + b B (standard states)
  2. polymorphs: same composition, different structures, direct comparison
  3. decomposition into another compound + elements (general reaction)

Case 1 and 3 share the same underlying primitive (reaction_delta_icohp) --
case 1 is not special-cased, it is simply the case where every product
happens to be a pure element. Case 2 needs no reaction balancing at all
(icohp_per_atom is directly comparable between same-composition polymorphs).

Core intensive quantity: icohp_per_atom(compound_dir), the sum of every
bond ICOHP in ICOHPLIST.lobster (LOBSTER's own accounting of every
symmetry-inequivalent periodic bond in the cell, unfiltered -- the same
raw source percolation_path.py's icohp_sum column already uses for the
whole-cell case, bond_pair=None) divided by the number of atoms in that
cell. This is a bulk/converged property, assumed transferable between a
compound's own cell and its constituent references' cells, exactly the
same assumption formation_energy_per_atom already makes for total energy.

Sign convention (established throughout this project, see
cohp_extraction.py's sign_convention_check): more negative ICOHP = more
bonding. So a negative delta_icohp_per_atom for a decomposition reaction
means the compound has MORE net bonding per atom than its reference
products -- the ICOHP analog of a negative (stabilizing) formation energy.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Sequence

from pymatgen.core import Structure
from pymatgen.io.lobster.outputs import Icohplist


def icohp_per_atom(compound_dir: Path) -> Dict[str, object]:
    """Sum every bond in ICOHPLIST.lobster (unfiltered -- every
    symmetry-inequivalent periodic bond LOBSTER reports for the cell) and
    normalize per atom. Uses CONTCAR (relaxed structure), same file
    cohp_extraction.py's CompleteCohp loader uses, for the atom count."""
    icohplist = Icohplist(filename=str(compound_dir / "ICOHPLIST.lobster"))
    values = [
        icohplist.icohpcollection.get_icohp_by_label(label)
        for label in icohplist.icohplist
    ]
    structure = Structure.from_file(compound_dir / "CONTCAR")
    n_sites = len(structure)
    icohp_total = sum(values)
    return {
        "compound_dir": str(compound_dir),
        "n_bonds": len(values),
        "n_sites": n_sites,
        "icohp_total": icohp_total,
        "icohp_per_atom": icohp_total / n_sites,
        "composition": structure.composition,
    }


def reaction_delta_icohp(reactant_dir: Path, product_dirs: Sequence[Path]) -> Dict[str, object]:
    """Balance reactant -> products by composition alone (pymatgen's
    Reaction, reduced-formula basis) and compute the ICOHP-of-reaction on
    the reactant's own cell (per atom), the same normalization
    formation_energy_per_atom uses.

    Works identically for case 1 (every product a pure element) and case 3
    (a mix of another compound and elements) -- whatever product_dirs are
    supplied get balanced together, no special-casing by product type.

    delta_icohp_total = icohp_total(reactant's own cell)
                         - sum_p [ k * rxn.get_coeff(product_p) ]
                                 * icohp_per_atom(product_p)
                                 * atoms_per_formula_unit(product_p)
    delta_icohp_per_atom = delta_icohp_total / n_sites(reactant's own cell)

    where k scales the balanced reaction (normally quoted per 1 reduced
    formula unit of reactant) up to the reactant's ACTUAL cell composition
    (n_sites(reactant) atoms), so the products' coefficients are directly
    in units of "formula units consumed per reactant cell".
    """
    from pymatgen.analysis.reaction_calculator import Reaction

    reactant = icohp_per_atom(reactant_dir)
    products = [icohp_per_atom(d) for d in product_dirs]

    reactant_comp = reactant["composition"]
    reduced_reactant, factor_reactant = reactant_comp.get_reduced_composition_and_factor()
    product_reduced = [p["composition"].get_reduced_composition_and_factor() for p in products]

    rxn = Reaction([reduced_reactant], [rc for rc, _ in product_reduced])
    # k scales "1 reduced-formula-unit of reactant" up to the reactant's
    # actual cell (n_sites atoms) -- factor_reactant = n_sites / (atoms per
    # reduced formula unit).
    k = factor_reactant / abs(rxn.get_coeff(reduced_reactant))

    predicted_total = 0.0
    per_product = []
    for p, (rc, _) in zip(products, product_reduced):
        coeff_formula_units = k * rxn.get_coeff(rc)
        atoms_per_fu = rc.num_atoms
        contribution = coeff_formula_units * atoms_per_fu * p["icohp_per_atom"]
        predicted_total += contribution
        per_product.append({
            "compound_dir": p["compound_dir"],
            "reduced_formula": rc.reduced_formula,
            "coeff_formula_units": coeff_formula_units,
            "icohp_per_atom": p["icohp_per_atom"],
            "contribution_to_total": contribution,
        })

    delta_total = reactant["icohp_total"] - predicted_total
    delta_per_atom = delta_total / reactant["n_sites"]

    return {
        "reactant_dir": str(reactant_dir),
        "reactant_formula": reactant_comp.reduced_formula,
        "reactant_n_sites": reactant["n_sites"],
        "reactant_icohp_total": reactant["icohp_total"],
        "reaction_string": str(rxn),
        "products": per_product,
        "predicted_icohp_total_from_products": predicted_total,
        "delta_icohp_total": delta_total,
        "delta_icohp_per_atom": delta_per_atom,
    }


def compare_polymorphs(compound_dirs: Sequence[Path]) -> Dict[str, object]:
    """Same composition, different structures -- no reaction balancing
    needed, icohp_per_atom is directly comparable. Returns a table sorted
    by icohp_per_atom (most negative = most bonding first)."""
    rows = []
    formulas = set()
    for d in compound_dirs:
        r = icohp_per_atom(d)
        formulas.add(r["composition"].reduced_formula)
        rows.append({
            "compound_dir": str(d),
            "reduced_formula": r["composition"].reduced_formula,
            "n_sites": r["n_sites"],
            "icohp_per_atom": r["icohp_per_atom"],
        })
    if len(formulas) > 1:
        raise ValueError(
            f"compare_polymorphs got mixed compositions {formulas} -- "
            "these are not polymorphs of each other."
        )
    rows.sort(key=lambda r: r["icohp_per_atom"])
    for row in rows:
        row["delta_icohp_per_atom_vs_most_bonding"] = (
            row["icohp_per_atom"] - rows[0]["icohp_per_atom"]
        )
    return {"reduced_formula": formulas.pop(), "polymorphs": rows}
