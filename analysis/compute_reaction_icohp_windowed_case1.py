"""Reaction Delta(ICOHP) restricted to a window near the Fermi level --
the extension proposed in the report's Conclusion: "a reaction
Delta(ICOHP) restricted to a window near the Fermi level, modeled on the
antibonding population... rather than on the fully integrated ICOHP."

Distinct from BOTH existing quantities:
  - reaction_icohp.reaction_delta_icohp() (mission #5): full occupied-range
    ICOHP sum (ICOHPLIST.lobster), no energy restriction.
  - compute_delta_antibonding_case1.py (this session): ANTIBONDING-ONLY
    part of the window, atomic-fraction-weighted average (intensive
    convention, since the antibonding-population metric itself is
    intensive -- built from LOBSTER's "average" trace, not summed).

This one: the SIGNED (bonding+antibonding, not clipped) ICOHP integrated
only within the same one-sided near-E_F window
(cohp_extraction.icohp_windowed_per_atom(), summed over every bond
label -- extensive, like the full ICOHP sum), combined with proper
reaction stoichiometry (pymatgen's Reaction balancer) the SAME way
reaction_icohp.reaction_delta_icohp() does -- not an atomic-fraction
average, because this quantity IS extensive.

Sign convention: reactant-minus-predicted-products, same as
reaction_icohp.py (NOT reaction_analysis's products-minus-reactants) --
this reuses reaction_icohp.py's own Reaction-balancing machinery
directly, so the sign must match it, not the other module.

Needs is_metal for the reactant AND every elemental-reference product
(E_ref selection is per-compound) -- sourced via
compute_icobi_antibonding_all.build_is_metal_map(), NOT
meta.get("is_metal") directly (186 main-campaign compounds' own
mp_metadata.json lacks it; this is the same gap fixed twice already this
session, see project memory).

Writes analysis/reaction_icohp_windowed_case1.csv.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import pandas as pd
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.core import Composition

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).parent))

import cohp_extraction as ce  # noqa: E402
import compute_reaction_icohp_case1 as ri_case1  # noqa: E402
from compute_icobi_antibonding_all import build_is_metal_map  # noqa: E402

STRUCTURES_ROOT = ri_case1.STRUCTURES_ROOT
ELEMENT_REFERENCE = ri_case1.ELEMENT_REFERENCE
OUT_CSV = Path(__file__).parent / "reaction_icohp_windowed_case1.csv"
OUT_JSON = Path(__file__).parent / "reaction_icohp_windowed_case1.json"

DELTA_E = 1.0  # eV, matches the antibonding-population metric's primary window


def is_computed_windowed(compound_dir: Path) -> bool:
    return (
        (compound_dir / "COHPCAR.lobster").exists()
        and (compound_dir / "CONTCAR").exists()
        and (compound_dir / "vasprun.xml").exists()
    )


def reaction_delta_icohp_windowed(
    reactant_dir: Path, reactant_is_metal: bool,
    product_dirs: list[Path], product_is_metal: list[bool],
    delta_e: float,
) -> dict:
    """Same balancing/scaling logic as reaction_icohp.reaction_delta_icohp(),
    with icohp_windowed_per_atom() substituted for icohp_per_atom()
    throughout -- see that function's docstring for the k-scaling and
    sign-convention explanation, unchanged here."""
    reactant = ce.icohp_windowed_per_atom(reactant_dir, is_metal=reactant_is_metal, delta_e=delta_e)
    products = [
        ce.icohp_windowed_per_atom(d, is_metal=im, delta_e=delta_e)
        for d, im in zip(product_dirs, product_is_metal)
    ]

    reactant_comp = reactant["composition"]
    reduced_reactant, factor_reactant = reactant_comp.get_reduced_composition_and_factor()
    product_reduced = [p["composition"].get_reduced_composition_and_factor() for p in products]

    rxn = Reaction([reduced_reactant], [rc for rc, _ in product_reduced])
    k = factor_reactant / abs(rxn.get_coeff(reduced_reactant))

    predicted_total = 0.0
    per_product = []
    for p, (rc, _) in zip(products, product_reduced):
        coeff_formula_units = k * rxn.get_coeff(rc)
        atoms_per_fu = rc.num_atoms
        contribution = coeff_formula_units * atoms_per_fu * p["icohp_windowed_per_atom"]
        predicted_total += contribution
        per_product.append({
            "compound_dir": p["compound_dir"],
            "reduced_formula": rc.reduced_formula,
            "coeff_formula_units": coeff_formula_units,
            "icohp_windowed_per_atom": p["icohp_windowed_per_atom"],
            "contribution_to_total": contribution,
        })

    delta_total = reactant["icohp_windowed_total"] - predicted_total
    delta_per_atom = delta_total / reactant["n_sites"]

    return {
        "reactant_dir": str(reactant_dir),
        "reactant_formula": reactant_comp.reduced_formula,
        "reactant_n_sites": reactant["n_sites"],
        "reactant_icohp_windowed_total": reactant["icohp_windowed_total"],
        "reaction_string": str(rxn),
        "products": per_product,
        "predicted_icohp_windowed_total_from_products": predicted_total,
        "delta_icohp_windowed_total": delta_total,
        "delta_icohp_windowed_per_atom": delta_per_atom,
        "delta_e": delta_e,
    }


def main() -> None:
    is_metal_map = build_is_metal_map()

    results = {}
    rows = []
    n_ok = n_skipped_element = n_skipped_single = n_not_ready = n_no_is_metal = n_failed = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
            meta_path = compound_dir / "mp_metadata.json"
            if not meta_path.exists() or not is_computed_windowed(compound_dir):
                continue
            compound_id = compound_dir.name
            meta = json.loads(meta_path.read_text())
            formula = meta.get("formula")
            if not formula:
                continue
            try:
                comp = Composition(formula)
            except Exception:
                continue
            elements = [str(e) for e in comp.elements]

            if len(elements) == 1:
                n_skipped_single += 1
                continue

            missing_refs = [e for e in elements if e not in ELEMENT_REFERENCE]
            if missing_refs:
                n_skipped_element += 1
                continue

            ref_dirs = [STRUCTURES_ROOT / ELEMENT_REFERENCE[e] for e in elements]
            not_ready = [d.name for d in ref_dirs if not is_computed_windowed(d)]
            if not_ready:
                n_not_ready += 1
                continue

            reactant_is_metal = is_metal_map.get(compound_id)
            product_is_metal = [is_metal_map.get(d.name) for d in ref_dirs]
            if reactant_is_metal is None or any(im is None for im in product_is_metal):
                n_no_is_metal += 1
                continue

            try:
                r = reaction_delta_icohp_windowed(
                    compound_dir, reactant_is_metal, ref_dirs, product_is_metal, DELTA_E
                )
                results[compound_id] = r
                n_ok += 1
                rows.append({
                    "compound_id": compound_id,
                    "mp_id": meta.get("mp_id"),
                    "formula": formula,
                    "family": meta.get("family"),
                    "bond_type": meta.get("bond_type"),
                    "is_metal": reactant_is_metal,
                    "theoretical": meta.get("theoretical"),
                    "energy_above_hull_eV_per_atom": meta.get("energy_above_hull_eV_per_atom"),
                    "reaction_string": r["reaction_string"],
                    "delta_icohp_windowed_total": r["delta_icohp_windowed_total"],
                    "delta_icohp_windowed_per_atom": r["delta_icohp_windowed_per_atom"],
                })
            except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
                results[compound_id] = {"error": f"{type(exc).__name__}: {exc}"}
                n_failed += 1
                print(f"FAILED {compound_id}: {type(exc).__name__}: {exc}")

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))
    print(
        f"\n{n_ok} ok, skipped_single_element={n_skipped_single}, "
        f"skipped_no_reference={n_skipped_element}, reference_not_ready={n_not_ready}, "
        f"no_is_metal={n_no_is_metal}, failed={n_failed}"
    )

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
