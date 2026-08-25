"""Full case-1 (decomposition-to-elements) reaction list SI table: every
compound in analysis/reaction_analysis_case1_full.csv (517 reactions, the
project's whole case-1 population -- NOT the family=="extension" subset
appendix_extension_{endo,exo}bondic_{lang}.tex already cover), with the
balanced reaction (pymatgen.analysis.reaction_calculator.Reaction, same
convention as reaction_icohp.py's own reaction_delta_icohp()), the SPACE
GROUP of the compound AND of every element reference structure on the
products side, the ViabilityLabel (analysis/case1_viability.csv), Delta
(ICOHP)/atom (mission #5) and Delta(ICOHP antibonding)/atom (mission
#4b, "--" where the antibonding descriptor lacks COHPCAR.lobster
coverage for the compound or one of its element references, see
analysis/compute_delta_antibonding_case1.py's docstring).

Writes appendix_full_reactions_{fr,en}.tex.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.core import Composition

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent
sys.path.insert(0, str(REPO_ROOT / "analysis"))
from compute_reaction_icohp_case1 import ELEMENT_REFERENCE, STRUCTURES_ROOT  # noqa: E402

REACTION_CSV = REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv"
ANTIBOND_CSV = REPO_ROOT / "analysis" / "delta_antibonding_case1.csv"
VIABILITY_CSV = REPO_ROOT / "analysis" / "case1_viability.csv"
ELEMENTS_CSV = REPO_ROOT / "analysis" / "element_reference_energetics.csv"

_SUBSCRIPT_RE = re.compile(r"(?<=[A-Za-z\)])(\d+)")


def _subscript(text: str) -> str:
    return _SUBSCRIPT_RE.sub(r"$_{\1}$", text)


def _esc(s: str) -> str:
    return str(s).replace("_", "\\_").replace("#", "\\#")


def _star_markup(theoretical, energy_above_hull) -> str:
    if theoretical:
        return ""
    if energy_above_hull is not None and pd.notna(energy_above_hull) and energy_above_hull <= 1e-6:
        return "$^*$"
    return "\\textcolor{red}{$^*$}"


def _fmt(x, digits=4) -> str:
    if x is None or pd.isna(x):
        return "--"
    return f"{x:.{digits}f}"


def _fmt_fr(x, digits=4) -> str:
    s = _fmt(x, digits)
    return s if s == "--" else s.replace(".", "{,}").replace("-", "$-$")


def _elements_formula_map() -> dict[str, str]:
    """element symbol -> its own reference structure's formula (e.g. O -> O2)."""
    out = {}
    for el, dirname in ELEMENT_REFERENCE.items():
        meta_path = STRUCTURES_ROOT / dirname / "mp_metadata.json"
        if meta_path.exists():
            out[el] = json.loads(meta_path.read_text()).get("formula", el)
        else:
            out[el] = el
    return out


def _balanced_reaction(formula: str, el_formula: dict[str, str]):
    """Returns (reduced_reactant, [(product_composition, element_symbol, coeff), ...]) or None."""
    try:
        reactant = Composition(formula)
    except Exception:
        return None
    elements = [str(e) for e in reactant.elements]
    if len(elements) < 2 or any(e not in ELEMENT_REFERENCE for e in elements):
        return None
    reduced_reactant, _ = reactant.get_reduced_composition_and_factor()
    product_comps = []
    for e in elements:
        pc = Composition(el_formula.get(e, e))
        reduced_pc, _ = pc.get_reduced_composition_and_factor()
        product_comps.append((reduced_pc, e))
    try:
        rxn = Reaction([reduced_reactant], [pc for pc, _ in product_comps])
    except Exception:
        return None
    products = [(pc, e, abs(rxn.get_coeff(pc))) for pc, e in product_comps]
    return reduced_reactant, products


VIABILITY_LABEL_TEXT = {
    "fr": {
        "stable_on_hull": "stable", "metastable_viable": "métastable viable",
        "unstable_nonexistent": "instable/inexistant", "insufficient_data": "donnée insuffisante",
    },
    "en": {
        "stable_on_hull": "stable", "metastable_viable": "metastable viable",
        "unstable_nonexistent": "unstable/nonexistent", "insufficient_data": "insufficient data",
    },
}

TEXT = {
    "fr": {
        "header": ["Composé", "mp-id", "Groupe d'espace", "Réaction (produits et leur groupe d'espace)",
                   "Viabilité", "$\\Delta$(ICOHP)/at.", "$\\Delta$(ICOHP antil.)/at."],
        "caption": (
            "Les 517 réactions case-1 (décomposition en éléments) de la population complète du projet "
            "(compléments aux Tableaux~H/I, restreints à la campagne \\og{}extension\\fg{}, 208 composés) --- "
            "composé (formule\\'e\\'etoil\\'ee, mp-id, groupe d'espace), réaction équilibrée avec le groupe "
            "d'espace de chaque référence élémentaire entre parenthèses, étiquette de viabilité "
            "(\\texttt{classify\\_viability()}), $\\Delta$(ICOHP)/atome (mission~5) et "
            "$\\Delta$(ICOHP antiliant)/atome (mission~4b, \\og{}--\\fg{} si la trace COHP complète manque "
            "pour le composé ou l'une de ses références)."
        ),
    },
    "en": {
        "header": ["Compound", "mp-id", "Space group", "Reaction (products and their space group)",
                   "Viability", "$\\Delta$(ICOHP)/at.", "$\\Delta$(ICOHP antib.)/at."],
        "caption": (
            "The project's complete 517-reaction case-1 (decomposition-to-elements) population (a "
            "superset of Tables~H/I, which are restricted to the ``extension'' campaign, 208 compounds) --- "
            "compound (starred formula, mp-id, space group), balanced reaction with each elemental "
            "reference's space group in parentheses, viability label (\\texttt{classify\\_viability()}), "
            "$\\Delta$(ICOHP)/atom (mission~5) and $\\Delta$(ICOHP antibonding)/atom (mission~4b, "
            "``--'' where the full COHP trace is missing for the compound or one of its references)."
        ),
    },
}


def _reaction_string_with_spacegroups(formula: str, el_formula: dict[str, str], el_spacegroup: dict[str, str]) -> str | None:
    balanced = _balanced_reaction(formula, el_formula)
    if balanced is None:
        return None
    reduced_reactant, products = balanced
    parts = []
    for pc, e, coeff in products:
        coeff_str = "" if abs(coeff - 1.0) < 1e-6 else (f"{coeff:.4g}\\,")
        sg = el_spacegroup.get(e)
        sg_str = sg if sg and pd.notna(sg) else "--"
        parts.append(f"{coeff_str}{_subscript(_esc(pc.reduced_formula))} ({_esc(sg_str)})")
    return f"{_subscript(_esc(reduced_reactant.reduced_formula))} $\\to$ " + " $+$ ".join(parts)


def main() -> None:
    rxn_df = pd.read_csv(REACTION_CSV)
    antibond_df = pd.read_csv(ANTIBOND_CSV)[["compound_id", "delta_icohp_antibond"]]
    via_df = pd.read_csv(VIABILITY_CSV)[["compound_id", "viability_label"]]
    elem_df = pd.read_csv(ELEMENTS_CSV)

    el_formula = _elements_formula_map()
    el_spacegroup = dict(zip(elem_df["element"], elem_df["spacegroup"]))

    m = rxn_df.merge(antibond_df, on="compound_id", how="left").merge(via_df, on="compound_id", how="left")

    spacegroups = {}
    for cid in m["compound_id"]:
        meta_path = STRUCTURES_ROOT / cid / "mp_metadata.json"
        if meta_path.exists():
            meta = json.loads(meta_path.read_text())
            spacegroups[cid] = meta.get("spacegroup") or meta.get("expected_spacegroup")
        else:
            spacegroups[cid] = None
    m["compound_spacegroup"] = m["compound_id"].map(spacegroups)

    m = m.sort_values("formula").reset_index(drop=True)

    for lang in ("fr", "en"):
        L = TEXT[lang]
        vlabel_map = VIABILITY_LABEL_TEXT[lang]
        fmt = _fmt_fr if lang == "fr" else _fmt
        body = []
        n_no_reaction = 0
        for _, r in m.iterrows():
            rxn_str = _reaction_string_with_spacegroups(r["formula"], el_formula, el_spacegroup)
            if rxn_str is None:
                n_no_reaction += 1
                continue
            star = _star_markup(r["theoretical"], r["energy_above_hull_eV_per_atom"])
            formula_disp = _subscript(_esc(r["formula"])) + star
            sg = r["compound_spacegroup"] if pd.notna(r["compound_spacegroup"]) else "--"
            vlabel = vlabel_map.get(r["viability_label"], "--") if pd.notna(r["viability_label"]) else "--"
            body.append(
                f"{formula_disp} & {_esc(r['mp_id']) if pd.notna(r['mp_id']) else '--'} & {_esc(sg)} & "
                f"{rxn_str} & {vlabel} & {fmt(r['delta_per_atom_eV'])} & {fmt(r['delta_icohp_antibond'])} \\\\"
            )
        ncol = len(L["header"])
        lines = [
            "\\begin{longtable}{@{}lllp{5.5cm}lrr@{}}",
            f"\\caption{{{L['caption']}}}\\\\",
            "\\toprule",
            " & ".join(L["header"]) + " \\\\",
            "\\midrule",
            "\\endfirsthead",
            "\\toprule",
            " & ".join(L["header"]) + " \\\\",
            "\\midrule",
            "\\endhead",
            "\\bottomrule",
            "\\endfoot",
            "\\bottomrule",
            "\\endlastfoot",
            *body,
            "\\end{longtable}",
        ]
        (HERE / f"appendix_full_reactions_{lang}.tex").write_text("\n".join(lines) + "\n")
        print(f"{lang}: {len(body)} reactions written (skipped {n_no_reaction} with no balanceable reaction)")


if __name__ == "__main__":
    main()
