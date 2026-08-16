"""Generate LaTeX appendix fragments for the FULL "extension" family
(163 compounds: extension batches 1-4, download_extension.py through
download_extension4.py -- elemental references, carbon allotropes, TiO2
high-pressure polymorphs, ionic/covalent extension compounds, and the
89-compound extension4 batch), per the user's explicit table-layout
request (2026-08-16): every compound and element listed, reactions
detailed, endobondic/exobondic split into separate tables, hull distance
shown, and experimental compounds marked with a superscript star on
their formula (theoretical compounds unmarked).

Correction (same day): an earlier version of this script scoped to
extension4's 89 JSON-listed compounds only, which silently dropped the
carbon allotropes and TiO2 high-pressure polymorphs (extension batch 2,
2026-08-14) from the report even though the user had asked for them to
be considered -- this version scopes to every mp_metadata.json with
family=="extension" instead, keyed off nothing but that field.

Writes, per language (fr/en), four fragments:
  appendix_extension_elements_{lang}.tex     -- single-element compounds
    (formula*, mp-id/COD-id, space group, E_hull). No Delta-ICOHP/label
    column: an element is not itself a decomposition reaction.
  appendix_extension_compounds_{lang}.tex    -- multi-element compounds
    (formula*, mp-id/COD-id, space group, E_hull, Delta-ICOHP/atom,
    endo/exo label where a case-1 reaction exists; "--" where it
    doesn't, e.g. missing elemental reference).
  appendix_extension_endobondic_{lang}.tex   -- the subset of case-1
    reactions with bonding_label=endobondic (balanced reaction string,
    target formula starred if experimental, E_hull, Delta-ICOHP/atom).
  appendix_extension_exobondic_{lang}.tex    -- same, exobondic subset.

Delta-ICOHP/atom and bonding_label come from reaction_analysis_case1_full.csv
(products-minus-reactants convention, classify.py's sign -- endobondic iff
delta_per_atom_eV>=0). The balanced reaction string comes from
reaction_icohp_case1.csv (reactant-minus-products convention -- only the
text is reused, not its own delta column, to avoid mixing sign
conventions in one table).

Run manually after any dataset change; not part of the automated analysis
pipeline, same convention as gen_appendix.py / gen_appendix_dft.py.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"

TEXT = {
    "fr": {
        "elements_caption": (
            "Éléments et allotropes de la campagne \\texttt{extension} (lots 1--4) : "
            "références utilisées pour la décomposition des composés multi-éléments, "
            "plus les allotropes supplémentaires (carbone, TiO2 n'est pas ici -- "
            "composé binaire, voir table des composés). $^*$ = référence expérimentale "
            "(non marqué = théorique). Pas de $\\Delta$ICOHP : un élément n'est pas "
            "lui-même une réaction de décomposition."
        ),
        "elements_header": ["Formule", "mp-id / COD", "Groupe d'espace", "$E_{hull}$ (eV/at)"],
        "compounds_caption": (
            "Composés cibles et polymorphes de la campagne \\texttt{extension} (lots "
            "1--4, 163 composés dont 89 d'\\texttt{extension4}) -- inclut les polymorphes "
            "haute-pression de TiO2 et tout composé multi-élément des lots précédents. "
            "$^*$ = composé expérimental (non marqué = théorique). $\\Delta$ICOHP par "
            "atome en convention produits $-$ réactifs (\\texttt{reaction\\_analysis}) ; "
            "``--'' si aucune réaction de cas 1 n'a pu être calculée (référence "
            "élémentaire manquante)."
        ),
        "compounds_header": ["Formule", "mp-id / COD", "Groupe d'espace", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)", "Étiquette"],
        "endo_caption": (
            "Réactions de décomposition en éléments (cas 1) classées \\textbf{endobondic} "
            "($\\Delta$ICOHP $\\geq 0$, produits $-$ réactifs) parmi les 163 composés de "
            "la campagne \\texttt{extension}. $^*$ = composé de départ expérimental."
        ),
        "exo_caption": (
            "Réactions de décomposition en éléments (cas 1) classées \\textbf{exobondic} "
            "($\\Delta$ICOHP $< 0$, produits $-$ réactifs) parmi les 163 composés de la "
            "campagne \\texttt{extension}. $^*$ = composé de départ expérimental."
        ),
        "reactions_header": ["Réaction (balancée)", "mp-id / COD", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)"],
        "endo": "endobondic", "exo": "exobondic",
    },
    "en": {
        "elements_caption": (
            "Elements and allotropes of the \\texttt{extension} campaign (batches 1--4): "
            "references used to decompose multi-element compounds, plus additional "
            "carbon allotropes (TiO2's polymorphs are a binary compound, see the "
            "compounds table instead). $^*$ = experimental reference (unmarked = "
            "theoretical). No $\\Delta$ICOHP column: an element is not itself a "
            "decomposition reaction."
        ),
        "elements_header": ["Formula", "mp-id / COD", "Space group", "$E_{hull}$ (eV/at)"],
        "compounds_caption": (
            "Target compounds and polymorphs of the \\texttt{extension} campaign "
            "(batches 1--4, 163 compounds, 89 from \\texttt{extension4}) -- includes "
            "TiO2's high-pressure polymorphs and every multi-element compound from "
            "earlier batches. $^*$ = experimental compound (unmarked = theoretical). "
            "$\\Delta$ICOHP per atom, products $-$ reactants convention "
            "(\\texttt{reaction\\_analysis}); ``--'' where no case-1 reaction could be "
            "computed (missing elemental reference)."
        ),
        "compounds_header": ["Formula", "mp-id / COD", "Space group", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)", "Label"],
        "endo_caption": (
            "Decomposition-into-elements (case 1) reactions labeled \\textbf{endobondic} "
            "($\\Delta$ICOHP $\\geq 0$, products $-$ reactants) among the 163 compounds "
            "of the \\texttt{extension} campaign. $^*$ = experimental starting compound."
        ),
        "exo_caption": (
            "Decomposition-into-elements (case 1) reactions labeled \\textbf{exobondic} "
            "($\\Delta$ICOHP $< 0$, products $-$ reactants) among the 163 compounds of "
            "the \\texttt{extension} campaign. $^*$ = experimental starting compound."
        ),
        "reactions_header": ["Reaction (balanced)", "mp-id / COD", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)"],
        "endo": "endobondic", "exo": "exobondic",
    },
}


def _esc(s: str) -> str:
    return s.replace("_", "\\_").replace("#", "\\#")


def _star(formula: str, theoretical: bool) -> str:
    f = _esc(formula)
    return f if theoretical else f + "$^*$"


def _id_label(meta: dict) -> str:
    mp_id = meta.get("mp_id") or meta.get("material_id")
    if mp_id:
        return mp_id
    cod_id = meta.get("cod_id")
    return f"COD:{cod_id}" if cod_id else "--"


def load_all_rows():
    element_rows, compound_rows = [], []
    for d in sorted(STRUCTURES_ROOT.iterdir()):
        meta_path = d / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("family") != "extension":
            continue
        formula = meta.get("formula")
        if not formula:
            continue
        try:
            n_elements = len(Composition(formula).elements)
        except Exception:
            continue

        row = {
            "compound_id": d.name, "formula": formula, "id_label": _id_label(meta),
            "theoretical": bool(meta.get("theoretical")), "spacegroup": meta.get("spacegroup") or "--",
            "eah": meta.get("energy_above_hull_eV_per_atom"),
        }
        (element_rows if n_elements == 1 else compound_rows).append(row)

    ri = pd.read_csv(REPO_ROOT / "analysis" / "reaction_icohp_case1.csv")[["compound_id", "reaction_string"]]
    ra = pd.read_csv(REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv")[["compound_id", "delta_per_atom_eV"]]
    reaction_lookup = ri.merge(ra, on="compound_id", how="inner").set_index("compound_id").to_dict("index")

    for row in compound_rows:
        r = reaction_lookup.get(row["compound_id"])
        row["reaction_string"] = r["reaction_string"] if r else None
        row["delta_per_atom_eV"] = r["delta_per_atom_eV"] if r else None

    return element_rows, compound_rows


def _longtable(caption: str, header: list[str], colspec: str, body_lines: list[str]) -> str:
    lines = [
        f"\\begin{{longtable}}{{@{{}}{colspec}@{{}}}}",
        f"\\caption{{{caption}}}\\\\",
        "\\toprule",
        " & ".join(header) + " \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        " & ".join(header) + " \\\\",
        "\\midrule",
        "\\endhead",
        "\\bottomrule",
        "\\endfoot",
        *body_lines,
        "\\end{longtable}",
    ]
    return "\n".join(lines) + "\n"


def write_elements(lang: str, rows: list[dict]) -> None:
    L = TEXT[lang]
    body = []
    for r in sorted(rows, key=lambda r: r["formula"]):
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        body.append(f"{_star(r['formula'], r['theoretical'])} & {_esc(r['id_label'])} & {_esc(r['spacegroup'])} & {eah} \\\\")
    out = _longtable(L["elements_caption"], L["elements_header"], "llll", body)
    (HERE / f"appendix_extension_elements_{lang}.tex").write_text(out)


def write_compounds(lang: str, rows: list[dict]) -> None:
    L = TEXT[lang]
    body = []
    for r in sorted(rows, key=lambda r: r["formula"]):
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        has_reaction = r["delta_per_atom_eV"] is not None
        delta = f"{r['delta_per_atom_eV']:.4f}" if has_reaction else "--"
        label = (L["endo"] if r["delta_per_atom_eV"] >= 0 else L["exo"]) if has_reaction else "--"
        body.append(
            f"{_star(r['formula'], r['theoretical'])} & {_esc(r['id_label'])} & {_esc(r['spacegroup'])} "
            f"& {eah} & {delta} & {label} \\\\"
        )
    out = _longtable(L["compounds_caption"], L["compounds_header"], "lllrrl", body)
    (HERE / f"appendix_extension_compounds_{lang}.tex").write_text(out)


def _starred_reaction_string(r: dict) -> str:
    rxn = r["reaction_string"].replace("->", "$\\rightarrow$")
    if not r["theoretical"]:
        rxn = re.sub(r"\b" + re.escape(r["formula"]) + r"\b", r["formula"] + "$^*$", rxn, count=1)
    return rxn


def write_reactions(lang: str, rows: list[dict], endobondic: bool) -> None:
    L = TEXT[lang]
    subset = [r for r in rows if r["delta_per_atom_eV"] is not None and (r["delta_per_atom_eV"] >= 0) == endobondic]
    body = []
    for r in sorted(subset, key=lambda r: r["formula"]):
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        body.append(
            f"{_starred_reaction_string(r)} & {_esc(r['id_label'])} & {eah} & {r['delta_per_atom_eV']:.4f} \\\\"
        )
    caption = L["endo_caption"] if endobondic else L["exo_caption"]
    out = _longtable(caption, L["reactions_header"], "lrrr", body)
    suffix = "endobondic" if endobondic else "exobondic"
    (HERE / f"appendix_extension_{suffix}_{lang}.tex").write_text(out)


def main():
    element_rows, compound_rows = load_all_rows()
    print(f"{len(element_rows)} element/allotrope rows")
    print(f"{len(compound_rows)} compound/polymorph rows")
    n_with_reaction = sum(1 for r in compound_rows if r["delta_per_atom_eV"] is not None)
    n_endo = sum(1 for r in compound_rows if r["delta_per_atom_eV"] is not None and r["delta_per_atom_eV"] >= 0)
    n_exo = sum(1 for r in compound_rows if r["delta_per_atom_eV"] is not None and r["delta_per_atom_eV"] < 0)
    print(f"with a case-1 reaction: {n_with_reaction} (endobondic={n_endo}, exobondic={n_exo})")

    for lang in ("fr", "en"):
        write_elements(lang, element_rows)
        write_compounds(lang, compound_rows)
        write_reactions(lang, compound_rows, endobondic=True)
        write_reactions(lang, compound_rows, endobondic=False)

    for stale in list(HERE.glob("appendix_extension4_*.tex")):
        stale.unlink()
        print(f"removed superseded {stale.name}")

    print("wrote appendix_extension_{elements,compounds,endobondic,exobondic}_{fr,en}.tex")


if __name__ == "__main__":
    main()
