"""Generate LaTeX appendix fragments for the FULL "extension" family
(163 compounds: extension batches 1-4, download_extension.py through
download_extension4.py -- elemental references, carbon allotropes, TiO2
high-pressure polymorphs, ionic/covalent extension compounds, and the
89-compound extension4 batch). Batch 5 (maxhull_binaries_stress_test,
2026-08-16, still computing) is deliberately excluded here -- these
fragments describe the already-analyzed 163, not a moving target; rerun
this script once that batch finishes to fold it in as its own update.

Per the user's table-layout requests (2026-08-16 and 2026-08-16 follow-up):
  - every compound and element listed, reactions detailed, endobondic/
    exobondic split into separate tables, hull distance shown;
  - stoichiometric numbers in every formula are LaTeX-subscripted
    (_subscript()), in both table cells and reaction strings;
  - every experimental entry is starred; the star is BLACK ($^*$) if the
    entry is the thermodynamically stable phase (energy_above_hull <=
    ~0) and RED (\\textcolor{red}{$^*$}) if it is an experimentally
    known but metastable polymorph/allotrope (energy_above_hull > 0);
    theoretical (never-synthesized) entries are never starred;
  - any formula (element or compound) with >=2 computed entries in this
    family is a polymorph/allotrope group -- these are pulled OUT of the
    plain elements/compounds tables and listed together in one dedicated
    table (appendix_extension_polymorphs_{lang}.tex), because the
    natural comparison for a polymorph group is polymorph-to-polymorph,
    not the decomposition-to-elements reaction the compounds/endo/exo
    tables are built around. The endobondic/exobondic reaction tables
    are NOT filtered this way -- every compound's own case-1 reaction
    (grouped or not) is still a valid, independent result.

Writes, per language (fr/en), five fragments:
  appendix_extension_elements_{lang}.tex     -- single-element compounds
    with no allotrope sibling in this dataset (formula*, mp-id/COD-id,
    space group, E_hull). No Delta-ICOHP/label column: an element is not
    itself a decomposition reaction.
  appendix_extension_compounds_{lang}.tex    -- multi-element compounds
    with no polymorph sibling in this dataset (formula*, mp-id/COD-id,
    space group, E_hull, Delta-ICOHP/atom, endo/exo label where a case-1
    reaction exists; "--" where it doesn't, e.g. missing elemental
    reference).
  appendix_extension_polymorphs_{lang}.tex   -- every element or compound
    that has >=2 entries sharing a formula, grouped by formula (a
    \\midrule separates groups), same 4 structural columns as the
    elements table -- no Delta-ICOHP column, deliberately (see above).
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
from collections import Counter
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"

EXCLUDED_BATCHES = {"maxhull_binaries_stress_test"}

# A digit run immediately after a letter or ")" is a stoichiometric
# subscript (TiO2 -> Ti, O2; Ca3N2 -> Ca3, N2); a digit run at the start
# of a token or after a space/"." is a standalone reaction coefficient
# (e.g. "0.3333 Ca3N2", "3 N2") and must NOT be subscripted.
_SUBSCRIPT_RE = re.compile(r"(?<=[A-Za-z\)])(\d+)")
_EPS_HULL = 1e-6


def _subscript(text: str) -> str:
    return _SUBSCRIPT_RE.sub(r"$_{\1}$", text)


def _esc(s: str) -> str:
    return s.replace("_", "\\_").replace("#", "\\#")


def _star_markup(theoretical: bool, energy_above_hull: float | None) -> str:
    if theoretical:
        return ""
    if energy_above_hull is not None and energy_above_hull <= _EPS_HULL:
        return "$^*$"
    return "\\textcolor{red}{$^*$}"


def _star(formula: str, theoretical: bool, energy_above_hull: float | None) -> str:
    return _subscript(_esc(formula)) + _star_markup(theoretical, energy_above_hull)


def _id_label(meta: dict) -> str:
    mp_id = meta.get("mp_id") or meta.get("material_id")
    if mp_id:
        return mp_id
    cod_id = meta.get("cod_id")
    return f"COD:{cod_id}" if cod_id else "--"


TEXT = {
    "fr": {
        "elements_caption": (
            "Éléments et allotropes sans polymorphe apparié de la campagne "
            "\\texttt{extension} : références utilisées pour la "
            "décomposition des composés multi-éléments. Les allotropes ayant "
            "un ou plusieurs autres allotropes calculés dans ce jeu de "
            "données (p.~ex. le carbone) sont regroupés dans le "
            "Tableau~\\ref{tab:ext-polymorphs}, pas ici. $^*$ noir = référence "
            "expérimentale stable (phase thermodynamique) ; "
            "\\textcolor{red}{$^*$} rouge = référence expérimentale "
            "métastable ; non marqué = théorique. Pas de $\\Delta$ICOHP : un "
            "élément n'est pas lui-même une réaction de décomposition."
        ),
        "elements_header": ["Formule", "mp-id / COD", "Groupe d'espace", "$E_{hull}$ (eV/at)"],
        "compounds_caption": (
            "Composés cibles sans polymorphe apparié de la campagne "
            "\\texttt{extension}. Les formules ayant plusieurs "
            "entrées dans ce jeu de données sont regroupées dans le "
            "Tableau~\\ref{tab:ext-polymorphs}, pas ici. $^*$ noir = composé "
            "expérimental stable (phase thermodynamique) ; "
            "\\textcolor{red}{$^*$} rouge = composé expérimental métastable ; "
            "non marqué = théorique. $\\Delta$ICOHP par atome en convention "
            "produits $-$ réactifs (\\texttt{reaction\\_analysis}) ; ``--'' si "
            "aucune réaction de cas 1 n'a pu être calculée (référence "
            "élémentaire manquante)."
        ),
        "compounds_header": ["Formule", "mp-id / COD", "Groupe d'espace", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)", "Étiquette"],
        "polymorphs_caption": (
            "Éléments (allotropes) et composés (polymorphes) ayant "
            "\\textbf{plusieurs entrées de même formule} dans la campagne "
            "\\texttt{extension} (39 groupes, 84 entrées) --- regroupés ici "
            "plutôt que dispersés dans les tableaux d'éléments/composés, car "
            "la comparaison naturelle pour un groupe de polymorphes est "
            "polymorphe-à-polymorphe, pas la réaction de décomposition en "
            "éléments (d'où l'absence de colonne $\\Delta$ICOHP). Un trait "
            "sépare chaque groupe ; au sein d'un groupe, tri par $E_{hull}$ "
            "croissant. $^*$ noir = expérimental stable ; "
            "\\textcolor{red}{$^*$} rouge = expérimental métastable ; non "
            "marqué = théorique."
        ),
        "endo_caption": (
            "Réactions de décomposition en éléments (cas 1) classées \\textbf{endobondic} "
            "($\\Delta$ICOHP $\\geq 0$, produits $-$ réactifs) parmi les 163 composés de "
            "la campagne \\texttt{extension}. $^*$ noir = composé de départ expérimental "
            "stable ; \\textcolor{red}{$^*$} rouge = expérimental métastable."
        ),
        "exo_caption": (
            "Réactions de décomposition en éléments (cas 1) classées \\textbf{exobondic} "
            "($\\Delta$ICOHP $< 0$, produits $-$ réactifs) parmi les 163 composés de la "
            "campagne \\texttt{extension}. $^*$ noir = composé de départ expérimental "
            "stable ; \\textcolor{red}{$^*$} rouge = expérimental métastable."
        ),
        "reactions_header": ["Réaction (balancée)", "mp-id / COD", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)"],
        "endo": "endobondic", "exo": "exobondic",
    },
    "en": {
        "elements_caption": (
            "Elements and allotropes with no matched polymorph in the "
            "\\texttt{extension} campaign: references used to "
            "decompose multi-element compounds. Allotropes with one or more "
            "other allotropes computed in this dataset (e.g. carbon) are "
            "grouped in Table~\\ref{tab:ext-polymorphs} instead. Black $^*$ = "
            "stable experimental reference (thermodynamic phase); red "
            "\\textcolor{red}{$^*$} = metastable experimental reference; "
            "unmarked = theoretical. No $\\Delta$ICOHP column: an element is "
            "not itself a decomposition reaction."
        ),
        "elements_header": ["Formula", "mp-id / COD", "Space group", "$E_{hull}$ (eV/at)"],
        "compounds_caption": (
            "Target compounds with no matched polymorph in the "
            "\\texttt{extension} campaign. Formulas with "
            "several entries in this dataset are grouped in "
            "Table~\\ref{tab:ext-polymorphs} instead. Black $^*$ = stable "
            "experimental compound (thermodynamic phase); red "
            "\\textcolor{red}{$^*$} = metastable experimental compound; "
            "unmarked = theoretical. $\\Delta$ICOHP per atom, products $-$ "
            "reactants convention (\\texttt{reaction\\_analysis}); ``--'' "
            "where no case-1 reaction could be computed (missing elemental "
            "reference)."
        ),
        "compounds_header": ["Formula", "mp-id / COD", "Space group", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)", "Label"],
        "polymorphs_caption": (
            "Elements (allotropes) and compounds (polymorphs) with "
            "\\textbf{several entries sharing a formula} in the "
            "\\texttt{extension} campaign (39 groups, 84 entries) --- "
            "grouped here rather than scattered across the elements/"
            "compounds tables, since the natural comparison for a polymorph "
            "group is polymorph-to-polymorph, not the decomposition-to-"
            "elements reaction (hence no $\\Delta$ICOHP column). A rule "
            "separates each group; within a group, sorted by increasing "
            "$E_{hull}$. Black $^*$ = stable experimental; red "
            "\\textcolor{red}{$^*$} = metastable experimental; unmarked = "
            "theoretical."
        ),
        "endo_caption": (
            "Decomposition-into-elements (case 1) reactions labeled \\textbf{endobondic} "
            "($\\Delta$ICOHP $\\geq 0$, products $-$ reactants) among the 163 compounds "
            "of the \\texttt{extension} campaign. Black $^*$ = stable experimental "
            "starting compound; red \\textcolor{red}{$^*$} = metastable experimental."
        ),
        "exo_caption": (
            "Decomposition-into-elements (case 1) reactions labeled \\textbf{exobondic} "
            "($\\Delta$ICOHP $< 0$, products $-$ reactants) among the 163 compounds of "
            "the \\texttt{extension} campaign. Black $^*$ = stable experimental starting "
            "compound; red \\textcolor{red}{$^*$} = metastable experimental."
        ),
        "reactions_header": ["Reaction (balanced)", "mp-id / COD", "$E_{hull}$ (eV/at)", "$\\Delta$ICOHP/at. (eV)"],
        "endo": "endobondic", "exo": "exobondic",
    },
}


def load_all_rows():
    element_rows, compound_rows = [], []
    for d in sorted(STRUCTURES_ROOT.iterdir()):
        meta_path = d / "mp_metadata.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("family") != "extension":
            continue
        if meta.get("batch") in EXCLUDED_BATCHES:
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

    formula_counts = Counter(r["formula"] for r in element_rows + compound_rows)

    ri = pd.read_csv(REPO_ROOT / "analysis" / "reaction_icohp_case1.csv")[["compound_id", "reaction_string"]]
    ra = pd.read_csv(REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv")[["compound_id", "delta_per_atom_eV"]]
    reaction_lookup = ri.merge(ra, on="compound_id", how="inner").set_index("compound_id").to_dict("index")

    for row in compound_rows:
        r = reaction_lookup.get(row["compound_id"])
        row["reaction_string"] = r["reaction_string"] if r else None
        row["delta_per_atom_eV"] = r["delta_per_atom_eV"] if r else None

    return element_rows, compound_rows, formula_counts


def _longtable(caption: str, header: list[str], colspec: str, body_lines: list[str], label: str | None = None) -> str:
    lines = [
        f"\\begin{{longtable}}{{@{{}}{colspec}@{{}}}}",
        f"\\caption{{{caption}}}" + (f"\\label{{{label}}}" if label else "") + "\\\\",
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


def write_elements(lang: str, rows: list[dict], formula_counts: Counter) -> None:
    L = TEXT[lang]
    singles = [r for r in rows if formula_counts[r["formula"]] == 1]
    body = []
    for r in sorted(singles, key=lambda r: r["formula"]):
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        body.append(f"{_star(r['formula'], r['theoretical'], r['eah'])} & {_esc(r['id_label'])} & {_esc(r['spacegroup'])} & {eah} \\\\")
    out = _longtable(L["elements_caption"], L["elements_header"], "llll", body)
    (HERE / f"appendix_extension_elements_{lang}.tex").write_text(out)


def write_compounds(lang: str, rows: list[dict], formula_counts: Counter) -> None:
    L = TEXT[lang]
    singles = [r for r in rows if formula_counts[r["formula"]] == 1]
    body = []
    for r in sorted(singles, key=lambda r: r["formula"]):
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        has_reaction = r["delta_per_atom_eV"] is not None
        delta = f"{r['delta_per_atom_eV']:.4f}" if has_reaction else "--"
        label = (L["endo"] if r["delta_per_atom_eV"] >= 0 else L["exo"]) if has_reaction else "--"
        body.append(
            f"{_star(r['formula'], r['theoretical'], r['eah'])} & {_esc(r['id_label'])} & {_esc(r['spacegroup'])} "
            f"& {eah} & {delta} & {label} \\\\"
        )
    out = _longtable(L["compounds_caption"], L["compounds_header"], "lllrrl", body)
    (HERE / f"appendix_extension_compounds_{lang}.tex").write_text(out)


def write_polymorphs(lang: str, element_rows: list[dict], compound_rows: list[dict], formula_counts: Counter) -> None:
    L = TEXT[lang]
    grouped = [r for r in element_rows + compound_rows if formula_counts[r["formula"]] >= 2]
    grouped.sort(key=lambda r: (r["formula"], r["eah"] if r["eah"] is not None else 0.0))

    body = []
    prev_formula = None
    for r in grouped:
        if prev_formula is not None and r["formula"] != prev_formula:
            body.append("\\midrule")
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        body.append(f"{_star(r['formula'], r['theoretical'], r['eah'])} & {_esc(r['id_label'])} & {_esc(r['spacegroup'])} & {eah} \\\\")
        prev_formula = r["formula"]

    out = _longtable(L["polymorphs_caption"], L["elements_header"], "llll", body, label="tab:ext-polymorphs")
    (HERE / f"appendix_extension_polymorphs_{lang}.tex").write_text(out)


def _starred_reaction_string(r: dict) -> str:
    rxn = _subscript(r["reaction_string"]).replace("->", "$\\rightarrow$")
    star = _star_markup(r["theoretical"], r["eah"])
    if star:
        target_sub = _subscript(r["formula"])
        rxn = re.sub(re.escape(target_sub), target_sub + star, rxn, count=1)
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
    element_rows, compound_rows, formula_counts = load_all_rows()
    n_grouped = sum(1 for r in element_rows + compound_rows if formula_counts[r["formula"]] >= 2)
    print(f"{len(element_rows)} element/allotrope rows, {len(compound_rows)} compound/polymorph rows")
    print(f"{n_grouped} rows belong to a polymorph/allotrope group ({sum(1 for c in formula_counts.values() if c >= 2)} groups)")
    n_with_reaction = sum(1 for r in compound_rows if r["delta_per_atom_eV"] is not None)
    n_endo = sum(1 for r in compound_rows if r["delta_per_atom_eV"] is not None and r["delta_per_atom_eV"] >= 0)
    n_exo = sum(1 for r in compound_rows if r["delta_per_atom_eV"] is not None and r["delta_per_atom_eV"] < 0)
    print(f"with a case-1 reaction: {n_with_reaction} (endobondic={n_endo}, exobondic={n_exo})")

    for lang in ("fr", "en"):
        write_elements(lang, element_rows, formula_counts)
        write_compounds(lang, compound_rows, formula_counts)
        write_polymorphs(lang, element_rows, compound_rows, formula_counts)
        write_reactions(lang, compound_rows, endobondic=True)
        write_reactions(lang, compound_rows, endobondic=False)

    for stale in list(HERE.glob("appendix_extension4_*.tex")):
        stale.unlink()
        print(f"removed superseded {stale.name}")

    print("wrote appendix_extension_{elements,compounds,polymorphs,endobondic,exobondic}_{fr,en}.tex")


if __name__ == "__main__":
    main()
