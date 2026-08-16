"""Generate LaTeX appendix fragments for the extension4 campaign (89
alkali/alkaline-earth binary compounds + their 16 constituent elemental
references), and for its reaction_icohp/reaction_analysis case-1
decomposition reactions -- requested explicitly by the user (list every
compound and element computed, detail the reactions considered).

Writes, per language (fr/en):
  appendix_extension4_compounds_{lang}.tex   -- 89 targets + 16 elements
  appendix_extension4_reactions_{lang}.tex   -- 89 decomposition reactions,
    reaction string from reaction_icohp_case1.csv (reactant-minus-products
    convention, human-readable), delta_per_atom_eV and bonding_label from
    reaction_analysis_case1_full.csv's products-minus-reactants convention
    (same sign classify.py uses -- endobondic if delta_per_atom_eV>=0).

Run manually after any dataset change; not part of the automated analysis
pipeline, same convention as gen_appendix.py / gen_appendix_dft.py.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"

sys.path.insert(0, str(REPO_ROOT / "analysis"))
import compute_reaction_icohp_case1 as ri_case1  # noqa: E402  (reused unmodified, for ELEMENT_REFERENCE)

EXTENSION4_ELEMENTS = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "N", "O", "F", "P", "S", "Cl"]

LABEL = {
    "fr": {
        "compounds_caption": (
            "Composés cibles de la campagne \\texttt{extension4} (89) et références "
            "élémentaires utilisées pour leur décomposition (16) -- \\texttt{mp\\_dataset/"
            "download\\_extension4.py} / \\texttt{ELEMENT\\_REFERENCE}."
        ),
        "compounds_header": ["Formule", "mp-id", "Type", "Groupe d'espace", "$E_{hull}$ (eV/at)", "Rôle"],
        "target": "cible", "element": "élément",
        "exp": "exp.", "theo": "théo.",
        "reactions_caption": (
            "Réactions de décomposition en éléments (case 1) calculées pour les 89 composés "
            "\\texttt{extension4}, avec $\\Delta$ICOHP par atome (convention produits $-$ "
            "réactifs, \\texttt{reaction\\_analysis}) et l'étiquette endobondic/exobondic "
            "correspondante."
        ),
        "reactions_header": ["Réaction (balancée)", "$\\Delta$ICOHP/atome (eV)", "Étiquette"],
        "endo": "endobondic", "exo": "exobondic",
    },
    "en": {
        "compounds_caption": (
            "Target compounds of the \\texttt{extension4} campaign (89) and elemental "
            "references used for their decomposition (16) -- \\texttt{mp\\_dataset/"
            "download\\_extension4.py} / \\texttt{ELEMENT\\_REFERENCE}."
        ),
        "compounds_header": ["Formula", "mp-id", "Type", "Space group", "$E_{hull}$ (eV/at)", "Role"],
        "target": "target", "element": "element",
        "exp": "exp.", "theo": "theo.",
        "reactions_caption": (
            "Decomposition-into-elements (case 1) reactions computed for the 89 "
            "\\texttt{extension4} compounds, with $\\Delta$ICOHP per atom (products $-$ "
            "reactants convention, \\texttt{reaction\\_analysis}) and the corresponding "
            "endobondic/exobondic label."
        ),
        "reactions_header": ["Reaction (balanced)", "$\\Delta$ICOHP/atom (eV)", "Label"],
        "endo": "endobondic", "exo": "exobondic",
    },
}


def _esc(s: str) -> str:
    return s.replace("_", "\\_").replace("#", "\\#")


def load_compound_rows():
    picks = json.loads((REPO_ROOT / "mp_dataset" / "extension4_candidates.json").read_text())
    kind_suffix = {"exp_polymorph": "exp", "theo_far_hull": "theo"}
    rows = []
    for p in picks:
        cid = f"extension_{p['formula']}_{kind_suffix[p['kind']]}_{p['material_id']}"
        meta = json.loads((STRUCTURES_ROOT / cid / "mp_metadata.json").read_text())
        rows.append({
            "formula": meta["formula"], "mp_id": meta["mp_id"],
            "theoretical": meta["theoretical"], "spacegroup": meta.get("spacegroup") or "--",
            "eah": meta.get("energy_above_hull_eV_per_atom"), "role": "target",
        })
    for el in EXTENSION4_ELEMENTS:
        cid = ri_case1.ELEMENT_REFERENCE[el]
        meta = json.loads((STRUCTURES_ROOT / cid / "mp_metadata.json").read_text())
        rows.append({
            "formula": meta["formula"], "mp_id": meta["mp_id"],
            "theoretical": meta["theoretical"], "spacegroup": meta.get("spacegroup") or "--",
            "eah": meta.get("energy_above_hull_eV_per_atom"), "role": "element",
        })
    return rows


def write_compounds_tex(lang: str, rows: list[dict]) -> None:
    L = LABEL[lang]
    lines = [
        "\\begin{longtable}{@{}lllllr@{}}",
        f"\\caption{{{L['compounds_caption']}}}\\\\",
        "\\toprule",
        " & ".join(L["compounds_header"]) + " \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        " & ".join(L["compounds_header"]) + " \\\\",
        "\\midrule",
        "\\endhead",
        "\\bottomrule",
        "\\endfoot",
    ]
    for r in sorted(rows, key=lambda r: (r["role"], r["formula"])):
        typ = L["theo"] if r["theoretical"] else L["exp"]
        role = L["target"] if r["role"] == "target" else L["element"]
        eah = f"{r['eah']:.4f}" if r["eah"] is not None else "--"
        lines.append(
            f"{_esc(r['formula'])} & {_esc(r['mp_id'])} & {typ} & {_esc(r['spacegroup'])} & {eah} & {role} \\\\"
        )
    lines.append("\\end{longtable}")
    (HERE / f"appendix_extension4_compounds_{lang}.tex").write_text("\n".join(lines) + "\n")


def load_reaction_rows():
    ri = pd.read_csv(REPO_ROOT / "analysis" / "reaction_icohp_case1.csv")
    ra = pd.read_csv(REPO_ROOT / "analysis" / "reaction_analysis_case1_full.csv")
    picks = json.loads((REPO_ROOT / "mp_dataset" / "extension4_candidates.json").read_text())
    kind_suffix = {"exp_polymorph": "exp", "theo_far_hull": "theo"}
    ext4_ids = {f"extension_{p['formula']}_{kind_suffix[p['kind']]}_{p['material_id']}" for p in picks}

    ri = ri[ri["compound_id"].isin(ext4_ids)][["compound_id", "reaction_string"]]
    ra = ra[ra["compound_id"].isin(ext4_ids)][["compound_id", "delta_per_atom_eV"]]
    merged = ri.merge(ra, on="compound_id", how="inner")
    return merged.to_dict("records")


def write_reactions_tex(lang: str, rows: list[dict]) -> None:
    L = LABEL[lang]
    lines = [
        "\\begin{longtable}{@{}lrl@{}}",
        f"\\caption{{{L['reactions_caption']}}}\\\\",
        "\\toprule",
        " & ".join(L["reactions_header"]) + " \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        " & ".join(L["reactions_header"]) + " \\\\",
        "\\midrule",
        "\\endhead",
        "\\bottomrule",
        "\\endfoot",
    ]
    for r in sorted(rows, key=lambda r: r["compound_id"]):
        label = L["endo"] if r["delta_per_atom_eV"] >= 0 else L["exo"]
        rxn = r["reaction_string"].replace("->", "$\\rightarrow$")
        lines.append(f"{rxn} & {r['delta_per_atom_eV']:.4f} & {label} \\\\")
    lines.append("\\end{longtable}")
    (HERE / f"appendix_extension4_reactions_{lang}.tex").write_text("\n".join(lines) + "\n")


def main():
    compound_rows = load_compound_rows()
    print(f"{len(compound_rows)} compound rows (89 targets + 16 elements expected)")
    reaction_rows = load_reaction_rows()
    print(f"{len(reaction_rows)} reaction rows (89 expected)")
    for lang in ("fr", "en"):
        write_compounds_tex(lang, compound_rows)
        write_reactions_tex(lang, reaction_rows)
    print("wrote appendix_extension4_{compounds,reactions}_{fr,en}.tex")


if __name__ == "__main__":
    main()
