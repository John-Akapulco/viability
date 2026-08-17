"""Generate the LaTeX appendix (longtable) listing every compound in the
campaign-2 dataset, grouped by family, with mp-id, PBE hull distance, and
the percolation-weight metric. Writes appendix_compounds_fr.tex and
appendix_compounds_en.tex (included via \\input{} from the two report
.tex files). Run manually after any dataset change; not part of the
automated analysis pipeline.
"""

import re
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
CSV_PATH = HERE.parent / "analysis" / "percolation_vs_hull.csv"

_SUBSCRIPT_RE = re.compile(r"(?<=[A-Za-z\)])(\d+)")


def _subscript(text: str) -> str:
    return _SUBSCRIPT_RE.sub(r"$_{\1}$", text)

FAMILY_ORDER = ["exp_stable", "exp_metastable", "theo_metastable"]

FAMILY_LABEL = {
    "fr": {
        "exp_stable": "Exp\\'erimental, stable (hull)",
        "exp_metastable": "Exp\\'erimental, m\\'etastable",
        "theo_metastable": "Th\\'eorique, m\\'etastable",
    },
    "en": {
        "exp_stable": "Experimental, stable (hull)",
        "exp_metastable": "Experimental, metastable",
        "theo_metastable": "Theoretical, metastable",
    },
}

BOND_LABEL_FR = {"ionic": "ionique", "covalent": "covalent", "metallic": "m\\'etallique", None: "--"}
BOND_LABEL_EN = {"ionic": "ionic", "covalent": "covalent", "metallic": "metallic", None: "--"}


def fmt_e(x: float) -> str:
    return f"{x:.4f}"


def fmt_w(x) -> str:
    if pd.isna(x):
        return "--"
    return f"{x:.2e}"


def build(lang: str) -> str:
    df = pd.read_csv(CSV_PATH)
    df = df.sort_values(["family", "energy_above_hull_eV_at"])
    bond_label = BOND_LABEL_FR if lang == "fr" else BOND_LABEL_EN

    header = {
        "fr": ("Formule", "mp-id", "Liaison", r"$E_\text{hull}$ (eV/at)", "Poids percolation (eV)"),
        "en": ("Formula", "mp-id", "Bonding", r"$E_\text{hull}$ (eV/at)", "Percolation weight (eV)"),
    }[lang]
    caption = {
        "fr": "Liste compl\\`ete des 186 compos\\'es, class\\'es par cat\\'egorie, avec identifiant Materials Project, distance au hull (PBE), et poids de percolation ICOHP minimal.",
        "en": "Complete list of the 186 compounds, grouped by category, with their Materials Project identifier, PBE hull distance, and minimum ICOHP percolation weight.",
    }[lang]
    label = "tab:appendix-compounds"
    continued = {"fr": "(suite)", "en": "(continued)"}[lang]

    lines = []
    lines.append(r"\begin{longtable}{@{}l l l r r@{}}")
    lines.append(rf"\caption{{{caption}}}\label{{{label}}}\\")
    lines.append(r"\toprule")
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\midrule")
    lines.append(r"\endfirsthead")
    lines.append(rf"\multicolumn{{5}}{{@{{}}l}}{{\tablename\ \thetable\ {continued}}}\\")
    lines.append(r"\toprule")
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\midrule")
    lines.append(r"\endhead")
    lines.append(r"\bottomrule")
    lines.append(r"\endfoot")

    for family in FAMILY_ORDER:
        sub = df[df["family"] == family]
        flabel = FAMILY_LABEL[lang][family]
        lines.append(rf"\multicolumn{{5}}{{@{{}}l}}{{\textbf{{{flabel}}} ($n={len(sub)}$)}} \\")
        lines.append(r"\midrule")
        for _, row in sub.iterrows():
            formula = _subscript(str(row["formula"]))
            mp_id = str(row["mp_id"])
            bt = row["bond_type"] if pd.notna(row["bond_type"]) else None
            bt_str = bond_label.get(bt, bt)
            ehull = fmt_e(row["energy_above_hull_eV_at"])
            w = fmt_w(row["icohp_percolation_weight_min"])
            lines.append(f"{formula} & {mp_id} & {bt_str} & {ehull} & {w} \\\\")
        lines.append(r"\midrule")

    lines.append(r"\end{longtable}")
    return "\n".join(lines)


def main():
    (HERE / "appendix_compounds_fr.tex").write_text(build("fr"))
    (HERE / "appendix_compounds_en.tex").write_text(build("en"))
    print("wrote appendix_compounds_fr.tex and appendix_compounds_en.tex")


if __name__ == "__main__":
    main()
