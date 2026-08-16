"""Generate the master compound/element list for the SI's first section:
every one of the 392 compounds with a usable ICOBI-based classification
(analysis/icohp_icobi_bondtype.csv, see report Sec. 5 for the
classification method itself), one longtable per bonding category
(metallic / ionic / covalent / not classified) instead of a bond-type
column -- each table's caption states the ICOBI range observed in that
category. No percolation weight column here (kept only in the
Percolation section's own duplicated table, see gen_appendix.py).

Columns: Formula (starred), mp-id, space group, ICOBI (mean/bond), ICOHP
(mean/bond, eV), ICOHP antibonding (normalized, dE=1.0), ICOBI
antibonding (normalized, dE=1.0).

Writes, per language (fr/en): appendix_master_list_{category}_{lang}.tex
for category in (metallic, ionic, covalent, unclassified).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent

BONDTYPE_CSV = REPO_ROOT / "analysis" / "icohp_icobi_bondtype.csv"
ICOBI_ANTIBOND_CSV = REPO_ROOT / "analysis" / "icobi_antibonding_all.csv"
THRESHOLD_JSON = REPO_ROOT / "analysis" / "icohp_icobi_bondtype_threshold.json"

_SUBSCRIPT_RE = re.compile(r"(?<=[A-Za-z\)])(\d+)")
_EPS_HULL = 1e-6

CATEGORIES = ["metallic", "ionic", "covalent", "mixed", "not_classified"]
CATEGORY_KEY = {"metallic": "metallic", "ionic": "ionic", "covalent": "covalent", "mixed": "mixed", "not_classified": "unclassified"}


def _subscript(text: str) -> str:
    return _SUBSCRIPT_RE.sub(r"$_{\1}$", text)


def _esc(s: str) -> str:
    return str(s).replace("_", "\\_").replace("#", "\\#")


def _star_markup(theoretical: bool, energy_above_hull) -> str:
    if theoretical:
        return ""
    if energy_above_hull is not None and pd.notna(energy_above_hull) and energy_above_hull <= _EPS_HULL:
        return "$^*$"
    return "\\textcolor{red}{$^*$}"


def _star(formula: str, theoretical: bool, energy_above_hull) -> str:
    return _subscript(_esc(formula)) + _star_markup(theoretical, energy_above_hull)


TEXT = {
    "fr": {
        "header": ["Formule", "mp-id", "Groupe d'espace", "ICOBI", "ICOHP (eV)", "ICOHP antiliant", "ICOBI antiliant"],
        "caption": {
            "metallic": "Composés et éléments métalliques (\\texttt{is\\_metal=True}, ICOBI observé de __LO__ à __HI__) --- formule, mp-id, groupe d'espace, ICOBI de la liaison dominante (premier voisinage), ICOHP de la liaison dominante, populations antiliantes ICOHP et ICOBI normalisées ($\\Delta E=1{,}0$~eV).",
            "ionic": "Composés ioniques (non métalliques, ICOBI __THRESH__, plage observée __LO__ à __HI__) --- mêmes colonnes.",
            "covalent": "Composés covalents (non métalliques, ICOBI __THRESH2__, plage observée __LO__ à __HI__) --- mêmes colonnes.",
            "mixed": "Composés à liaison mixte, type Zintl (une paire homoatomique fortement covalente coexiste avec une paire hétéroatomique nettement plus faible, p.~ex. l'azoture N$_3^-$ lié ioniquement au cation, plage observée __LO__ à __HI__) --- mêmes colonnes.",
            "not_classified": "Composés non classifiés (données ICOBI indisponibles) --- mêmes colonnes, ICOBI/populations antiliantes en \\og{}--\\fg{} le cas échéant.",
        },
        "legend": (
            "$^*$ noir = composé expérimental stable (phase thermodynamique) ; "
            "\\textcolor{red}{$^*$} rouge = composé expérimental métastable ; "
            "non marqué = théorique. ICOBI/ICOHP : valeur de la paire "
            "d'espèces dominante au premier voisinage (pas une moyenne sur "
            "tous les voisins LOBSTER, voir \\S2.4)."
        ),
        "na": "--",
    },
    "en": {
        "header": ["Formula", "mp-id", "Space group", "ICOBI", "ICOHP (eV)", "ICOHP antibonding", "ICOBI antibonding"],
        "caption": {
            "metallic": "Metallic compounds and elements (\\texttt{is\\_metal=True}, observed ICOBI __LO__ to __HI__) --- formula, mp-id, space group, dominant (first-shell) bond ICOBI, dominant bond ICOHP, normalized ICOHP and ICOBI antibonding populations ($\\Delta E=1.0$~eV).",
            "ionic": "Ionic compounds (non-metallic, ICOBI __THRESH__, observed range __LO__ to __HI__) --- same columns.",
            "covalent": "Covalent compounds (non-metallic, ICOBI __THRESH2__, observed range __LO__ to __HI__) --- same columns.",
            "mixed": "Mixed-bonding (Zintl-type) compounds (a strongly covalent homoatomic pair coexists with a markedly weaker heteroatomic pair, e.g. the azide N$_3^-$ ion ionically bound to its cation; observed range __LO__ to __HI__) --- same columns.",
            "not_classified": "Unclassified compounds (no ICOBI data available) --- same columns, ICOBI/antibonding shown as ``--'' where missing.",
        },
        "legend": (
            "Black $^*$ = stable experimental compound (thermodynamic "
            "phase); red \\textcolor{red}{$^*$} = metastable experimental "
            "compound; unmarked = theoretical. ICOBI/ICOHP: value of the "
            "dominant first-shell species pair (not a mean over every "
            "LOBSTER-reported neighbor, see \\S2.4)."
        ),
        "na": "--",
    },
}


def _fmt(x, digits=4) -> str:
    if x is None or pd.isna(x):
        return "--"
    return f"{x:.{digits}f}"


def load_data() -> pd.DataFrame:
    df = pd.read_csv(BONDTYPE_CSV)
    antibond = pd.read_csv(ICOBI_ANTIBOND_CSV)[["compound_id", "icobi_antibond_w_normalized"]]
    df = df.merge(antibond, on="compound_id", how="left")
    return df


def _longtable(caption: str, header: list[str], body: list[str], legend: str) -> str:
    ncol = len(header)
    lines = [
        "\\begin{longtable}{@{}lllrrrr@{}}",
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
        "\\bottomrule",
        f"\\multicolumn{{{ncol}}}{{@{{}}p{{0.95\\linewidth}}}}{{\\footnotesize {legend}}} \\\\",
        "\\endlastfoot",
        *body,
        "\\end{longtable}",
    ]
    return "\n".join(lines) + "\n"


def write_category(lang: str, df: pd.DataFrame, category: str, threshold: float) -> None:
    L = TEXT[lang]
    sub = df[df["icobi_label"] == category]
    sub = sub.sort_values("formula")

    lo = _fmt(sub["icobi_mean"].min()) if not sub.empty else "--"
    hi = _fmt(sub["icobi_mean"].max()) if not sub.empty else "--"
    thresh_str = f"{threshold:.4f}".replace(".", "{,}") if lang == "fr" else f"{threshold:.4f}"
    caption = (
        L["caption"][category]
        .replace("__LO__", lo).replace("__HI__", hi)
        .replace("__THRESH2__", f"$\\geq {thresh_str}$")
        .replace("__THRESH__", f"$<{thresh_str}$")
    )

    body = []
    for _, r in sub.iterrows():
        star = _star(r["formula"], bool(r["theoretical"]), r.get("energy_above_hull_eV_at"))
        mp_id = r["mp_id"] if pd.notna(r["mp_id"]) else "--"
        sg = _esc(r["spacegroup_symbol"]) if pd.notna(r["spacegroup_symbol"]) else "--"
        body.append(
            f"{star} & {_esc(mp_id)} & {sg} & {_fmt(r['icobi_mean'])} & {_fmt(r['icohp_mean'])} "
            f"& {_fmt(r['antibond_w_normalized'])} & {_fmt(r['icobi_antibond_w_normalized'])} \\\\"
        )

    out = _longtable(caption, L["header"], body, L["legend"])
    key = CATEGORY_KEY[category]
    (HERE / f"appendix_master_list_{key}_{lang}.tex").write_text(out)
    print(f"{lang} {key}: {len(sub)} rows, ICOBI [{lo}, {hi}]")


def main() -> None:
    df = load_data()
    print(f"{len(df)} total rows in {BONDTYPE_CSV.name}")
    threshold = json.loads(THRESHOLD_JSON.read_text())["covalent_threshold_icobi_primary_mean"]
    for lang in ("fr", "en"):
        for category in CATEGORIES:
            write_category(lang, df, category, threshold)
    print("wrote appendix_master_list_{metallic,ionic,covalent,mixed,unclassified}_{fr,en}.tex")


if __name__ == "__main__":
    main()
