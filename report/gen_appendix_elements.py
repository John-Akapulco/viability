"""Generate the full elemental-reference-structure SI table: all 62
entries of analysis/compute_reaction_icohp_case1.ELEMENT_REFERENCE (the
"products" side of every case-1 decomposition reaction), from
analysis/element_reference_energetics.csv
(analysis/compute_element_reference_energetics.py).

Columns: Element, mp-id, Space group, MP E/atom (eV), Campaign E/atom
(eV, this project's own VASP OSZICAR E0), Delta (meV/atom), ICOHP (eV),
ICOHP antibonding (normalized).

The Delta column is large (up to tens of eV/atom) for elements using a
PAW potential with more valence electrons treated explicitly (e.g. a
"_d"/"_pv"/"_sv" semi-core variant) -- this is an EXPECTED, well-known
consequence of comparing raw absolute VASP total energies across
different pseudopotential valence-electron conventions, not a DFT-
pipeline error: absolute total energies are only meaningful within one
consistent potential set. The caption says so explicitly.

Writes appendix_elements_{fr,en}.tex.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
REPO_ROOT = HERE.parent
CSV_PATH = REPO_ROOT / "analysis" / "element_reference_energetics.csv"


def _fmt(x, digits=4) -> str:
    if x is None or pd.isna(x):
        return "--"
    return f"{x:.{digits}f}"


def _fmt_fr(x, digits=4) -> str:
    s = _fmt(x, digits)
    return s if s == "--" else s.replace(".", "{,}").replace("-", "$-$")


TEXT = {
    "fr": {
        "header": ["Élém.", "mp-id", "GE", "$E_\\text{MP}$ (eV)",
                   "$E_\\text{camp}$ (eV)", "$\\Delta$ (meV)", "ICOHP", "ICOHP antil."],
        "caption": (
            "Les 62 structures de référence élémentaires utilisées comme "
            "produits de chaque réaction case-1 (décomposition en "
            "éléments ; GE = groupe d'espace). "
            "$\\Delta = E_\\text{campagne} - E_\\text{MP}$ par "
            "atome : les grands écarts (jusqu'à plusieurs dizaines "
            "d'eV/atome pour les métaux lourds) reflètent une différence "
            "de convention de pseudopotentiel PAW (nombre d'électrons de "
            "semi-c\\oe{}ur trait\\'es explicitement comme \\'electrons de "
            "valence, variantes \\texttt{\\_d}/\\texttt{\\_pv}/\\texttt{\\_sv}) "
            "et NON une erreur du pipeline DFT : les \\'energies totales "
            "absolues ne sont comparables qu'au sein d'un m\\^eme jeu de "
            "pseudopotentiels. N n'a pas "
            "d'entr\\'ee MP (bo\\^ite-dim\\`ere construite \\`a la main, "
            "\\S\\ref{sec:appendix-master-list})."
        ),
        "na": "--",
    },
    "en": {
        "header": ["Elem.", "mp-id", "SG", "$E_\\text{MP}$ (eV)",
                   "$E_\\text{camp}$ (eV)", "$\\Delta$ (meV)", "ICOHP", "ICOHP antib."],
        "caption": (
            "The 62 elemental reference structures used as the products "
            "side of every case-1 (decomposition-to-elements) reaction "
            "(SG = space group). "
            "$\\Delta = E_\\text{campaign} - E_\\text{MP}$ per atom: the "
            "large offsets (up to several tens of eV/atom for heavy "
            "metals) reflect a PAW pseudopotential valence-electron "
            "convention difference (semi-core electrons explicitly "
            "included, \\texttt{\\_d}/\\texttt{\\_pv}/\\texttt{\\_sv} "
            "variants), NOT a DFT-pipeline error: absolute total "
            "energies are only meaningful within one consistent "
            "pseudopotential set. N has no MP entry (hand-built "
            "dimer box, \\S\\ref{sec:appendix-master-list})."
        ),
        "na": "--",
    },
}


def _longtable(caption: str, header: list[str], body: list[str]) -> str:
    ncol = len(header)
    lines = [
        "\\begin{longtable}{@{}ll l@{\\hspace{4pt}} rrrrr@{}}",
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
        "\\endlastfoot",
        *body,
        "\\end{longtable}",
    ]
    return "\n".join(lines) + "\n"


def write(lang: str, df: pd.DataFrame) -> None:
    L = TEXT[lang]
    fmt = _fmt_fr if lang == "fr" else _fmt
    body = []
    for _, r in df.iterrows():
        delta_meV = r["mp_vs_campaign_diff_meV_per_atom"]
        body.append(
            f"{r['element']} & {r['mp_id'] if pd.notna(r['mp_id']) else L['na']} & "
            f"{r['spacegroup'] if pd.notna(r['spacegroup']) else L['na']} & "
            f"{fmt(r['mp_energy_per_atom_eV'], 2)} & {fmt(r['campaign_energy_per_atom_eV'], 2)} & "
            f"{fmt(delta_meV, 0)} & {fmt(r['icohp_mean'], 3)} & {fmt(r['antibond_w_normalized'], 3)} \\\\"
        )
    out = _longtable(L["caption"], L["header"], body)
    (HERE / f"appendix_elements_{lang}.tex").write_text(out)
    print(f"{lang}: {len(df)} elements written to appendix_elements_{lang}.tex")


def main() -> None:
    df = pd.read_csv(CSV_PATH).sort_values("element").reset_index(drop=True)
    for lang in ("fr", "en"):
        write(lang, df)


if __name__ == "__main__":
    main()
