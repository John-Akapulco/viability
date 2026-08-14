"""Generate the LaTeX appendix fragments describing the VASP DFT
calculation parameters: fixed INCAR settings (identical across all 186
compounds), the PAW-PBE potential chosen per element, and the per-compound
k-points grid + NBANDS (which do vary by compound). Reads directly from
mp_dataset/structures/*/{INCAR,KPOINTS,POTCAR} on disk (POTCAR is not
tracked in git -- VASP license -- but is present locally).

Writes appendix_dft_common_{fr,en}.tex, appendix_dft_potcar_{fr,en}.tex,
appendix_dft_kpoints_{fr,en}.tex. Run manually after any dataset change;
not part of the automated analysis pipeline.
"""

import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
STRUCTURES_ROOT = HERE.parent / "mp_dataset" / "structures"
CSV_PATH = HERE.parent / "analysis" / "percolation_vs_hull.csv"

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

# INCAR keys that are identical across the whole dataset (everything except
# NBANDS, which is derived per compound from the LOBSTER basis size).
COMMON_KEYS = ["PREC", "ENCUT", "EDIFF", "IBRION", "NSW", "ISMEAR", "SIGMA",
                "LASPH", "ISYM", "LWAVE", "LCHARG", "NPAR", "KPAR"]

DESCRIPTION = {
    "fr": {
        "PREC": "Pr\\'ecision globale (grille de charge, projecteurs)",
        "ENCUT": "\\'Energie de coupure des ondes planes (eV)",
        "EDIFF": "Crit\\`ere de convergence \\'electronique (eV)",
        "IBRION": "Type de dynamique ionique ($-1$ = statique, aucune relaxation)",
        "NSW": "Nombre de pas ioniques (0 = calcul statique sur la structure MP)",
        "ISMEAR": "Type d'\\'elargissement ($0$ = gaussien)",
        "SIGMA": "Largeur d'\\'elargissement (eV)",
        "LASPH": "Corrections non sph\\'eriques au potentiel PAW",
        "ISYM": "Sym\\'etrie ($-1$ = d\\'esactiv\\'ee, requis par LOBSTER)",
        "LWAVE": "\\'Ecriture du WAVECAR (requis par LOBSTER)",
        "LCHARG": "\\'Ecriture du CHGCAR",
        "NPAR": "Parall\\'elisation sur les bandes",
        "KPAR": "Parall\\'elisation sur les points k",
    },
    "en": {
        "PREC": "Overall precision (charge grid, projectors)",
        "ENCUT": "Plane-wave cutoff energy (eV)",
        "EDIFF": "Electronic convergence criterion (eV)",
        "IBRION": "Ionic dynamics type ($-1$ = static, no relaxation)",
        "NSW": "Number of ionic steps (0 = static run on the MP structure)",
        "ISMEAR": "Smearing type ($0$ = Gaussian)",
        "SIGMA": "Smearing width (eV)",
        "LASPH": "Non-spherical corrections to the PAW potential",
        "ISYM": "Symmetry ($-1$ = disabled, required by LOBSTER)",
        "LWAVE": "Write WAVECAR (required by LOBSTER)",
        "LCHARG": "Write CHGCAR",
        "NPAR": "Band-level parallelization",
        "KPAR": "K-point-level parallelization",
    },
}


def parse_incar(path: Path) -> dict:
    values = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        values[key.strip()] = val.strip()
    return values


def parse_kpoints(path: Path) -> tuple[str, str]:
    lines = path.read_text().splitlines()
    mesh_type = lines[2].strip()
    grid = lines[3].strip()
    return mesh_type, grid


def gather_potcar_titels() -> dict[str, str]:
    """element -> TITEL string, collected once per element (identical across
    every compound that uses it, since the mapping is fixed per element)."""
    out: dict[str, str] = {}
    for compound_dir in sorted(STRUCTURES_ROOT.iterdir()):
        potcar = compound_dir / "POTCAR"
        if not potcar.exists():
            continue
        text = potcar.read_text(errors="replace")
        titels = re.findall(r"TITEL\s*=\s*PAW_PBE\s+(\S+)\s+(\S+)", text)
        elements = re.findall(r"VRHFIN\s*=\s*([A-Za-z]+)", text)
        for (variant, date), el in zip(titels, elements):
            escaped_variant = variant.replace("_", r"\_")
            out.setdefault(el, f"{escaped_variant} ({date})")
    return out


def build_common_table(lang: str) -> str:
    sample_incar = parse_incar(next(STRUCTURES_ROOT.iterdir()) / "INCAR")
    header = {"fr": ("Param\\`etre", "Valeur", "Description"),
              "en": ("Parameter", "Value", "Description")}[lang]
    lines = [r"\begin{tabular}{@{}l l p{7.5cm}@{}}", r"\toprule",
             " & ".join(header) + r" \\", r"\midrule"]
    for key in COMMON_KEYS:
        val = sample_incar.get(key, "?")
        desc = DESCRIPTION[lang][key]
        lines.append(f"\\texttt{{{key}}} & {val} & {desc} \\\\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def build_potcar_table(lang: str) -> str:
    titels = gather_potcar_titels()
    header = {"fr": ("\\'El\\'ement", "Potentiel PAW-PBE"),
              "en": ("Element", "PAW-PBE potential")}[lang]
    caption = {
        "fr": "Potentiels PAW-PBE utilis\\'es, un par \\'el\\'ement pr\\'esent dans le jeu de donn\\'ees (biblioth\\`eque MPRelaxSet de pymatgen, avec une exception document\\'ee : W substitu\\'e \\`a W\\_pv/W\\_sv, voir \\S7).",
        "en": "PAW-PBE potentials used, one per element present in the dataset (pymatgen's MPRelaxSet library, with one documented exception: W substituted for W\\_pv/W\\_sv, see \\S7).",
    }[lang]
    lines = [r"\begin{longtable}{@{}l l@{}}",
             rf"\caption{{{caption}}}\\",
             r"\toprule", " & ".join(header) + r" \\", r"\midrule",
             r"\endfirsthead", r"\toprule", " & ".join(header) + r" \\", r"\midrule", r"\endhead",
             r"\bottomrule", r"\endfoot"]
    for el in sorted(titels):
        lines.append(f"{el} & \\texttt{{{titels[el]}}} \\\\")
    lines.append(r"\end{longtable}")
    return "\n".join(lines)


def build_kpoints_table(lang: str) -> str:
    df = pd.read_csv(CSV_PATH).sort_values(["family", "energy_above_hull_eV_at"])
    header = {"fr": ("Formule", "mp-id", "Grille k-points", "NBANDS"),
              "en": ("Formula", "mp-id", "k-point grid", "NBANDS")}[lang]
    caption = {
        "fr": "Grille de points k (densit\\'e automatique \\texttt{pymatgen}, \\texttt{kppvol=100}) et nombre de bandes (d\\'eriv\\'e de la base LOBSTER) par compos\\'e.",
        "en": "k-point grid (\\texttt{pymatgen} automatic density, \\texttt{kppvol=100}) and number of bands (derived from the LOBSTER basis) per compound.",
    }[lang]
    continued = {"fr": "(suite)", "en": "(continued)"}[lang]

    lines = [r"\begin{longtable}{@{}l l l r@{}}",
             rf"\caption{{{caption}}}\\",
             r"\toprule", " & ".join(header) + r" \\", r"\midrule",
             r"\endfirsthead",
             rf"\multicolumn{{4}}{{@{{}}l}}{{\tablename\ \thetable\ {continued}}}\\",
             r"\toprule", " & ".join(header) + r" \\", r"\midrule", r"\endhead",
             r"\bottomrule", r"\endfoot"]

    for family in FAMILY_ORDER:
        sub = df[df["family"] == family]
        flabel = FAMILY_LABEL[lang][family]
        lines.append(rf"\multicolumn{{4}}{{@{{}}l}}{{\textbf{{{flabel}}} ($n={len(sub)}$)}} \\")
        lines.append(r"\midrule")
        for _, row in sub.iterrows():
            compound_dir = STRUCTURES_ROOT / row["compound_id"]
            incar = parse_incar(compound_dir / "INCAR")
            mesh_type, grid = parse_kpoints(compound_dir / "KPOINTS")
            nbands = incar.get("NBANDS", "?")
            grid_str = f"{grid} ({mesh_type[0]})"  # (G)amma or (M)onkhorst
            lines.append(f"{row['formula']} & {row['mp_id']} & {grid_str} & {nbands} \\\\")
        lines.append(r"\midrule")
    lines.append(r"\end{longtable}")
    return "\n".join(lines)


def main():
    for lang in ("fr", "en"):
        (HERE / f"appendix_dft_common_{lang}.tex").write_text(build_common_table(lang))
        (HERE / f"appendix_dft_potcar_{lang}.tex").write_text(build_potcar_table(lang))
        (HERE / f"appendix_dft_kpoints_{lang}.tex").write_text(build_kpoints_table(lang))
    print("wrote appendix_dft_{common,potcar,kpoints}_{fr,en}.tex")


if __name__ == "__main__":
    main()
