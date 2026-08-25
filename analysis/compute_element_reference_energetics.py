"""Elemental reference structures (analysis/compute_reaction_icohp_case1.py's
ELEMENT_REFERENCE, 62 entries -- one structure per element used as the
"products" side of every case-1 decomposition reaction) with, for each:
space group, MP's own tabulated total energy per atom (fetched live via
the Materials Project API, same mp_api.client.MPRester convention as
mp_dataset/download_*.py), this campaign's own VASP total energy per atom
(parsed from OSZICAR's final ionic step, E0 line -- the extrapolated
T->0 energy, same convention VASP itself recommends for comparing
total energies), and the element's own ICOHP/ICOBI and antibonding-ICOHP
descriptors (already computed dataset-wide in icohp_icobi_bondtype.csv /
icohp_antibonding_full.csv).

This is a DFT-pipeline consistency check (does our own relaxation of the
same structure converge to an energy consistent with MP's tabulated
value?), not a formation-energy comparison -- an element's own formation
enthalpy is 0 by definition (it is its own reference state), so that
comparison would be trivial and uninformative.

2 of the 62 references (N: gasref_N2_dimerbox; and, if present, any
other hand-built gasref_* box) have no Materials Project entry at all
(isolated dimer built by hand to match a specific manuscript methodology,
not a periodic solid with a meaningful space group) -- reported with
MP columns as "--".

Writes analysis/element_reference_energetics.csv and, via
report/gen_appendix_elements.py (separate script), the SI table.
"""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

import pandas as pd
from mp_api.client import MPRester

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "analysis"))
from compute_reaction_icohp_case1 import ELEMENT_REFERENCE  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
API_KEY_PATH = Path(os.path.expanduser("~/.mp_api_key"))
OUT_CSV = Path(__file__).parent / "element_reference_energetics.csv"

_OSZICAR_E0_RE = re.compile(r"E0=\s*([-0-9.Ee+]+)")


def _parse_oszicar_final_e0(oszicar_path: Path) -> float | None:
    """Last ionic-step line's E0 (extrapolated T->0 total energy, eV, whole cell)."""
    last = None
    for line in oszicar_path.read_text(errors="replace").splitlines():
        m = _OSZICAR_E0_RE.search(line)
        if m:
            last = float(m.group(1))
    return last


def _nsites(meta: dict, contcar_path: Path) -> int | None:
    for key in ("nsites", "num_sites", "num_sites_structure"):
        if meta.get(key):
            return int(meta[key])
    if contcar_path.exists():
        lines = contcar_path.read_text().splitlines()
        counts = [int(x) for x in lines[6].split()]
        return sum(counts)
    return None


def main() -> None:
    rows = []
    mp_ids_needed = []
    for element, dirname in sorted(ELEMENT_REFERENCE.items()):
        d = STRUCTURES_ROOT / dirname
        meta_path = d / "mp_metadata.json"
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        mp_id = meta.get("mp_id") or meta.get("material_id")
        nsites = _nsites(meta, d / "CONTCAR")
        e0_total = _parse_oszicar_final_e0(d / "OSZICAR")
        row = {
            "element": element,
            "compound_id": dirname,
            "mp_id": mp_id,
            "spacegroup": meta.get("spacegroup") or meta.get("expected_spacegroup"),
            "nsites": nsites,
            "campaign_energy_per_atom_eV": (e0_total / nsites) if (e0_total is not None and nsites) else None,
            "mp_energy_per_atom_eV": None,
        }
        rows.append(row)
        if mp_id:
            mp_ids_needed.append(mp_id)

    print(f"Fetching MP energy_per_atom for {len(mp_ids_needed)} material ids...")
    mp_energy = {}
    if mp_ids_needed and API_KEY_PATH.exists():
        api_key = API_KEY_PATH.read_text().strip()
        with MPRester(api_key) as mpr:
            docs = mpr.materials.summary.search(
                material_ids=mp_ids_needed,
                fields=["material_id", "energy_per_atom"],
            )
        mp_energy = {str(dd.material_id): dd.energy_per_atom for dd in docs}
    else:
        print("WARNING: no API key found or no mp_ids to fetch; MP energies left blank.")

    for row in rows:
        if row["mp_id"]:
            row["mp_energy_per_atom_eV"] = mp_energy.get(row["mp_id"])

    df = pd.DataFrame(rows)

    bondtype = pd.read_csv(REPO_ROOT / "analysis" / "icohp_icobi_bondtype.csv")[
        ["compound_id", "icohp_mean", "icobi_mean"]
    ]
    antibond = pd.read_csv(REPO_ROOT / "analysis" / "icohp_antibonding_full.csv")[
        ["compound_id", "antibond_w_raw", "antibond_w_normalized"]
    ]
    df = df.merge(bondtype, on="compound_id", how="left").merge(antibond, on="compound_id", how="left")

    df["mp_vs_campaign_diff_meV_per_atom"] = (
        (df["campaign_energy_per_atom_eV"] - df["mp_energy_per_atom_eV"]) * 1000.0
    )

    df = df.sort_values("element").reset_index(drop=True)
    df.to_csv(OUT_CSV, index=False)

    n_mp_matched = df["mp_energy_per_atom_eV"].notna().sum()
    n_icohp = df["icohp_mean"].notna().sum()
    n_antibond = df["antibond_w_raw"].notna().sum()
    print(f"{len(df)} elements written to {OUT_CSV}")
    print(f"  with MP energy_per_atom: {n_mp_matched}")
    print(f"  with ICOHP data: {n_icohp}")
    print(f"  with antibonding data: {n_antibond}")
    if n_mp_matched:
        diffs = df["mp_vs_campaign_diff_meV_per_atom"].dropna()
        print(f"  MP-vs-campaign energy diff: mean={diffs.mean():.2f} meV/at, "
              f"max|.|={diffs.abs().max():.2f} meV/at, median={diffs.median():.2f} meV/at")


if __name__ == "__main__":
    main()
