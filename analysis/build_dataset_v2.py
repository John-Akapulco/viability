"""Extend the dataset with the two new mission-#3 descriptors (network
dimensionality, periodic min-cut) and formation_energy_per_atom, joined
onto the existing percolation_vs_hull.csv columns.

Separate script rather than editing build_dataset.py in place, since the
column set is different enough (dimensionality at 3 thetas, min-cut per
direction x metric) that extending in place would make both harder to
read. Reuses percolation_path.load_bonds and the mission-#3 modules
(network_dimensionality, periodic_mincut) unchanged -- no duplicated
parsing or algorithm logic.

Output: analysis/percolation_vs_formation_energy.csv
"""

import json
import sys
from pathlib import Path

import pandas as pd
from pymatgen.core import Structure

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "Percolation_viability"))
from percolation_path import load_bonds, max_translation_extent  # noqa: E402
import network_dimensionality as nd  # noqa: E402
import periodic_mincut as pmc  # noqa: E402

STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
PERCOLATION_CSV = Path(__file__).parent / "percolation_vs_hull.csv"
FORMATION_ENERGY_JSON = REPO_ROOT / "mp_dataset" / "formation_energies.json"
OUT_CSV = Path(__file__).parent / "percolation_vs_formation_energy.csv"

THETAS = {"dimensionality": 0.10, "dimensionality_theta05": 0.05, "dimensionality_theta20": 0.20}


def compute_row_extras(compound_id: str) -> dict:
    compound_dir = STRUCTURES_ROOT / compound_id
    contcar = compound_dir / "CONTCAR"
    icohp_path = compound_dir / "ICOHPLIST.lobster"
    icobi_path = compound_dir / "ICOBILIST.lobster"

    structure = Structure.from_file(str(contcar))
    num_sites = len(structure)
    bonds = load_bonds(
        icohp_path if icohp_path.exists() else None,
        icobi_path if icobi_path.exists() else None,
        num_sites=num_sites,
    )
    coord_bound = max(3, 2 * max_translation_extent(bonds) + 1) if bonds else 3

    out = {}
    for col, theta in THETAS.items():
        dim = nd.network_dimensionality(bonds, num_sites, "icohp", theta=theta, coord_bound=coord_bound)
        out[col] = dim["dimensionality"]

    strongest_bond_magnitude = None
    for metric in ("icohp", "icobi"):
        cuts = pmc.all_directions_min_cut(bonds, num_sites, metric, coord_bound=coord_bound)
        for d in ("a", "b", "c"):
            out[f"mincut_{metric}_{d}"] = cuts[d]
        finite = [v for v in cuts.values() if v is not None]
        out[f"mincut_{metric}_min"] = min(finite) if finite else None
        if metric == "icohp":
            # |strongest bond| = max(|icohp|) = |min(raw icohp)| (most negative
            # raw value is the strongest bond, matching percolation_path.py's
            # compute_aggregates() "min" convention) -- NOT min(|icohp|), which
            # would be the *weakest* bond and silently invert the normalization.
            vals = [abs(b.icohp) for b in bonds if b.icohp is not None]
            strongest_bond_magnitude = max(vals) if vals else None

    if out.get("mincut_icohp_min") is not None and strongest_bond_magnitude:
        out["mincut_icohp_min_normalized"] = out["mincut_icohp_min"] / strongest_bond_magnitude
    else:
        out["mincut_icohp_min_normalized"] = None

    return out


def main():
    base = pd.read_csv(PERCOLATION_CSV)
    formation_energies = json.loads(FORMATION_ENERGY_JSON.read_text())
    base["formation_energy_per_atom"] = base["mp_id"].map(formation_energies)

    n_missing_fe = base["formation_energy_per_atom"].isna().sum()
    if n_missing_fe:
        print(f"WARNING: {n_missing_fe} compounds missing formation_energy_per_atom")

    extra_rows = []
    n_failed = 0
    for compound_id in base["compound_id"]:
        try:
            extra_rows.append({"compound_id": compound_id, **compute_row_extras(compound_id)})
        except Exception as exc:  # noqa: BLE001 - batch must not die on one bad compound
            print(f"FAILED {compound_id}: {type(exc).__name__}: {exc}")
            extra_rows.append({"compound_id": compound_id})
            n_failed += 1

    extra_df = pd.DataFrame(extra_rows)
    merged = base.merge(extra_df, on="compound_id", how="left")
    merged.to_csv(OUT_CSV, index=False)

    print(f"\nWrote {len(merged)} rows to {OUT_CSV} ({n_failed} failed extras)")
    print(merged["dimensionality"].value_counts(dropna=False).sort_index())


if __name__ == "__main__":
    main()
