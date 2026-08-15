"""LOBSTER band-overlap quality check across the 74 extension_* compounds --
the check flagged as not-yet-done in METRIC_DEFINITION_reaction_icohp.md
section 6.

Parses bandOverlaps.lobster's "maxDeviation is: <value>" lines (one per
k-point; only written by LOBSTER for k-points whose orthonormalized
projected-band overlap matrix deviates from identity beyond its internal
threshold -- most well-behaved calculations have no bandOverlaps.lobster
file at all). Reports the worst (max) maxDeviation per compound.

Threshold: >0.1 flagged as a real projection-quality concern, matching the
order of magnitude of the two previously-known problem cases
(extension_CaN_mp-1058549, extension_CaO_mp-2605, both ~15-17 -- see
project-viability-extension-campaign memory). No LOBSTER-recommended
official cutoff is assumed; 0.1 is this script's own threshold, chosen to
be well above LOBSTER's typical well-behaved-calculation range (<0.03,
seen in the 6-pilot dataset) and well below the worst known cases.

Writes analysis/lobster_quality_extension.json.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
STRUCTURES_ROOT = REPO_ROOT / "mp_dataset" / "structures"
OUT = Path(__file__).parent / "lobster_quality_extension.json"
FLAG_THRESHOLD = 0.1

# Elements this compound is used as the elemental reference for, per
# compute_reaction_icohp_case1.py's ELEMENT_REFERENCE dict -- duplicated
# here (not imported) since that dict is keyed by element symbol, not
# compound_dir name; small enough to invert by hand and worth keeping
# independent of that script's internals for this standalone check.
import sys  # noqa: E402
sys.path.insert(0, str(REPO_ROOT / "analysis"))
from compute_reaction_icohp_case1 import ELEMENT_REFERENCE  # noqa: E402

DIR_TO_ELEMENTS = {}
for element, dirname in ELEMENT_REFERENCE.items():
    DIR_TO_ELEMENTS.setdefault(dirname, []).append(element)


def main():
    results = []
    for d in sorted(STRUCTURES_ROOT.glob("extension_*")):
        bo = d / "bandOverlaps.lobster"
        if not bo.exists():
            results.append({"compound_dir": d.name, "max_deviation": None, "flagged": False,
                             "note": "no bandOverlaps.lobster file (LOBSTER only writes one when a k-point exceeds its internal threshold)"})
            continue
        text = bo.read_text()
        devs = [float(x) for x in re.findall(r"maxDeviation is:\s*([\d.eE+-]+)", text)]
        max_dev = max(devs) if devs else None
        results.append({
            "compound_dir": d.name,
            "max_deviation": max_dev,
            "n_flagged_kpoints": len(devs),
            "flagged": max_dev is not None and max_dev > FLAG_THRESHOLD,
            "used_as_elemental_reference_for": DIR_TO_ELEMENTS.get(d.name, []),
        })

    n_with_file = sum(1 for r in results if r["max_deviation"] is not None)
    n_flagged = sum(1 for r in results if r["flagged"])
    n_flagged_and_reference = sum(1 for r in results if r["flagged"] and r["used_as_elemental_reference_for"])

    summary = {
        "threshold": FLAG_THRESHOLD,
        "n_total_extension_dirs": len(results),
        "n_with_bandOverlaps_file": n_with_file,
        "n_flagged": n_flagged,
        "n_flagged_used_as_elemental_reference": n_flagged_and_reference,
        "flagged_compounds": sorted(
            [r for r in results if r["flagged"]], key=lambda r: -r["max_deviation"]
        ),
        "all_results": results,
    }
    OUT.write_text(json.dumps(summary, indent=2, default=str))
    print(f"n_total={len(results)} n_with_file={n_with_file} n_flagged(>{FLAG_THRESHOLD})={n_flagged} "
          f"n_flagged_and_used_as_reference={n_flagged_and_reference}")
    for r in summary["flagged_compounds"]:
        ref_note = f" [reference for {r['used_as_elemental_reference_for']}]" if r["used_as_elemental_reference_for"] else ""
        print(f"  {r['compound_dir']}: max_deviation={r['max_deviation']:.4f}{ref_note}")


if __name__ == "__main__":
    main()
