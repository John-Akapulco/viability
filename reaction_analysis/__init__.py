"""Schema-driven ΔICOHP/ΔICOBI reaction analysis.

Distinct from, and not a replacement for, the ad hoc `reaction_icohp.py`
module used in mission #5 (still valid, still used by
`analysis/compute_reaction_icohp_case1.py` and case 2). This package is a
from-scratch redesign around a formal Pydantic schema (`schema.py`) meant
to eventually cover all three reaction types mission #5 only partially
reached (case 3 was judged not tractable there without new DFT -- the
schema itself does not depend on that being solved).

No chemical roles, thresholds, or decision rules are hardcoded anywhere in
this package -- every judgment is driven by the input data
(`CompoundEntry.role`, `Reaction` members/coefficients), the same
principle already in place in the `chiral_Mat` project.

As of this commit: schema + parsing + balance-checking + the three
delta normalizations, validated on synthetic fixtures only. No real
LOBSTER production data has been run through this package yet -- see
each module's docstring for exactly what is and isn't covered.
"""
