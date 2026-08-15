"""Unit conversion between eV and kJ/mol. Every other module in this
package works exclusively in eV (matching the rest of `viability`); this
conversion is only needed at the boundary with literature values reported
in kJ/mol (e.g. Reitz & Dronskowski, ic-2026-04181q) -- the constant lives
here so it is never hardcoded a second time elsewhere.
"""

from __future__ import annotations

EV_TO_KJ_PER_MOL = 96.4853075


def ev_to_kj_per_mol(value_eV: float) -> float:
    return value_eV * EV_TO_KJ_PER_MOL


def kj_per_mol_to_ev(value_kJ_per_mol: float) -> float:
    return value_kJ_per_mol / EV_TO_KJ_PER_MOL
