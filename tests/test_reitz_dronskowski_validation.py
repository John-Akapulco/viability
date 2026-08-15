"""End-to-end validation of reaction_analysis (balance.py + delta.py +
classify.py) against the 7 worked examples in Reitz & Dronskowski,
ic-2026-04181q ("Calculus of Bonding Energetics" section) -- external,
published numbers independent of any viability CSP data, exercised
without touching parse_lobster.py (fixtures encode already-summed
IcohpSummary values directly; see fixtures/reitz_dronskowski_cases.yaml's
header for why).

Worked arithmetic for all 7 cases (products - reactants, eV, before
kJ/mol conversion) -- kept here rather than only in the YAML comments so
a reviewer can check the fixture's sum_total_eV values without cross-
referencing two files:

  1. Pb(N3)2 -> Pb + 3 N2
     reactant = -89.104; products = -5.688 + 3*(-23.161) = -75.171
     delta = -75.171 - (-89.104) = +13.933 eV = +1344.3 kJ/mol (manuscript: +1345)

  2. S4N2 -> 1/2 S8 + N2
     reactant = -49.24; products = 0.5*(-46.8) + (-23.161) = -46.561
     delta = -46.561 - (-49.24) = +2.679 eV = +258.5 kJ/mol (manuscript: +258)

  3. S4N4 -> 1/2 S8 + 2 N2
     reactant = -73.84; products = 0.5*(-46.8) + 2*(-23.161) = -69.722
     delta = -69.722 - (-73.84) = +4.118 eV = +397.3 kJ/mol (manuscript: +399)

  4. ZnSn -> Zn + Sn
     reactant = -16.869; products = -12.686 + (-7.672) = -20.358
     delta = -20.358 - (-16.869) = -3.489 eV = -336.7 kJ/mol (manuscript: -337)

  5. CaO[sphalerite] -> CaO[rocksalt]
     reactant = -3.46; product = -4.278
     delta = -4.278 - (-3.46) = -0.818 eV = -78.9 kJ/mol (manuscript: -79)

  6. CaN -> 1/3 Ca3N2 + 1/6 N2
     reactant = -4.302; products = (1/3)*(-7.692) + (1/6)*(-23.161) = -2.564 - 3.860 = -6.424
     delta = -6.424 - (-4.302) = -2.122 eV = -204.8 kJ/mol (manuscript: -205)

  7. Mn2O7 -> 2 MnO2 + 3/2 O2
     reactant = -62.21; products = 2*(-18.54) + 1.5*(-18.05) = -37.08 - 27.075 = -64.155
     delta = -64.155 - (-62.21) = -1.945 eV = -187.7 kJ/mol (manuscript: -186)

Residuals (up to ~1.7 kJ/mol) are the manuscript's own rounding of its
reported per-bond eV values to 2-3 decimals propagating through the
reconstruction -- not slack for an implementation bug; each test's
tolerance is set from tolerance_kJ_per_mol in the YAML, not padded ad hoc.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import yaml

from reaction_analysis.classify import BondingLabel
from reaction_analysis.delta import compute_delta
from reaction_analysis.schema import CompoundEntry, IcohpSummary, Reaction, ReactionMember
from reaction_analysis.units import ev_to_kj_per_mol

FIXTURE_PATH = Path(__file__).parent / "fixtures" / "reitz_dronskowski_cases.yaml"


def _load_fixture() -> dict:
    return yaml.safe_load(FIXTURE_PATH.read_text())


def _build_entry(compound: dict) -> CompoundEntry:
    composition = compound["composition"]
    n_atoms = sum(composition.values())
    sum_total = compound["sum_total_eV"]
    icohp = IcohpSummary(
        sum_total_eV=sum_total,
        sum_per_atom_eV=sum_total / n_atoms,
        sum_per_formula_unit_eV=sum_total,  # Z == 1 for every fixture compound
        mean_per_bond_eV=sum_total,
        n_bonds=1,
    )
    return CompoundEntry(
        compound_id=compound["id"],
        formula=compound.get("formula", compound["id"]),
        composition=composition,
        Z=1,
        space_group_symbol=compound.get("space_group_symbol", "unknown"),
        space_group_number=compound.get("space_group_number", 0),
        role=compound["role"],
        icohp=icohp,
        source_path=f"<fixture:{FIXTURE_PATH.name}#{compound['id']}>",
    )


class TestReitzDronskowskiCases(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = _load_fixture()
        cls.entries = {c["id"]: _build_entry(c) for c in cls.fixture["compounds"]}

    @staticmethod
    def _member(m: dict) -> ReactionMember:
        return ReactionMember(compound_id=m["id"], coefficient=m["coefficient"])

    def _run(self, reaction_fixture: dict):
        reaction = Reaction(
            reaction_id=reaction_fixture["id"],
            type=reaction_fixture["type"],
            reactants=[self._member(m) for m in reaction_fixture["reactants"]],
            products=[self._member(m) for m in reaction_fixture["products"]],
        )
        results = compute_delta(reaction, self.entries)
        self.assertEqual(len(results), 1)
        return results[0]

    def test_all_seven_cases(self):
        for reaction_fixture in self.fixture["reactions"]:
            with self.subTest(case=reaction_fixture["id"]):
                result = self._run(reaction_fixture)

                self.assertIsNone(result.error, f"{reaction_fixture['id']}: {result.error}")

                # (a) numeric ΔICOHP, in kJ/mol, within the fixture's own tolerance
                got_kJ = ev_to_kj_per_mol(result.delta_per_formula_unit_eV)
                expected_kJ = reaction_fixture["expected_delta_kJ_per_mol"]
                tol = reaction_fixture["tolerance_kJ_per_mol"]
                self.assertAlmostEqual(
                    got_kJ, expected_kJ, delta=tol,
                    msg=f"{reaction_fixture['id']}: got {got_kJ:.2f} kJ/mol, expected {expected_kJ} +/- {tol}",
                )

                # (b) BondingLabel
                expected_label = BondingLabel(reaction_fixture["expected_label"])
                self.assertEqual(
                    result.bonding_label, expected_label.value,
                    f"{reaction_fixture['id']}: expected {expected_label}, got {result.bonding_label}",
                )

                # (c) balance.py accepted the reaction as written (no
                # balance error was raised -- already implied by
                # result.error being None above, asserted explicitly here
                # for readability of intent)
                self.assertIsNone(result.error)


if __name__ == "__main__":
    unittest.main()
