# -*- coding: utf-8 -*-

import os
import unittest

from obadb.analyzer.chempotdiag.chem_pot_diag import ChemPotDiag

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.corrections.corrections import ExtendedFnvCorrection
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectEnergiesTest(unittest.TestCase):

    def setUp(self):
        """ """
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        unitcell = UnitcellCalcResults.load_json(unitcell_file)
        perfect_file = os.path.join(test_dir,
                                    "MgO/defects/perfect/dft_results.json")
        perfect = SupercellCalcResults.load_json(perfect_file)

        defect_dirs = ["Va_O1_1", "Va_O1_2", "Va_O1_0"]
        defects = []
        for dd in defect_dirs:
            d = os.path.join(test_dir, "MgO/defects", dd)
            defect_entry = \
                DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
            dft_results = \
                SupercellCalcResults.load_json(
                    os.path.join(d, "dft_results.json"))
            correction = \
                ExtendedFnvCorrection.load_json(os.path.join(d, "correction.json"))

            defect = Defect(defect_entry=defect_entry,
                            dft_results=dft_results,
                            correction=correction)

            defects.append(defect)

        # temporary insert values
        chem_pot = ChemPotDiag.load_vertices_yaml(
            os.path.join(test_dir, "MgO/vertices_MgO.yaml"))

        chem_pot_label = "A"

        self.defect_energies = \
            DefectEnergies.from_files(unitcell=unitcell,
                                      perfect=perfect,
                                      defects=defects,
                                      chem_pot=chem_pot,
                                      chem_pot_label=chem_pot_label,
                                      system="MgO")

    def test_print(self):
        print(self.defect_energies)

    def test_energies(self):
        d = self.defect_energies.as_dict()
        de = DefectEnergies.from_dict(d)
        dd = de.as_dict()
        self.assertEqual(d, dd)

    def test_multiplicity(self):
        actual = self.defect_energies.multiplicity["Va_O1"][2][0]
        expected = 8
        self.assertEqual(actual, expected)

    def test_U(self):
        actual = self.defect_energies.u(name="Va_O1", charge=[0, 1, 2])
        expected = 1.82926181856907
        self.assertAlmostEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()

