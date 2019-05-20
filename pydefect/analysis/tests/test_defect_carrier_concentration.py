# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.analysis.defect_carrier_concentration import DefectConcentration
from pydefect.analysis.chempotdiag.chem_pot_diag import ChemPotDiag
from pydefect.corrections.corrections import ExtendedFnvCorrection
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "Feb. 25, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectConcentrationTest(unittest.TestCase):

    def setUp(self):
        """ """
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        self.unitcell = UnitcellCalcResults.load_json(unitcell_file)
        perfect_file = os.path.join(test_dir,
                                    "MgO/defects/perfect/dft_results.json")
        perfect = SupercellCalcResults.load_json(perfect_file)

        defect_dirs = ["Mg_O1_0", "Mg_O1_1", "Mg_O1_2", "Mg_O1_3", "Mg_O1_4",
                       "Mg_i1_0", "Mg_i1_1", "Mg_i1_2", "Va_O1_1", "Va_O1_2",
                       "Va_O1_0"]
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
            DefectEnergies.from_files(unitcell=self.unitcell,
                                      perfect=perfect,
                                      defects=defects,
                                      chem_pot=chem_pot,
                                      chem_pot_label=chem_pot_label,
                                      system="MgO")

    def test_from_defect_energies(self):
        temperature = 10000
        num_sites_filename = os.path.join(test_dir,
                                          "MgO/defects/num_sites.yaml")

        dc1 = DefectConcentration.from_defect_energies(
            defect_energies=self.defect_energies,
            temperature=temperature,
            unitcell=self.unitcell,
            num_sites_filename=num_sites_filename)

        print(dc1.energies)
        print(dc1.temperature)
        print(dc1.e_f)
        print(p)
        print(n)
        print(dc1.concentration)

        temperature2 = 1000

        # dc2 = DefectConcentration.from_defect_energies(
        #     defect_energies=self.defect_energies,
        #     temperature=temperature2,
        #     unitcell=self.unitcell,
        #     num_sites_filename=num_sites_filename,
        #     previous_concentration=dc1,
        #     verbose=True)

        # print("-------------------------------")
        # print(dc2.energies)
        # print(dc2.temperature)
        # print(dc2.e_f)
        # print(dc2.p)
        # print(dc2.n)
        # print(dc2.concentration)


if __name__ == "__main__":
    unittest.main()

