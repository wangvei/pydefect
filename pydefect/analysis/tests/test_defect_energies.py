# -*- coding: utf-8 -*-

from collections import namedtuple
import numpy as np
import os
import tempfile
import unittest

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.core.correction import Correction
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "Feb. 25, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectEnergiesTest(unittest.TestCase):

    def setUp(self):
        """ """
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        unitcell = UnitcellDftResults.load_json(unitcell_file)
        perfect_file = os.path.join(test_dir,
                                    "MgO/defects/perfect/dft_results.json")
        perfect = SupercellDftResults.load_json(perfect_file)

        # defect_dirs = ["Va_O1_1", "Va_O1_2", "Va_O1_0"]
        defect_dirs = ["Mg_O1_0", "Mg_O1_1", "Mg_O1_2", "Mg_O1_3", "Mg_O1_4",
                       "Mg_i1_0", "Mg_i1_1", "Mg_i1_2", "Va_O1_1", "Va_O1_2",
                       "Va_O1_0"]
        defects = []
        for dd in defect_dirs:
            d = os.path.join(test_dir, "MgO/defects", dd)
            defect_entry = \
                DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
            dft_results = \
                SupercellDftResults.load_json(
                    os.path.join(d, "dft_results.json"))
            correction = \
                Correction.load_json(os.path.join(d, "correction.json"))

            defect = Defect(defect_entry=defect_entry,
                            dft_results=dft_results,
                            correction=correction)

            defects.append(defect)

        # temporary insert values
        chem_pot = {"A": {"Mg": -2.1, "O": 0}, "B": {"Mg": 0, "O": -2.1}}
        chem_pot_label = "A"

        self.defect_energies = DefectEnergies(unitcell=unitcell,
                                              perfect=perfect,
                                              defects=defects,
                                              chem_pot=chem_pot,
                                              chem_pot_label=chem_pot_label,
                                              system_name="MgO")
#        filtering_words=["Va_O"],

    def test_energies(self):
        print(self.defect_energies._defect_energies)
        print(self.defect_energies.vbm)
        print(self.defect_energies.cbm)
        print(self.defect_energies.band_gap)

    def test_calc_transition_levels(self):
        self.defect_energies.calc_transition_levels()
        print(self.defect_energies._transition_levels)
        self.defect_energies.plot_energy()
#        self.defect_energies.plot_energy(x_range=[-0.5, 5], y_range=[-5, 20])
#        self.defect_energies.plot_energy(file_name="test.eps")


if __name__ == "__main__":
    unittest.main()

