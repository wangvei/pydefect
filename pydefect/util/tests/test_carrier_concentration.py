# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.analysis.defect_carrier_concentration import CarrierConcentration
from pydefect.core.unitcell_calc_results import UnitcellCalcResults

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


# class CarrierConcentrationTest(unittest.TestCase):

    # def setUp(self):
    #     self.unitcell = UnitcellCalcResults.load_json("unitcell_MgSe.json")
    #     print(self.unitcell)

    # def test(self):
    #     t = 1298
    #     cc = CarrierConcentration.from_unitcell(t, self.unitcell)

# #        print(cc)
#         cc.get_plot()
