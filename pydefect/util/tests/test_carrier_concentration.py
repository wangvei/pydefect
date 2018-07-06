# -*- coding: utf-8 -*-

import numpy as np
import os
import unittest

from pydefect.util.carrier_concentration import fermi_dirac_dist, \
    bose_einstein_dist, maxwell_boltzmann_dist, CarrierConcentration
from pydefect.core.unitcell_dft_results import UnitcellDftResults

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "July 4, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class CarrierConcentrationTest(unittest.TestCase):

    def setUp(self):
        self.unitcell = UnitcellDftResults.load_json("unitcell_MgSe.json")
        print(self.unitcell)

    def test(self):
        t = 1298
        cc = CarrierConcentration.from_unitcell(t, self.unitcell)

#        print(cc)
        cc.get_plot()
