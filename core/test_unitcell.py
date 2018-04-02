#!/usr/bin/env python

import numpy as np
import unittest

from pymatgen.electronic_structure.core import Spin

from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.core.unitcell import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

TEST_DIRECTORY = "../examples/MgO/unitcell"
DEFECT_JSON = DIRNAME_VAC + "/unitcell.json"


class UnitcellTest(unittest.TestCase):

    def setUp(self):
        """ """
        self._u = Unitcell.from_vasp_results(TEST_DIRECTORY)
        self._u_set = Unitcell()

    def test_set_vasp_results(self):
        self._u_set.set_vasp_results(TEST_DIRECTORY)
        self.assertEqual(self._u.final_structure, self._u_set.final_structure)
        self.assertEqual(self._u.total_energy, self._u_set.total_energy)
        np.testing.assert_equal(self._u.eigenvalues[Spin.up],
                                self._u_set.eigenvalues[Spin.up])

    def test_json(self):
        self._p.to_json_file("test_unitcell.json")
        self._p_json = Perfect.json_load("test_unitcell.json")


class DefectTest(unittest.TestCase):

    def setUp(self):

        self._p = Perfect.from_vasp_results(DIRNAME_PERFECT)
        self._d = Defect.from_vasp_results(DIRNAME_VAC)
        defect_entry = DefectEntry.json_load(DEFECT_JSON)
        self._d_set = Defect(defect_entry)

    def test_set_vasp_results(self):
        self._d_set.set_vasp_results(DIRNAME_VAC)
        print(self._d_set.charge)
        print(self._d_set.total_energy)
        self._d_set.set_relative_values(self._p)
        print(self._d_set.relative_total_energy)

    

#        self.assertEqual(self._d.final_structure, self._d_set.final_structure)
#        self.assertEqual(self._d.total_energy, self._d_set.total_energy)
#        np.testing.assert_equal(self._d.eigenvalues[Spin.up],
#                                self._d_set.eigenvalues[Spin.up])
#        self.assertEqual(self._d.electrostatic_potential,
#                         self._d_set.electrostatic_potential)

#    def test_set_relative_values(self):


if __name__ == "__main__":
    unittest.main()

