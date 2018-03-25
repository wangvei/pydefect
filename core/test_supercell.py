#!/usr/bin/env python

import unittest
import numpy as np
from supercell import *
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.electronic_structure.core import Spin
from pydefect.input_maker.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


TEST_DIRECTORY = "../examples/MgO"
DIRNAME_PERFECT = TEST_DIRECTORY + "/defects/perfect"
DIRNAME_VAC = TEST_DIRECTORY + "/defects/Va_O1_2"


class PerfectTest(unittest.TestCase):

    def setUp(self):
        """ """
        self._p = Perfect.from_vasp_results(DIRNAME_PERFECT)
        self._p_set = Perfect()

    def test_set_vasp_results(self):
        self._p_set.set_vasp_results(DIRNAME_PERFECT)
        self.assertEqual(self._p.final_structure, self._p_set.final_structure)
        self.assertEqual(self._p.total_energy, self._p_set.total_energy)
        np.testing.assert_equal(self._p.eigenvalues[Spin.up],
                                self._p_set.eigenvalues[Spin.up])
        self.assertEqual(self._p.electrostatic_potential,
                         self._p_set.electrostatic_potential)

    def test_json(self):
        self._p.to_json_file("test_supercell.json")
        self._p_json = Perfect.json_load("test_supercell.json")


class DefectTest(unittest.TestCase):

    def setUp(self):

        self._d = Defect.from_vasp_results(DIRNAME_VAC)
        self._d_set = Defect()
        self.defect_entry = DefectEntry.json_load()

    def test_set_vasp_results(self):
        self._d_set.set_vasp_results(DIRNAME_VAC)
        self.assertEqual(self._d.final_structure, self._d_set.final_structure)
        self.assertEqual(self._d.total_energy, self._d_set.total_energy)
        np.testing.assert_equal(self._d.eigenvalues[Spin.up],
                                self._d_set.eigenvalues[Spin.up])
        self.assertEqual(self._d.electrostatic_potential,
                         self._d_set.electrostatic_potential)


if __name__ == "__main__":
    unittest.main()

