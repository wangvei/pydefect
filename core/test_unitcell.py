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

TEST_DIRECTORY = "../examples/MgO/unitcell/structure_optimization"
TEST_DIELECTRICS_DIRECTORY = "../examples/MgO/unitcell/dielectric_constants"
DEFECT_JSON = TEST_DIRECTORY + "/unitcell.json"


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
        self._u.to_json_file("test_unitcell.json")
        self._u_json = Unitcell.json_load("test_unitcell.json")
        self.assertEqual(self._u.final_structure, self._u_json.final_structure)
        self.assertEqual(self._u.total_energy, self._u_json.total_energy)
        np.testing.assert_equal(self._u.eigenvalues[Spin.up],
                                self._u_json.eigenvalues[Spin.up])

    def test_set_dielectric_constants(self):
        expected_static = np.array([[3.166727, 0, 0],
                                    [0, 3.166727, 0],
                                    [0, 0, 3.166727]])
        expected_ionic = np.array([[9.102401, 0, 0],
                                   [0, 9.102448, 0],
                                   [0, 0, 9.102542]])
        expected_total = np.array([[12.269128, 0, 0],
                                   [0, 12.269175, 0],
                                   [0, 0, 12.269269]])
        self._u.set_dielectric_constants_from_outcar(TEST_DIELECTRICS_DIRECTORY)

        np.testing.assert_equal(self._u.static_dielectric_tensor,
                                expected_static)
        np.testing.assert_equal(self._u.ionic_dielectric_tensor, expected_ionic)
        np.testing.assert_equal(self._u.total_dielectric_tensor, expected_total)

