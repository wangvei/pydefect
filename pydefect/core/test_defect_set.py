#!/usr/bin/env python

import unittest

from pydefect.core.defect_set import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


TEST_DIRECTORY = "../examples/MgO"

TEST_UNITCELL_DIRECTORY = TEST_DIRECTORY + "/unitcell/structure_optimization"
TEST_DIELECTRICS_DIRECTORY = "../examples/MgO/unitcell/dielectric_constants"
UNITCELL_JSON = TEST_UNITCELL_DIRECTORY + "/unitcell.json"

TEST_PERFECT_DIRECTORY = TEST_DIRECTORY + "/defects/perfect"
TEST_VAC_DIRECTORY = TEST_DIRECTORY + "/defects/Va_O1_2"
VAC_DEFECT_ENTRY_JSON = TEST_VAC_DIRECTORY + "/defect_entry_Va_O1_2.json"


class DefectSetTest(unittest.TestCase):

    def setUp(self):
        """ """
        _d = [DefectSupercell.from_vasp_files(TEST_VAC_DIRECTORY)]
        _u = Unitcell.from_vasp_files(TEST_UNITCELL_DIRECTORY)
        _p = PerfectSupercell.from_vasp_files(TEST_PERFECT_DIRECTORY)
        self._defect_set = DefectSet(unitcell=_u, perfect=_p, defects=_d)

    def test_make_dft_results_json_files(self):
        DefectSet.make_dft_results_json_files(TEST_UNITCELL_DIRECTORY,
                                              TEST_DIELECTRICS_DIRECTORY,
                                              TEST_PERFECT_DIRECTORY,
                                              [TEST_VAC_DIRECTORY])

    #     self._u_set.set_vasp_results(TEST_UNITCELL_DIRECTORY)
    #     self.assertEqual(self._u.final_structure, self._u_set.final_structure)
    #     self.assertEqual(self._u.total_energy, self._u_set.total_energy)
    #     np.testing.assert_equal(self._u.eigenvalues[Spin.up],
    #                             self._u_set.eigenvalues[Spin.up])

    # def test_json(self):
    #     self._u.to_json_file("test_unitcell.json")
    #     self._u_json = Unitcell.json_load("test_unitcell.json")
    #     self.assertEqual(self._u.final_structure, self._u_json.final_structure)
    #     self.assertEqual(self._u.total_energy, self._u_json.total_energy)
    #     np.testing.assert_equal(self._u.eigenvalues[Spin.up],
    #                             self._u_json.eigenvalues[Spin.up])

    # def test_set_dielectric_constants(self):
    #     expected_static = np.array([[3.166727, 0, 0],
    #                                 [0, 3.166727, 0],
    #                                 [0, 0, 3.166727]])
    #     expected_ionic = np.array([[9.102401, 0, 0],
    #                                [0, 9.102448, 0],
    #                                [0, 0, 9.102542]])
    #     expected_total = np.array([[12.269128, 0, 0],
    #                                [0, 12.269175, 0],
    #                                [0, 0, 12.269269]])
    #     self._u.set_dielectric_constants_from_outcar(TEST_DIELECTRICS_DIRECTORY)

        # np.testing.assert_equal(self._u.static_dielectric_tensor,
        #                         expected_static)
        # np.testing.assert_equal(self._u.ionic_dielectric_tensor, expected_ionic)
        # np.testing.assert_equal(self._u.total_dielectric_tensor, expected_total)


# class PerfectSupercellTest(unittest.TestCase):

    # def setUp(self):
    #     """ """
    #     self._p = PerfectSupercell.from_vasp_results(TEST_PERFECT_DIRECTORY)
    #     self._p_set = PerfectSupercell()

    # def test_set_vasp_results(self):
    #     self._p_set.set_vasp_results(TEST_PERFECT_DIRECTORY)
    #     self.assertEqual(self._p.final_structure, self._p_set.final_structure)
    #     self.assertEqual(self._p.total_energy, self._p_set.total_energy)
    #     np.testing.assert_equal(self._p.eigenvalues[Spin.up],
    #                             self._p_set.eigenvalues[Spin.up])
    #     self.assertEqual(self._p.electrostatic_potential,
    #                      self._p_set.electrostatic_potential)

    # def test_json(self):
    #     self._p.to_json_file("test_supercell.json")
    #     self._p_json = PerfectSupercell.json_load("test_supercell.json")


# class DefectTest(unittest.TestCase):

    # def setUp(self):

        # self._p = PerfectSupercell.from_vasp_results(TEST_PERFECT_DIRECTORY)
        # self._d = DefectSupercell.from_vasp_results(TEST_VAC_DIRECTORY)
        # defect_entry = DefectEntry.json_load(VAC_DEFECT_ENTRY_JSON)
        # self._d_set = DefectSupercell(defect_entry)

    # def test_set_vasp_results(self):
    #     self._d_set.set_vasp_results(TEST_VAC_DIRECTORY)
    #     print(self._d_set.charge)
    #     print(self._d_set.total_energy)
    #     self._d_set.set_relative_values(self._p)
    #     print(self._d_set.relative_total_energy)

#

# #        self.assertEqual(self._d.final_structure, self._d_set.final_structure)
# #        self.assertEqual(self._d.total_energy, self._d_set.total_energy)
# #        np.testing.assert_equal(self._d.eigenvalues[Spin.up],
# #                                self._d_set.eigenvalues[Spin.up])
# #        self.assertEqual(self._d.electrostatic_potential,
# #                         self._d_set.electrostatic_potential)

# #    def test_set_relative_values(self):


# if __name__ == "__main__":
#     unittest.main()

