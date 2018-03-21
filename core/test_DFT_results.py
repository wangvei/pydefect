#!/usr/bin/env python

import unittest
import numpy as np
from DFT_results import *
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.electronic_structure.core import Spin

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "Feb. 25, 2018"

#FILENAME_TO_JSON_FILE_VAC = "examples/Va_Mg1_-2.json"
#FILENAME_TO_JSON_FILE_INT = "examples/Mg_i1_1.json"

TEST_DIRECTORY = "../examples/MgO"
DIRNAME_VAC = TEST_DIRECTORY + "/defects/Va_O1_2"
DIRNAME_UNITCELL = TEST_DIRECTORY + "/unitcell/structure_optimization"
DIRNAME_DIELECTRIC = TEST_DIRECTORY + "/unitcell/dielectric_constants"


class SupercellDftResultsTest(unittest.TestCase):

    def setUp(self):
        """ """
        # CAUTION: When constructing Structure object from Structure.from_file
        #          structure = Structure.from_file(DIRNAME_VAC + “/CONTCAR”)
        #          velocities are not stored.
        #          Therefore, equality becomes False.
        contcar = Poscar.from_file(DIRNAME_VAC + "/CONTCAR")
        final_structure = contcar.structure
        total_energy = -93.76904720
        eigenvalues = {Spin.up: np.array(
            [[[-14.2806, 1.], [-13.4696, 1.], [-13.1066, 1.], [-12.9398, 1.],
              [-12.9398, 1.], [-12.7681, 1.], [-12.7681, 1.], [-1.4322, 1.],
              [-1.4322, 1.], [-1.2961, 1.], [-0.9877, 1.], [-0.667, 1.],
              [-0.3208, 1.], [-0.3208, 1.], [0.9452, 1.], [1.2223, 1.],
              [1.2223, 1.], [1.4722, 1.], [1.4722, 1.], [1.6674, 1.],
              [1.7079, 1.], [1.8786, 1.], [1.8786, 1.], [2.1577, 1.],
              [2.3723, 1.], [2.3723, 1.], [2.5667, 1.], [2.5667, 1.],
              [4.3061, 0.], [8.9622, 0.], [10.0048, 0.], [10.5871, 0.],
              [10.5871, 0.], [11.374, 0.], [11.606, 0.], [11.606, 0.]],
             [[-14.1444, 1.], [-13.4227, 1.], [-13.2385, 1.], [-12.9832, 1.],
              [-12.9615, 1.], [-12.8634, 1.], [-12.7333, 1.], [-0.9824, 1.],
              [-0.9814, 1.], [-0.8391, 1.], [-0.5031, 1.], [-0.3256, 1.],
              [-0.1488, 1.], [0.1553, 1.], [0.408, 1.], [0.6119, 1.],
              [0.6527, 1.], [1.0225, 1.], [1.2301, 1.], [1.288, 1.],
              [1.5418, 1.], [1.5436, 1.], [1.7448, 1.], [1.9201, 1.],
              [2.0783, 1.], [2.2655, 1.], [2.3217, 1.], [2.4294, 1.],
              [5.3997, 0.], [8.5505, 0.], [9.4856, 0.], [9.9455, 0.],
              [11.049, 0.], [11.9159, 0.], [12.5617, 0.], [12.8315, 0.]]])}

        electrostatic_potential = \
            [-34.59, -35.5244, -35.5244, -35.5244, -35.5244, -35.5244, -35.5244,
             -34.59, -70.0739, -70.0739, -70.0739, -70.0739, -70.0739, -70.0739,
             -70.4981]

        self._MgO_Va_O1_2 = SupercellDftResults(
            final_structure, total_energy, eigenvalues, electrostatic_potential)
        self._MgO_Va_O1_2_from_vasp_files = \
            SupercellDftResults.from_vasp_files(DIRNAME_VAC)

        self.d = self._MgO_Va_O1_2.as_dict()
        self.d_from_vasp_files = self._MgO_Va_O1_2_from_vasp_files.as_dict()

    def test_from_vasp_files(self):

        self.assertTrue(self.d["final_structure"] == self.d_from_vasp_files["final_structure"])
        self.assertEqual(self.d["total_energy"], self.d_from_vasp_files["total_energy"])
        self.assertTrue(
            (self.d["eigenvalues"]["1"] == self.d_from_vasp_files["eigenvalues"]["1"]))
        self.assertTrue(
            self.d["electrostatic_potential"] == self.d_from_vasp_files["electrostatic_potential"])

    def test_dict(self):
        MgO_Va_O1_2_fd = SupercellDftResults.from_dict(self.d_from_vasp_files)
        self.assertTrue(
            self._MgO_Va_O1_2.final_structure == MgO_Va_O1_2_fd.final_structure)

    def test_json(self):
        self._MgO_Va_O1_2.to_json_file("test_DFT_results.json")

    def test_json_load(self):
        self._MgO_Va_O1_2.json_load("test_DFT_results.json")


class UnitcellDftResultsTest(unittest.TestCase):

    def setUp(self):
        """ """
        contcar = Poscar.from_file(DIRNAME_UNITCELL + "/CONTCAR")
        final_structure = contcar.structure
        total_energy = -11.91129199
        self.static_dielectric_tensor = np.array(
            [[3.166727, 0.0, 0.0], [0.0, 3.166727, 0.0], [0.0, 0.0, 3.166727]])
        self.ionic_dielectric_tensor = np.array(
            [[9.102401, -0.0, -0.0],
             [-0.0, 9.102448, 0.0],
             [-0.0, 0.0, 9.102542]])
        eigenvalues = {Spin.up: np.array(
            [[[-14.1364, 1], [2.9978, 1], [2.9978, 1], [2.9978, 1], [7.479, 0],
              [18.6525, 0], [18.6525, 0], [18.6525, 0]],
             [[-13.9805, 1], [2.0526, 1], [2.8993, 1], [2.8993, 1], [8.3965, 0],
              [16.9813, 0], [18.9958, 0], [18.9958, 0]],
             [[-13.5759, 1], [0.3868, 1], [2.6721, 1], [2.6721, 1], [9.8536, 0],
              [15.4092, 0], [19.4381, 0], [19.4381, 0]],
             [[-13.1221, 1], [-0.9774, 1], [2.4581, 1], [2.4581, 1],
              [10.6136, 0], [15.1306, 0], [19.2233, 0], [19.2233, 0]],
             [[-12.9153, 1], [-1.5075, 1], [2.3728, 1], [2.3728, 1],
              [10.7333, 0], [15.3147, 0], [19.0457, 0], [19.0457, 0]],
             [[-13.9306, 1], [1.9735, 1], [2.7729, 1], [2.7729, 1], [8.6915, 0],
              [16.5049, 0], [19.1306, 0], [19.1306, 0]],
             [[-13.6143, 1], [0.7711, 1], [2.5295, 1], [2.5332, 1], [9.9851, 0],
              [15.1627, 0], [19.1168, 0], [19.5995, 0]],
             [[-13.1791, 1], [-0.517, 1], [2.175, 1], [2.327, 1], [11.0794, 0],
              [14.7426, 0], [18.1219, 0], [19.6339, 0]],
             [[-12.8837, 1], [-1.247, 1], [1.8888, 1], [2.2631, 1],
              [11.4044, 0], [15.1238, 0], [17.611, 0], [19.4261, 0]],
             [[-12.9493, 1], [-1.0815, 1], [1.8251, 1], [2.3726, 1],
              [11.3855, 0], [15.0655, 0], [17.7046, 0], [19.496, 0]],
             [[-13.3232, 1], [-0.1161, 1], [2.0394, 1], [2.6021, 1],
              [10.8135, 0], [14.919, 0], [18.3212, 0], [19.99, 0]],
             [[-13.743, 1], [1.1538, 1], [2.4897, 1], [2.8252, 1], [9.473, 0],
              [15.7733, 0], [18.9498, 0], [19.6815, 0]],
             [[-13.4202, 1], [0.4155, 1], [2.2821, 1], [2.2821, 1], [10.808, 0],
              [14.0963, 0], [18.1994, 0], [20.1894, 0]],
             [[-13.0858, 1], [-0.3016, 1], [1.8392, 1], [2.0718, 1],
              [11.9148, 0], [13.5683, 0], [16.8649, 0], [20.4933, 0]],
             [[-12.8004, 1], [-0.725, 1], [1.1377, 1], [2.0118, 1],
              [12.5168, 0], [14.3505, 0], [16.0304, 0], [20.3257, 0]],
             [[-12.7554, 1], [-0.6608, 1], [0.6861, 1], [2.1312, 1],
              [12.7568, 0], [14.9282, 0], [15.8278, 0], [20.1413, 0]],
             [[-12.9825, 1], [-0.4211, 1], [1.0582, 1], [2.3723, 1],
              [12.2422, 0], [14.6272, 0], [16.7207, 0], [20.2032, 0]],
             [[-12.8905, 1], [-0.635, 1], [1.8501, 1], [1.8501, 1],
              [12.4067, 0], [12.5923, 0], [16.1254, 0], [21.1754, 0]],
             [[-12.7076, 1], [-0.7827, 1], [1.3214, 1], [1.7777, 1],
              [12.5145, 0], [13.1844, 0], [15.7109, 0], [21.2011, 0]],
             [[-12.6521, 1], [-0.4035, 1], [0.3589, 1], [1.8889, 1],
              [12.9598, 0], [14.4795, 0], [15.2963, 0], [20.7678, 0]],
             [[-12.6666, 1], [-0.9817, 1], [1.685, 1], [1.685, 1], [12.0862, 0],
              [12.6907, 0], [15.8269, 0], [21.5828, 0]],
             [[-13.2865, 1], [0.042, 1], [1.8809, 1], [2.4059, 1], [11.236, 0],
              [14.378, 0], [17.7916, 0], [20.2447, 0]],
             [[-12.9293, 1], [-0.6225, 1], [1.322, 1], [2.1614, 1],
              [12.2754, 0], [14.4725, 0], [16.6366, 0], [19.8304, 0]],
             [[-12.7711, 1], [-0.9032, 1], [1.1201, 1], [2.0452, 1],
              [12.4892, 0], [14.9483, 0], [16.2195, 0], [19.5838, 0]],
             [[-12.8293, 1], [-0.5655, 1], [1.3524, 1], [1.7898, 1],
              [13.0016, 0], [13.3066, 0], [15.9116, 0], [20.1318, 0]],
             [[-12.6711, 1], [-0.4187, 1], [0.5802, 1], [1.6609, 1],
              [13.4357, 0], [14.4675, 0], [15.2594, 0], [19.2653, 0]],
             [[-12.715, 1], [-0.3017, 1], [0.3761, 1], [1.8512, 1],
              [13.5396, 0], [14.5505, 0], [15.4167, 0], [19.3251, 0]],
             [[-12.642, 1], [-0.6877, 1], [1.1518, 1], [1.4743, 1],
              [12.9233, 0], [13.563, 0], [15.5076, 0], [19.9636, 0]],
             [[-12.6177, 1], [0.1231, 1], [0.1231, 1], [1.2747, 1],
              [14.4722, 0], [14.4722, 0], [14.8549, 0], [18.1715, 0]]])}

        self._MgO_unitcell = UnitcellDftResults(
            final_structure, total_energy, static_dielectric_tensor=None,
            ionic_dielectric_tensor=None, eigenvalues=eigenvalues)

        self._MgO_unitcell_from_vasp_files = \
            UnitcellDftResults.from_vasp_files(DIRNAME_UNITCELL)

        # go through a dict
        self.d = self._MgO_unitcell.as_dict()
        self.d_from_vasp_files = self._MgO_unitcell_from_vasp_files.as_dict()

    def test_from_vasp_files(self):
        self.assertEqual(self.d["final_structure"],
                         self.d_from_vasp_files["final_structure"])
        self.assertEqual(self.d["total_energy"],
                         self.d_from_vasp_files["total_energy"])
        self.assertEqual(self.d["eigenvalues"]["1"],
                         self.d_from_vasp_files["eigenvalues"]["1"])

    def test_set_dielectric_constants_from_outcar(self):
        self._MgO_unitcell.\
            set_dielectric_constants_from_outcar(DIRNAME_DIELECTRIC)
        np.testing.assert_equal(self._MgO_unitcell.static_dielectric_tensor,
                                self.static_dielectric_tensor)

    def test_set_static_dielectric_tensor(self):
        a = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self._MgO_unitcell.static_dielectric_tensor = a
        np.testing.assert_equal(self._MgO_unitcell.static_dielectric_tensor, a)

    def test_set_ionic_dielectric_tensor(self):
        b = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        self._MgO_unitcell.ionic_dielectric_tensor = b
        np.testing.assert_equal(self._MgO_unitcell.ionic_dielectric_tensor, b)


if __name__ == "__main__":
    unittest.main()

