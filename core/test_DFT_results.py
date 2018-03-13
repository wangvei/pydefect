#!/usr/bin/env python

import unittest
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
DIRNAME_UNITCELL = TEST_DIRECTORY + "/unitcell"


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
        self._MgO_Va_O1_2_fvf = \
            SupercellDftResults.from_vasp_files(DIRNAME_VAC)

        self.d = self._MgO_Va_O1_2.as_dict()
        self.d_fvf = self._MgO_Va_O1_2_fvf.as_dict()

    def test_from_vasp_files(self):

        self.assertTrue(self.d["final_structure"] == self.d_fvf["final_structure"])
        self.assertEqual(self.d["total_energy"], self.d_fvf["total_energy"])
        self.assertTrue(
            (self.d["eigenvalues"]["1"] == self.d_fvf["eigenvalues"]["1"]))
        self.assertTrue(
            self.d["electrostatic_potential"] == self.d_fvf["electrostatic_potential"])

    def test_dict(self):
        MgO_Va_O1_2_fd = SupercellDftResults.from_dict(self.d_fvf)
        self.assertTrue(
            self._MgO_Va_O1_2.final_structure == MgO_Va_O1_2_fd.final_structure)

    def test_json(self):
        self._MgO_Va_O1_2.to_json_file("test_DFT_results.json")

#    def test_json_load(self):
#        self._MgO_Va_O1_2.json_load("test_DFT_results.json")

#class UnitcellDftResultsTest(unittest.TestCase):
#
#    def setUp(self):
#        """ """


if __name__ == "__main__":
    unittest.main()

