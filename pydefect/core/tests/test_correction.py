import unittest
import os

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pydefect.core.correction import Ewald, Correction, CorrectionMethod
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 20, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core", "MgO")


dirname_unitcell = test_dir + "/unitcell/structure_optimization"
dirname_dielectric = test_dir + "/unitcell/dielectric_constants"

dirname_perfect = test_dir + "/defects/perfect"
# interstitial directory

vacancy_expected_alignment = "NEED_TO_SPECIFY!"
dirname_vacancy = test_dir + "/defects/Va_O1_2/"
vac_defect_entry_json = dirname_vacancy + "/defect_entry.json"

# calculated by shell script made by Dr. Kumagai
expected_ewald = 0.0411503184705
expected_num_real_vector = 11051
expected_num_reciprocal_vector = 10608
expected_vacancy_potential_difference = 0.2812066
expected_vacancy_alignment_like_term = -0.5624132
expected_vacancy_lattice_energy = -1.2670479
expected_vacancy_symbols = ["Mg"] * 8 + ["O"] * 7
expected_vacancy_model_pot = [-0.221616953329,
                              -0.0757569991956,
                              -0.0745712827233,
                              -0.0770321065632,
                              -0.0755424109893,
                              -0.0784910684584,
                              -0.0792252195274,
                              -0.22161641786,
                              -0.160982214428,
                              -0.160966818748,
                              -0.16098264708,
                              -0.160975656079,
                              -0.160984578901,
                              -0.160983417362,
                              -0.30115082168]
expected_vacancy_difference_electrostatic_pot = [
    -0.7019,
    0.2297,
    0.2336,
    0.2350,
    0.2338,
    0.2323,
    0.2259,
    -0.7017,
    0.2582,
    0.2575,
    0.2573,
    0.2580,
    0.2581,
    0.2584,
    0.6767
]
expected_vacancy_distances_list = [3.67147,
                                   2.25636,
                                   2.25288,
                                   2.26013,
                                   2.25573,
                                   2.26446,
                                   2.26664,
                                   3.67347,
                                   2.99998,
                                   2.99316,
                                   2.99870,
                                   2.99665,
                                   3.00100,
                                   2.99985,
                                   4.24097]


class EwaldTest(unittest.TestCase):

    def setUp(self):
        unitcell = UnitcellDftResults()
        unitcell.set_static_dielectric_tensor_from_vasp(dirname_dielectric)
        unitcell.set_ionic_dielectric_tensor_from_vasp(dirname_dielectric)
        perfect = SupercellDftResults.from_vasp_files(dirname_perfect)
        self._structure = perfect.final_structure
        self._dielectric_tensor = unitcell.total_dielectric_tensor

    def test_optimize(self):
        ewald =\
            Ewald.from_optimization(self._structure, self._dielectric_tensor)

        self.assertAlmostEqual(0, 1-expected_ewald/ewald.ewald_param, 1)
        self.assertAlmostEqual(
            0, 1-expected_num_real_vector/ewald.num_real_lattice, 1)
        self.assertAlmostEqual(
            0, 1-expected_num_reciprocal_vector/ewald.num_reciprocal_lattice, 1)
        # print(ewald.ewald_param)
        # print(ewald.num_real_lattice)
        # print(ewald.num_reciprocal_lattice)
        # print(ewald.dielectric_tensor)
        # print(ewald.num_real_lattice)


class CorrectionTest(unittest.TestCase):

    def setUp(self):
        self._unitcell = UnitcellDftResults()
        self._unitcell.set_static_dielectric_tensor_from_vasp(
            dirname_dielectric)
        self._unitcell.set_ionic_dielectric_tensor_from_vasp(dirname_dielectric)
        self._perfect = SupercellDftResults.from_vasp_files(dirname_perfect)
        self._vacancy_entry = DefectEntry.json_load(vac_defect_entry_json)
        self._vacancy = SupercellDftResults.from_vasp_files(dirname_vacancy)
        self._structure = self._perfect.final_structure
        self._dielectric_tensor = self._unitcell.total_dielectric_tensor
        self._ewald = Ewald(
            lattice_matrix=self._structure.lattice.matrix,
            dielectric_tensor=self._unitcell.total_dielectric_tensor,
            ewald_param=expected_ewald,
            num_real_lattice=expected_num_real_vector,
            num_reciprocal_lattice=expected_num_reciprocal_vector)
        # self._ewald = \
        #     Ewald.from_optimization(self._structure, self._dielectric_tensor)
        print("Correction_test setUp completed")

    def test_compute_extended_fnv(self):
        vacancy_correction = \
            Correction.compute_correction(self._vacancy_entry,
                                          self._vacancy,
                                          self._perfect,
                                          self._unitcell)
        self.assertAlmostEqual(vacancy_correction.lattice_energy,
                               expected_vacancy_lattice_energy, 4)
        self.assertAlmostEqual(vacancy_correction.diff_ave_pot,
                               expected_vacancy_potential_difference, 5)
        self.assertAlmostEqual(vacancy_correction.alignment,
                               expected_vacancy_alignment_like_term, 5)

        # TODO: Check case of irreducible sites like O1, O2
        assert_array_equal(vacancy_correction.symbols_without_defect,
                           expected_vacancy_symbols)
        assert_array_almost_equal(vacancy_correction.distances_from_defect,
                                  expected_vacancy_distances_list, 5)
        assert_array_almost_equal(vacancy_correction.model_pot,
                                  expected_vacancy_model_pot, 5)
        assert_array_almost_equal(vacancy_correction.difference_electrostatic_pot,
                                  expected_vacancy_difference_electrostatic_pot, 5)

    def test_plot_distance_vs_potential(self):

        vacancy_correction = \
            Correction(CorrectionMethod.extended_fnv,
                       self._ewald,
                       expected_vacancy_lattice_energy,
                       expected_vacancy_potential_difference,
                       expected_vacancy_alignment_like_term,
                       expected_vacancy_symbols,
                       expected_vacancy_distances_list,
                       expected_vacancy_difference_electrostatic_pot,
                       expected_vacancy_model_pot)
        expected_max_sphere_radius = 2.45194
        self.assertAlmostEqual(vacancy_correction.max_sphere_radius,
                               expected_max_sphere_radius, 5)
        vacancy_correction.plot_distance_vs_potential()


if __name__ == "__main__":
    unittest.main()
