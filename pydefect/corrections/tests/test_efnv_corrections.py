import unittest
import os
import tempfile
import numpy as np
from pydefect.util.testing import PydefectTest

from pydefect.corrections.efnv_corrections import *
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class CalcMaxSphereRadiusTest(PydefectTest):

    def setUp(self) -> None:
        lattice_vectors_1 = np.array([[5, 0, 0], [0, 10, 0], [0, 0, 20]])
        lattice_vectors_2 = np.array([[10, 0, 0], [0, 10, 0], [10, 10, 10]])
        self.radius_1 = calc_max_sphere_radius(lattice_vectors_1)
        self.radius_2 = calc_max_sphere_radius(lattice_vectors_2)

    def test_radii(self):
        self.assertEqual(10.0, self.radius_1)
        self.assertEqual(5.0, self.radius_2)


class CreateLatticeSetTest(PydefectTest):

    def setUp(self) -> None:
        lattice_vectors_1 = np.array([[10, 10, 0], [10, -10, 0], [0, 0, 14]])
        max_length = 15
        self.set = create_lattice_set(lattice_vectors_1, max_length)

    def test_lattice_set(self):
        expected = [[-10, -10, 0], [-10, 10, 0], [0, 0, -14], [0, 0, 0],
                    [0, 0, 14], [10, -10, 0], [10, 10, 0]]
        self.assertEqual(expected, self.set)


class CalcRelativePotentialTest(PydefectTest):

    def setUp(self) -> None:
        dirname = os.path.join("core", "MgO", "defects")
        defect_name = os.path.join(dirname, "Va_O1_2", "dft_results.json")
        defect = self.get_object_by_name(SupercellCalcResults.load_json,
                                         defect_name)
        entry_name = os.path.join(dirname, "Va_O1_2", "defect_entry.json")
        defect_entry = self.get_object_by_name(DefectEntry.load_json,
                                               entry_name)
        perfect_name = os.path.join(dirname, "perfect", "dft_results.json")
        perfect = self.get_object_by_name(SupercellCalcResults.load_json,
                                          perfect_name)

        self.relative_potential = \
            calc_relative_potential(defect=defect, perfect=perfect,
                                    defect_entry=defect_entry)

    def test(self):
        expected = [0.7231, -0.2318, -0.2361, -0.2374, -0.2363, -0.2343,
                    -0.2282, 0.7229, -0.2351, -0.2341, -0.2334, -0.2355,
                    -0.2352, -0.2353, -0.6934]
        self.assertArrayAlmostEqual(expected, self.relative_potential, 4)


class EwaldTest(PydefectTest):

    def setUp(self):
        unitcell = UnitcellCalcResults()
        unitcell_dir = self.TEST_FILES_DIR / "core" / "MgO" / "unitcell" / "dielectric_constants"
        unitcell.set_static_dielectric_tensor_from_vasp(unitcell_dir)
        unitcell.set_ionic_dielectric_tensor_from_vasp(unitcell_dir)
        perfect_dir = \
            self.TEST_FILES_DIR / "core" / "MgO" / "defects" / "perfect"
        perfect = SupercellCalcResults.from_vasp_files(perfect_dir)
        structure = perfect.final_structure
        dielectric_tensor = unitcell.total_dielectric_tensor
        ewald = 0.1
        prod_cutoff_fwhm = 0.1
        self.ewald = Ewald.from_optimization(structure, dielectric_tensor,
                                             prod_cutoff_fwhm=10)

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.ewald.to_json_file(tmp_file.name)
        ewald_from_json = DefectEntry.load_json(tmp_file.name)
        self.assertEqual(self.ewald.as_dict(), ewald_from_json.as_dict())

    def test_optimize(self):
        print(len(self.ewald.reciprocal_neighbor_lattices))
        print(self.ewald.real_neighbor_lattices)

        # expected_num_reciprocal_vector = 10693
        # self.assertEqual(len(ewald.reciprocal_neighbor_lattices),
        #                  expected_num_reciprocal_vector)
        # expected = np.array([16.987572, -63.703395, -29.728251])
        # self.assertArrayAlmostEqual(expected,
        #                             ewald.real_neighbor_lattices[100])
        # print(ewald.real_neighbor_lattices[100])

#        ewald.to_json_file("ewald.json")


class ExtendedFnvCorrectionTest(PydefectTest):

    def setUp(self):
        self._unitcell = UnitcellCalcResults()
        self._unitcell.set_static_dielectric_tensor_from_vasp(
            dirname_dielectric)
        self._unitcell.set_ionic_dielectric_tensor_from_vasp(dirname_dielectric)
        self._perfect = SupercellCalcResults.from_vasp_files(dirname_perfect)
        self._structure = self._perfect.final_structure
        self._dielectric_tensor = self._unitcell.total_dielectric_tensor

        self._lattice = self._perfect.final_structure.lattice
        self._dielectric_tensor = self._unitcell.total_dielectric_tensor
        _ewald = 0.1
        _prod_cutoff_fwhm = 0.1
        _real_neighbor_lattices = [[1, 1, 1], [2, 2, 2]]
        _reciprocal_neighbor_lattices = [[1, 1, 1], [2, 2, 2]]
        Ewald(self._lattice, self._dielectric_tensor, _ewald, _prod_cutoff_fwhm, _real_neighbor_lattices, _reciprocal_neighbor_lattices).to_json_file("ewald.json")
        self._ewald = "ewald.json"

        # self.ewald = \
        #     Ewald.from_optimization(self.perfect_structure, self.dielectric_tensor)
        self._vacancy_entry = DefectEntry.load_json(vac_defect_entry_json)
        self._vacancy = SupercellCalcResults.load_json(vac_dft_results_json)
#        self._interstitial_entry = DefectEntry.load_json(int_defect_entry_json)
#        self._interstitial = \
#            SupercellCalcResults.from_vasp_files(dirname_interstitial)
#        self._substitutional_entry = DefectEntry.load_json(sub_defect_entry_json)
#        self._substitutional = \
#            SupercellCalcResults.from_vasp_files(dirname_substitutional)

    def test_dict(self):
        expected_vacancy_lattice_energy = -1.2670479
        vacancy_correction = \
            ExtendedFnvCorrection(self._ewald,
                                  lattice_matrix=self._lattice.matrix,
                                  lattice_energy=expected_vacancy_lattice_energy,
                                  ave_pot_diff=expected_vacancy_potential_difference,
                                  alignment_correction_energy=expected_vacancy_alignment_like_term,
                                  symbols_without_defect=expected_vacancy_symbols,
                                  distances_from_defect=expected_vacancy_distances_list,
                                  difference_electrostatic_pot=expected_vacancy_difference_electrostatic_pot,
                                  model_pot=expected_vacancy_model_pot,
                                  manual_correction_energy=0.0)
        # object -> dict -> object
        d = vacancy_correction.as_dict()
        correction_from_dict = ExtendedFnvCorrection.from_dict(d)
        self.assertTrue(correction_from_dict.as_dict() == vacancy_correction.as_dict())

    # def test_json(self):
    #     tmp_file = tempfile.NamedTemporaryFile()
    #     vacancy_correction = \
    #         ExtendedFnvCorrection(self._ewald,
    #                               lattice_matrix=self._lattice.matrix,
    #                               lattice_energy=expected_vacancy_lattice_energy,
    #                               ave_pot_diff=expected_vacancy_potential_difference,
    #                               alignment_correction_energy=expected_vacancy_alignment_like_term,
    #                               symbols_without_defect=expected_vacancy_symbols,
    #                               distances_from_defect=expected_vacancy_distances_list,
    #                               difference_electrostatic_pot=expected_vacancy_difference_electrostatic_pot,
    #                               model_pot=expected_vacancy_model_pot,
    #                               manual_correction_energy=0.0)

        # vacancy_correction.to_json_file(tmp_file.name)
        # loaded = ExtendedFnvCorrection.load_json(tmp_file.name)
        # print(vacancy_correction.as_dict())
        # print(loaded.as_dict())

        # self.assertTrue(loaded.as_dict() == vacancy_correction.as_dict())

    def test_compute_extended_fnv(self):
        # vacancy
        vacancy_correction = \
            ExtendedFnvCorrection.compute_correction(self._vacancy_entry,
                                                     self._vacancy,
                                                     self._perfect,
                                                     self._unitcell,
                                                     self._ewald)
        expected_vacancy_lattice_energy = -1.2670479
        self.assertAlmostEqual(vacancy_correction.lattice_energy,
                               expected_vacancy_lattice_energy, 4)
        # self.assertAlmostEqual(vacancy_correction.ave_pot_diff,
        #                        expected_vacancy_potential_difference, 5)
#        self.assertAlmostEqual(vacancy_correction.alignment_correction_energy,
#                               expected_vacancy_alignment_like_term, 5)
        # # TODO: Check case of irreducible sites like O1, O2
        # assert_array_equal(vacancy_correction.symbols_without_defect,
        #                    expected_vacancy_symbols)
        # assert_array_almost_equal(vacancy_correction.distances_from_defect,
        #                           expected_vacancy_distances_list, 5)
#        print(vacancy_correction.model_pot)
        self.assertArrayAlmostEqual(vacancy_correction.model_pot,
                                    expected_vacancy_model_pot, 5)
        # assert_array_almost_equal(vacancy_correction.difference_electrostatic_pot,
        #                           expected_vacancy_difference_electrostatic_pot, 5)

        # # interstitial
        # interstitial_correction = \
        #     ExtendedFnvCorrection.compute_correction(self._interstitial_entry,
        #                                              self._interstitial,
        #                                              self._perfect,
        #                                              self._unitcell)
        # self.assertAlmostEqual(interstitial_correction.lattice_energy,
        #                        expected_interstitial_lattice_energy, 4)
        # self.assertAlmostEqual(interstitial_correction.ave_pot_diff,
        #                        expected_interstitial_potential_difference, 5)
        # self.assertAlmostEqual(interstitial_correction.alignment_energy,
        #                        expected_interstitial_alignment_like_term, 5)
        # # TODO: Check case of irreducible sites like O1, O2
        # assert_array_equal(interstitial_correction.symbols_without_defect,
        #                    expected_interstitial_symbols)
        # assert_array_almost_equal(interstitial_correction.distances_from_defect,
        #                           expected_interstitial_distances_list, 5)
        # assert_array_almost_equal(interstitial_correction.model_pot,
        #                           expected_interstitial_model_pot, 5)
        # assert_array_almost_equal(
        #     interstitial_correction.difference_electrostatic_pot,
        #     expected_interstitial_difference_electrostatic_pot, 5)

        # # substitutional
        # substitutional_correction = \
        #     ExtendedFnvCorrection.compute_correction(self._substitutional_entry,
        #                                              self._substitutional,
        #                                              self._perfect,
        #                                              self._unitcell)
        # self.assertAlmostEqual(substitutional_correction.lattice_energy,
        #                        expected_substitutional_lattice_energy, 4)
        # self.assertAlmostEqual(substitutional_correction.ave_pot_diff,
        #                        expected_substitutional_potential_difference, 5)
        # self.assertAlmostEqual(substitutional_correction.alignment_energy,
        #                        expected_substitutional_alignment_like_term, 5)
        # # TODO: Check case of irreducible sites like O1, O2
        # # Methods to get defect position slightly differs between that of shell
        # # script and pydefect module. Then, distances_from_defect also differ
        # # between shell and pydefect. However, as tested above, difference of
        # # correction energy is not critical at all.
        # assert_array_equal(substitutional_correction.symbols_without_defect,
        #                    expected_substitutional_symbols)
        # assert_array_almost_equal(
        #     substitutional_correction.distances_from_defect,
        #     expected_substitutional_distances_list, 1)
        # assert_array_almost_equal(substitutional_correction.model_pot,
        #                           expected_substitutional_model_pot, 2)
        # assert_array_almost_equal(
        #     substitutional_correction.difference_electrostatic_pot,
        #     expected_substitutional_difference_electrostatic_pot, 2)

    def test_plot_distance_vs_potential(self):

        expected_vacancy_lattice_energy = -1.2670479
        vacancy_correction = \
            ExtendedFnvCorrection(self._ewald,
                                  self._lattice.matrix,
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
        vacancy_correction.plot_distance_vs_potential("pot.pdf")

        # interstitial_correction = \
        #     ExtendedFnvCorrection(self._ewald,
        #                           expected_interstitial_lattice_energy,
        #                           expected_interstitial_potential_difference,
        #                           expected_interstitial_alignment_like_term,
        #                           expected_interstitial_symbols,
        #                           expected_interstitial_distances_list,
        #                           expected_interstitial_difference_electrostatic_pot,
        #                           expected_interstitial_model_pot)
        # interstitial_correction.plot_distance_vs_potential()

        # substitutional_correction = \
        #     ExtendedFnvCorrection(self._ewald,
        #                           expected_substitutional_lattice_energy,
        #                           expected_substitutional_potential_difference,
        #                           expected_substitutional_alignment_like_term,
        #                           expected_substitutional_symbols,
        #                           expected_substitutional_distances_list,
        #                           expected_substitutional_difference_electrostatic_pot,
        #                           expected_substitutional_model_pot)
        # substitutional_correction.plot_distance_vs_potential()


if __name__ == "__main__":
    unittest.main()
