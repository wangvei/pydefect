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
        unitcell_dir = (self.TEST_FILES_DIR / "core" / "MgO" / "unitcell"
                        / "dielectric_constants")
        unitcell.set_static_dielectric_tensor_from_vasp(unitcell_dir)
        unitcell.set_ionic_dielectric_tensor_from_vasp(unitcell_dir)
        perfect_dir = \
            self.TEST_FILES_DIR / "core" / "MgO" / "defects" / "perfect"
        perfect = SupercellCalcResults.from_vasp_files(perfect_dir)
        structure = perfect.final_structure
        dielectric_tensor = unitcell.total_dielectric_tensor
        self.ewald = Ewald.from_optimization(structure, dielectric_tensor,
                                             prod_cutoff_fwhm=25)

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self.ewald.to_json_file(tmp_file.name)
        ewald_from_json = DefectEntry.load_json(tmp_file.name)
        self.assertEqual(self.ewald.as_dict(), ewald_from_json.as_dict())

    def test_optimize(self):
        self.assertEqual(10693, len(self.ewald.reciprocal_neighbor_lattices))
        self.assertArrayAlmostEqual([16.987572, -63.703395, -29.728251],
                                    self.ewald.real_neighbor_lattices[100])


class ExtendedFnvCorrectionTest(PydefectTest):

    def setUp(self):
        unitcell = UnitcellCalcResults()
        unitcell_dir = (self.TEST_FILES_DIR / "core" / "MgO" / "unitcell"
                        / "dielectric_constants")
        unitcell.set_static_dielectric_tensor_from_vasp(unitcell_dir)
        unitcell.set_ionic_dielectric_tensor_from_vasp(unitcell_dir)
        perfect_dir = \
            self.TEST_FILES_DIR / "core" / "MgO" / "defects" / "perfect"
        perfect = SupercellCalcResults.from_vasp_files(perfect_dir)
        structure = perfect.final_structure
        dielectric_tensor = unitcell.total_dielectric_tensor
        ewald = Ewald.from_optimization(structure, dielectric_tensor,
                                        prod_cutoff_fwhm=25)
        tmp_file = tempfile.NamedTemporaryFile()
        ewald.to_json_file(tmp_file.name)

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
        self.vacancy_correction = \
            ExtendedFnvCorrection.compute_correction(
                defect_entry=defect_entry,
                defect_dft=defect,
                perfect_dft=perfect,
                unitcell_dft=unitcell,
                ewald_json=tmp_file.name)

    #        expected_vacancy_lattice_energy = -1.2670479

    def test_dict(self):
        d = self.vacancy_correction.as_dict()
        from_dict = ExtendedFnvCorrection.from_dict(d).as_dict()
        print(d)
        print(from_dict)
        self.assertEqual(d, from_dict)

    # def test_json(self):
    #     tmp_file = tempfile.NamedTemporaryFile()
    #     self.vacancy_correction.to_json_file(tmp_file.name)
    #     from_json = ExtendedFnvCorrection.load_json(tmp_file.name).as_dict()
    #     d = self.vacancy_correction.as_dict()
    #     print(d)
    #     print(from_json)
    #     self.assertEqual(d, from_json)

    def test_compute_extended_fnv(self):
        # vacancy
        expected_vacancy_lattice_energy = -1.2670479
        self.assertAlmostEqual(self.vacancy_correction.lattice_energy,
                               expected_vacancy_lattice_energy, 3)
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

        expected_vacancy_model_pot = [
            -0.221616107852, -0.0764266164062, -0.0752496904232,
            -0.0776770615374, -0.0761845615171, -0.0791754509674,
            -0.0798985084573, -0.221615505252, -0.160981719771,
            -0.16096545335, -0.160981815887, -0.160975167555,
            -0.160983850502, -0.160982833121, -0.30114977165]
        expected_vacancy_difference_electrostatic_pot = [
            -0.7019, 0.2297, 0.2336, 0.2350, 0.2338, 0.2323, 0.2259, -0.7017,
            0.2582, 0.2575, 0.2573, 0.2580, 0.2581, 0.2584, 0.6767 ]
        expected_vacancy_distances_list = [
            3.67147, 2.25636, 2.25288, 2.26013, 2.25573, 2.26446, 2.26664,
            3.67347, 2.99998, 2.99316, 2.99870, 2.99665, 3.00100, 2.99985,
            4.24097]

        self.assertArrayAlmostEqual(self.vacancy_correction.model_pot,
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
        expected_max_sphere_radius = 2.45194
        self.assertAlmostEqual(self.vacancy_correction.max_sphere_radius,
                               expected_max_sphere_radius, 5)
        self.vacancy_correction.plot_distance_vs_potential("pot.pdf")

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

