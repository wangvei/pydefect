# -*- coding: utf-8 -*-
import os
import tempfile
import numpy as np

from pydefect.corrections.efnv_corrections import (
    calc_max_sphere_radius, create_lattice_set, calc_relative_potential,
    Ewald, ExtendedFnvCorrection, point_charge_energy, calc_ewald_sum,
    constants_for_anisotropic_ewald_sum)
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_entry import DefectEntry
from pydefect.util.testing import PydefectTest


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
        filename = ["defects", "MgO", "Va_O1_2", "dft_results.json"]
        defect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "Va_O1_2", "defect_entry.json"]
        defect_entry = self.get_object_by_name(
            DefectEntry.load_json, filename)

        filename = ["defects", "MgO", "perfect", "dft_results.json"]
        perfect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        self.relative_potential = \
            calc_relative_potential(defect=defect, perfect=perfect,
                                    defect_entry=defect_entry)

    def test(self):
        expected = 63
        self.assertEqual(expected, len(self.relative_potential))

        expected = -34.8800 - -34.6537 # -0.2263
        self.assertAlmostEqual(expected, self.relative_potential[0], 4)


class EwaldTest(PydefectTest):

    def setUp(self):
        filename = (self.TEST_FILES_DIR / "defects" / "MgO" / "unitcell.json")
        unitcell = UnitcellCalcResults.load_json(filename)

        filename = ["defects", "MgO", "perfect", "dft_results.json"]
        perfect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        structure = perfect.final_structure
        dielectric_tensor = unitcell.total_dielectric_tensor

        self.ewald = Ewald.from_optimization(structure, dielectric_tensor,
                                             prod_cutoff_fwhm=25)

    def test_msonable(self):
        self.assertMSONable(self.ewald)

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self.ewald.to_json_file(tmp_file.name)
        ewald_from_json = DefectEntry.load_json(tmp_file.name)
        self.assertEqual(self.ewald.as_dict(), ewald_from_json.as_dict())

    def test_optimize(self):
        expected = 11753
        actual = len(self.ewald.reciprocal_neighbor_lattices)
        self.assertEqual(expected, actual)

        expected = [-109.45293, 33.67782, 25.25837]
        actual = self.ewald.real_neighbor_lattices[100]
        self.assertArrayAlmostEqual(expected, actual, 5)


class ExtendedFnvCorrectionTest(PydefectTest):

    def setUp(self):
        filename = (self.TEST_FILES_DIR / "defects" / "MgO" / "unitcell.json")
        unitcell = UnitcellCalcResults.load_json(filename)

        filename = ["defects", "MgO", "perfect", "dft_results.json"]
        perfect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        structure = perfect.final_structure
        dielectric_tensor = unitcell.total_dielectric_tensor

        ewald = Ewald.from_optimization(structure, dielectric_tensor,
                                        prod_cutoff_fwhm=20)

        filename = ["defects", "MgO", "Va_O1_2", "dft_results.json"]
        defect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "Va_O1_2", "defect_entry.json"]
        defect_entry = self.get_object_by_name(
            DefectEntry.load_json, filename)

        self.correction = \
            ExtendedFnvCorrection.compute_correction(
                defect_entry=defect_entry,
                defect_dft=defect,
                perfect_dft=perfect,
                dielectric_tensor=unitcell.total_dielectric_tensor,
                ewald=ewald)

        self.correction.manually_added_correction_energy = 0.1

    def test_dict(self):
        d = self.correction.as_dict()
        from_dict = ExtendedFnvCorrection.from_dict(d).as_dict()
        self.assertEqual(d, from_dict)

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self.correction.to_json_file(tmp_file.name)
        from_json = ExtendedFnvCorrection.load_json(tmp_file.name).as_dict()
        d = self.correction.as_dict()
        self.assertEqual(d.keys(), from_json.keys())

    def test_compute_extended_fnv(self):
        actual = self.correction.lattice_energy
        expected = -0.9734410724309901
        self.assertAlmostEqual(expected, actual, 3)

        actual = self.correction.ave_pot_diff
        expected = 0.15425136945427362
        self.assertAlmostEqual(expected, actual, 3)

        actual = self.correction.alignment_correction_energy
        expected = -0.30850273890854724
        self.assertAlmostEqual(expected, actual, 3)

        actual = len(self.correction.pc_pot)
        expected = 63
        self.assertEqual(expected, actual)

        actual = self.correction.electrostatic_pot[0]
        expected = 0.2263
        self.assertAlmostEqual(expected, actual, 7)

    def test_manual_energy(self):
        actual = self.correction.manually_added_correction_energy
        expected = 0.1
        self.assertAlmostEqual(expected, actual, 7)

    def test_plot_distance_vs_potential(self):
        actual = self.correction.max_sphere_radius
        expected = 8.419456 / 2
        self.assertAlmostEqual(actual, expected)
        self.correction.plot_potential("pot.pdf")
        os.remove("pot.pdf")

