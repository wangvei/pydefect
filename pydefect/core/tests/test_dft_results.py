# -*- coding: utf-8 -*-

import numpy as np
import os
import tempfile
import unittest

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.dft_results import defect_center, distances_from_point, \
    SupercellDftResults, UnitcellDftResults

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "Feb. 25, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectCenterTest(unittest.TestCase):
    def test_defect_center(self):

        name = "Va_O1"
        test_structure = Structure.from_file(
           os.path.join(test_dir, "POSCAR-MgO64atoms-Va_O1"))
        removed_atoms = {32: [0.25, 0.25, 0.25]}
        inserted_atoms = []
        element_diff = {"O": -1}
        charge = 2
        # Construct DefectEntry class object.
        vac_defect_entry = DefectEntry(name, test_structure, removed_atoms,
                                       inserted_atoms, element_diff, charge)

        # Test a single defect.
        actual1 = defect_center(vac_defect_entry)
        expected = [0.25, 0.25, 0.25]
        self.assertEqual(actual1, expected)

        # Test with relaxed structure.
        contcar = Structure.from_file(
            os.path.join(test_dir, "CONTCAR-MgO64atoms-Va_O1"))

        actual2 = defect_center(vac_defect_entry, contcar)
        self.assertEqual(actual2, expected)

        name_complex = "Va_O1+i_2H"
        complex_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO64atoms-Va_O1+i_2H"))
        removed_atoms = {32: [0.25, 0.25, 0.25]}
        inserted_atoms = [63, 64]
        element_diff = {"O": -1, "H": 2}
        charge = 0

        # Construct DefectEntry class object.
        complex_defect_entry = DefectEntry(name_complex, complex_structure,
                                           removed_atoms, inserted_atoms,
                                           element_diff, charge)
        actual3 = defect_center(complex_defect_entry)
        expected = [0.26, 0.26, 0.26]
        self.assertEqual(actual3, expected)


class DistancesFromPointTest(unittest.TestCase):
    def test_distances_from_point(self):
        name = "tmp"
        structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-min_distance_under_pbc"))
        removed_atoms = {}
        inserted_atoms = [0]
        element_diff = {"Fe": 1}
        charge = 0
        # Construct DefectEntry class object.
        defect_entry = DefectEntry(name, structure, removed_atoms,
                                   inserted_atoms, element_diff, charge)

        # (Fe1 - Co1) = 7.70552(0)
        # Fe1 0.25000 0.25000 0.25000(0, 0, 0) + x, y, z
        # Co1 0.50000 -0.12500 -0.12500(0, -1, -1) + x, y, z
        # (Fe1 - Ni1) = 4.67707(0)
        # Fe1  0.25000 0.25000 0.25000(0, 0, 0) + x, y, z
        # Ni1 -0.12500 0.25000 0.50000(-1, 0, 0) + x, y, z
        actual = distances_from_point(structure, defect_entry)
        expected = [0.0, 7.705517503711221, 4.677071733467427]
        self.assertEqual(actual, expected)


class SupercellDftResultsTest(unittest.TestCase):

    def setUp(self):
        """ """
        final_structure = Structure.from_file(
            os.path.join(test_dir, "MgO/defects/Va_O1_2", "CONTCAR"))
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
        # electrostatic_potential is a property because it is used for
        # test_relative_potential method.
        self.electrostatic_potential = \
            [-34.69, -35.5244, -35.5244, -35.5244, -35.5244, -35.5244, -35.5244,
             -34.59, -70.0739, -70.0739, -70.0739, -70.0739, -70.0739, -70.0739,
             -70.4981]

        self._MgO_Va_O1_2 = SupercellDftResults(
            final_structure, total_energy, eigenvalues,
            self.electrostatic_potential)

        self._MgO_Va_O1_2_from_vasp_files = \
            SupercellDftResults.from_vasp_files(
                os.path.join(test_dir, "MgO/defects/Va_O1_2"))

        self._MgO_perfect_from_vasp_files = \
            SupercellDftResults.from_vasp_files(
                os.path.join(test_dir, "MgO/defects/perfect"))

        self.d = self._MgO_Va_O1_2.as_dict()
        self.d_from_vasp_files = self._MgO_Va_O1_2_from_vasp_files.as_dict()

        name = "Va_O1"
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "MgO/defects/Va_O1_2", "POSCAR"))
        removed_atoms = {8: [0.25, 0.25, 0.25]}
        inserted_atoms = {}
        changes_of_num_elements = {"O": -1}
        charge = 2

        self._defect_entry_MgO_Va_O1_2 = \
            DefectEntry(name, initial_structure, removed_atoms, inserted_atoms,
                        changes_of_num_elements, charge)

    def test_from_vasp_files(self):
        # CAUTION: When constructing Structure object from Structure.from_file
        #          velocities are not stored.
        #          Therefore, equality check of Structure objects returns False.
        #          If the structure is converted via poscar file format, it may
        #          be solved.
        # contcar = Poscar.from_file(os.path.join(DIRNAME_UNITCELL, "CONTCAR"))
        # final_structure = contcar.structure

        # self.assertTrue(self.d["initial_structure"] ==
        #                 self.d_from_vasp_files["initial_structure"])
        # self.assertTrue(self.d["final_structure"] ==
        #                 self.d_from_vasp_files["final_structure"])
        self.assertEqual(self.d["total_energy"],
                         self.d_from_vasp_files["total_energy"])
        self.assertTrue((self.d["eigenvalues"]["1"] ==
                         self.d_from_vasp_files["eigenvalues"]["1"]))
        self.assertTrue(self.d["electrostatic_potential"] ==
                        self.d_from_vasp_files["electrostatic_potential"])

    def test_dict(self):
        MgO_Va_O1_2_fd = SupercellDftResults.from_dict(self.d_from_vasp_files)
        np.testing.assert_equal(MgO_Va_O1_2_fd.eigenvalues[Spin.up],
                                self._MgO_Va_O1_2.eigenvalues[Spin.up])

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self._MgO_Va_O1_2.to_json_file(tmp_file.name)
        MgO_Va_O1_2_from_json = SupercellDftResults.json_load(tmp_file.name)
        np.testing.assert_equal(MgO_Va_O1_2_from_json.eigenvalues[Spin.up],
                                self._MgO_Va_O1_2.eigenvalues[Spin.up])

    def test_relative_total_energy(self):
        actual = self._MgO_Va_O1_2_from_vasp_files.\
            relative_total_energy(self._MgO_perfect_from_vasp_files)

        expected = -93.76904720 - -95.46878101

        self.assertEqual(actual, expected)

    def test_relative_potential(self):
        actual = self._MgO_Va_O1_2_from_vasp_files.\
            relative_potential(self._MgO_perfect_from_vasp_files,
                               self._defect_entry_MgO_Va_O1_2)

        perfect_potential = [-35.2923, -35.2923, -35.2923, -35.2923, -35.2923,
                             -35.2923, -35.2923, -35.2923, -69.8160, -69.8160,
                             -69.8160, -69.8160, -69.8160, -69.8160, -69.8160]

        expected = [x - y for x, y in
                    zip(self.electrostatic_potential, perfect_potential)]

        self.assertTrue(actual == expected)


class UnitcellDftResultsTest(unittest.TestCase):

    def setUp(self):
        """ """
        self.unitcell = UnitcellDftResults(band_edge=None,
                                           band_edge2=None,
                                           static_dielectric_tensor=None,
                                           ionic_dielectric_tensor=None,
                                           total_dos=None)

    def test_band_edge(self):
        print(self.unitcell.band_edge)
        self.unitcell.band_edge = [0.0, 1.0]
        print(self.unitcell.band_edge)
        self.unitcell.set_band_edge_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/structure_optimization"))
        print(self.unitcell.band_edge)

    def test_band_edge2(self):
        print(self.unitcell.band_edge2)
        self.unitcell.band_edge2 = [0.0, 1.0]
        print(self.unitcell.band_edge2)
        self.unitcell.set_band_edge2_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/structure_optimization"))
        print(self.unitcell.band_edge2)

    def test_dielectric_constant(self):
        self.unitcell.set_dielectric_constants_from_outcar(
            os.path.join(test_dir, "MgO/unitcell/dielectric_constants"))

        np.testing.assert_equal(self._MgO_unitcell.static_dielectric_tensor,
                                self.static_dielectric_tensor)
        np.testing.assert_equal(self._MgO_unitcell.ionic_dielectric_tensor,
                                self.ionic_dielectric_tensor)





    def test_set_static_dielectric_tensor(self):
        a = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self._MgO_unitcell.static_dielectric_tensor = a
        np.testing.assert_equal(self._MgO_unitcell.static_dielectric_tensor, a)

    def test_set_ionic_dielectric_tensor(self):
        b = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        self._MgO_unitcell.ionic_dielectric_tensor = b
        np.testing.assert_equal(self._MgO_unitcell.ionic_dielectric_tensor, b)

    def test_set_static_dielectric_tensor_from_outcar(self):
        self._MgO_unitcell.set_static_dielectric_tensor_from_outcar(
            os.path.join(test_dir, "MgO/unitcell/dielectric_constants"))

        np.testing.assert_equal(self._MgO_unitcell.static_dielectric_tensor,
                                self.static_dielectric_tensor)

    def test_set_ionic_dielectric_tensor_from_outcar(self):
        self._MgO_unitcell.set_ionic_dielectric_tensor_from_outcar(
            os.path.join(test_dir, "MgO/unitcell/dielectric_constants"))

        np.testing.assert_equal(self._MgO_unitcell.ionic_dielectric_tensor,
                                self.ionic_dielectric_tensor)

    def test_total_dielectric_tensor(self):
        self._MgO_unitcell.set_dielectric_constants_from_outcar(
            os.path.join(test_dir, "MgO/unitcell/dielectric_constants"))

        d = self.static_dielectric_tensor + self.ionic_dielectric_tensor

        np.testing.assert_equal(self._MgO_unitcell.total_dielectric_tensor, d)

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self._MgO_unitcell.to_json_file(tmp_file.name)
        MgO_unitcell_from_json = UnitcellDftResults.json_load(tmp_file.name)
        np.testing.assert_equal(MgO_unitcell_from_json.eigenvalues[Spin.up],
                                self._MgO_unitcell.eigenvalues[Spin.up])


if __name__ == "__main__":
    unittest.main()

