# -*- coding: utf-8 -*-

from copy import deepcopy
import numpy as np
import os
import tempfile
import unittest

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.util.testing import PymatgenTest

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results \
    import SupercellCalcResults
from pydefect.util.structure_tools import defect_center, distances_from_point
from pydefect.core.unitcell_calc_results import UnitcellCalcResults

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

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

    def test_defect_center(self):
        # Test for periodic boundary condition (PBC).
        name = "MgO64atoms-Va_Mg1+Va_O1"
        test_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO64atoms-Va_Mg1+Va_O1"))
        removed_atoms = {0: [0, 0, 0], 39: [0.75, 0.75, 0.75]}
        inserted_atoms = []
        element_diff = {"Mg": -1, "O": -1}
        charge = 0
        # Construct DefectEntry class object.
        pbc_defect_entry = DefectEntry(name, test_structure, removed_atoms,
                                        inserted_atoms, element_diff, charge)

        # Test a single defect.
        actual4 = defect_center(pbc_defect_entry)
        expected = [-0.125, -0.125, -0.125]
        self.assertEqual(actual4, expected)


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

        self._MgO_Va_O1_2 = \
            SupercellCalcResults.from_vasp_files(
                os.path.join(test_dir, "MgO/defects/Va_O1_2"))

        self._MgO_perfect = \
            SupercellCalcResults.from_vasp_files(
                os.path.join(test_dir, "MgO/defects/perfect"))

        self.d_from_vasp_files = self._MgO_Va_O1_2.as_dict

        name = "Va_O1"
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "MgO/defects/Va_O1_2", "POSCAR"))
        perturbed_initial_structure = initial_structure.copy()
        removed_atoms = {8: [0.25, 0.25, 0.25]}
        inserted_atoms = {}
        changes_of_num_elements = {"O": -1}
        charge = 2
        initial_symmetry = "Oh"
        multiplicity = 4
        perturbed_sites = []

        self._defect_entry_MgO_Va_O1_2 = \
            DefectEntry(name, initial_structure, perturbed_initial_structure, removed_atoms, inserted_atoms,
                        changes_of_num_elements, charge,initial_symmetry,
                        multiplicity, perturbed_sites)

    def test_print(self):
        print(self._MgO_Va_O1_2)

    def test_from_vasp_files(self):
        # CAUTION: When constructing Structure object from Structure.from_file
        #          velocities are not stored.
        #          Therefore, equality check of Structure objects returns False.
        #          If the structure is converted via poscar file format, it may
        #          be solved.

        # energy
        expected = -93.64519527
        self.assertEqual(self.d_from_vasp_files["total_energy"], expected)

        # total_magnetization
        expected = 3.4e-06
        self.assertEqual(self.d_from_vasp_files["total_magnetization"], expected)

        # eigenvalue: test only a single point
        expected = [-1.43572e+01, 1.0]
        self.assertEqual(self.d_from_vasp_files["eigenvalues"]["1"][0][0],
                         expected)

        # electrostatic_potential
        expected = [-34.5752, -35.5301, -35.5344, -35.5357, -35.5346, -35.5326,
                    -35.5265, -34.5754, -70.027, -70.026, -70.0253, -70.0274,
                    -70.0271, -70.0272, -70.4853]
        self.assertEqual(self.d_from_vasp_files["electrostatic_potential"],
                        expected)

    def test_dict(self):
        MgO_Va_O1_2_fd = SupercellCalcResults.from_dict(self.d_from_vasp_files)
        print(MgO_Va_O1_2_fd)
        np.testing.assert_equal(MgO_Va_O1_2_fd.eigenvalues[Spin.up],
                                self._MgO_Va_O1_2.eigenvalues[Spin.up])

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self._MgO_Va_O1_2.to_json_file(tmp_file.name)
        MgO_Va_O1_2_from_json = SupercellCalcResults.load_json(tmp_file.name)
        np.testing.assert_equal(MgO_Va_O1_2_from_json.eigenvalues[Spin.up],
                                self._MgO_Va_O1_2.eigenvalues[Spin.up])

    def test_relative_total_energy(self):
        actual = self._MgO_Va_O1_2.\
            relative_total_energy(self._MgO_perfect)

        expected = -93.64519527 - -95.36395670

        self.assertEqual(actual, expected)

    def test_relative_potential(self):
        actual = self._MgO_Va_O1_2.\
            relative_potential(self._MgO_perfect,
                               self._defect_entry_MgO_Va_O1_2)

        perfect_potential = \
            [-35.2983, -35.2983, -35.2983, -35.2983, -35.2983, -35.2983,
             -35.2983, -35.2983, -69.7919, -69.7919, -69.7919, -69.7919,
             -69.7919, -69.7919, -69.7919]

        expected = [x - y for x, y in
                    zip(self._MgO_Va_O1_2.electrostatic_potential,
                        perfect_potential)]

        self.assertTrue(actual == expected)


if __name__ == "__main__":
    unittest.main()

