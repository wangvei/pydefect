# -*- coding: utf-8 -*-

import numpy as np
import os
import unittest

from pydefect.util.structure import structure_to_spglib_cell, \
    spglib_cell_to_structure, find_equivalent_sites,\
    find_hpkot_primitive, structure_to_seekpath, perturb_neighbors, \
    NotStandardizedPrimitiveError
from pymatgen.core.structure import Structure

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "April 4, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class StructureToSpglibCellTest(unittest.TestCase):

    def setUp(self):
        structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(structure)

    def test(self):
        expected_lattice = \
            [[4.2468930196252002, 0.0000000000000000, 0.0000000000000000],
             [0.0000000000000000, 4.2468930196252002, 0.0000000000000000],
             [0.0000000000000000, 0.0000000000000000, 4.2468930196252002]]
        expected_1st_position = [0, 0, 0]
        expected_atomic_numbers = [12, 12, 12, 12, 8, 8, 8, 8]

        self.assertTrue(expected_lattice, self.cell[0])
        self.assertTrue(expected_1st_position, self.cell[1][0])
        self.assertTrue(expected_atomic_numbers, self.cell[2])


class SpglibCellToStructureTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(self.structure)
        self.returned_structure = spglib_cell_to_structure(self.cell)

    def test(self):
        self.assertTrue(self.structure == self.returned_structure)


class SpglibCellToStructureTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("SPOSCAR")
#        self.structure = Structure.from_file("POSCAR-SiPO")

    def test(self):
        s, r, e = find_equivalent_sites(self.structure)
        # print(s)
        # print(r)
        # print(e)



class FindPrimitiveTest(unittest.TestCase):
    def test(self):
        expected = Structure.from_file("PPOSCAR-MgO")
        actual = find_hpkot_primitive(Structure.from_file("BPOSCAR-MgO"))

        from pymatgen.analysis.structure_matcher import StructureMatcher
        sm = StructureMatcher(ltol=0.0001, stol=0.0001, angle_tol=0.001,
                              primitive_cell=False, scale=False)
        self.assertTrue(sm.fit(expected, actual))


class StructureToSeekpathTest(unittest.TestCase):

    def setUp(self):
        self.conventional_cell_structure = \
            Structure.from_file("PPOSCAR-MgO-strange_lattice")
        self.primitive_cell_structure = Structure.from_file("PPOSCAR-MgO")

    def test_raise_error(self):
        with self.assertRaises(NotStandardizedPrimitiveError) as cm:
            structure_to_seekpath(self.conventional_cell_structure)

    def test_return_res(self):
        res = structure_to_seekpath(self.primitive_cell_structure)
        expected = np.array([0.5, 0, 0.5])
        np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
        self.assertEqual("X", res["explicit_kpoints_labels"][-1])


class PerturbNeighborsTest(unittest.TestCase):

    def test(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))
        center = [0.0, 0.0, 0.0]
        cutoff = 3.0
        distance = 0.2

        # TODO: test the displacement distances
        perturbed_defect_structure, perturbed_sites = \
            perturb_neighbors(structure, center, cutoff, distance)
        true_perturbed_sites = [0, 40, 44, 48, 50, 56, 57]
        self.assertEqual(perturbed_sites, true_perturbed_sites)


class StructureToSeekpathTest(unittest.TestCase):

    def setUp(self):
        self.conventional_cell_structure = \
            Structure.from_file("PPOSCAR-MgO-strange_lattice")
        self.primitive_cell_structure = Structure.from_file("PPOSCAR-MgO")

    def test_raise_error(self):
        with self.assertRaises(NotStandardizedPrimitiveError) as cm:
            structure_to_seekpath(self.conventional_cell_structure)

    def test_return_res(self):
        res = structure_to_seekpath(self.primitive_cell_structure)
        expected = np.array([0.5, 0, 0.5])
        np.testing.assert_equal(expected, res["explicit_kpoints_rel"][-1])
        self.assertEqual("X", res["explicit_kpoints_labels"][-1])


class PerturbNeighborsTest(unittest.TestCase):

    def test(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))
        center = [0.0, 0.0, 0.0]
        cutoff = 3.0
        distance = 0.2

        # TODO: test the displacement distances
        perturbed_defect_structure, perturbed_sites = \
            perturb_neighbors(structure, center, cutoff, distance)
        true_perturbed_sites = [0, 40, 44, 48, 50, 56, 57]
        self.assertEqual(perturbed_sites, true_perturbed_sites)


if __name__ == "__main__":
    unittest.main()