#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import unittest

from pydefect.util.structure import structure_to_spglib_cell, \
    spglib_cell_to_structure, find_primitive, structure2seekpath, \
    perturb_neighbors
from pymatgen import Structure

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

    def test(self):
        structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(structure)


class SpglibCellToStructureTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("BPOSCAR-MgO")
        self.cell = structure_to_spglib_cell(self.structure)

    def test(self):
        s = spglib_cell_to_structure(self.cell)
        self.assertTrue(s == self.structure)


class FindPrimitiveTest(unittest.TestCase):
    def test(self):
        actual = find_primitive(Structure.from_file("BPOSCAR-MgO"))
        print(actual)
        expected = Structure.from_file("PPOSCAR-MgO")
        print(expected)
        self.assertTrue(actual == expected)


class Structure2SeekpathTest(unittest.TestCase):

    def setUp(self):
#        self.structure = Structure.from_file("PPOSCAR-YMnO3")
        self.structure = Structure.from_file("PPOSCAR-YMnO3-bc_exchanged")

    def test_structure2seekpath(self):
#        structure2seekpath(self.structure)
        res = structure2seekpath(self.structure)
        print(res["primitive_lattice"])


class PerturbAroundAPointTest(unittest.TestCase):

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