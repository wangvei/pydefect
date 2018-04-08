# -*- coding: utf-8 -*-

from io import StringIO
import numpy as np
import os
import sys
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.input_maker import normalized_random_3d_vector, \
    random_vector, perturb_around_a_point, get_int_from_string, \
    parse_defect_name, print_already_exist, print_is_being_constructed, \
    filter_name_set, DefectMaker

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.irreducible_site import IrreducibleSite

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

DEFAULT_POTCAR_DIR = "/home/common/default_POTCAR"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class NormalizedRandom3dVectorTest(unittest.TestCase):

    def setUp(self): 
        self.v = normalized_random_3d_vector()

    def test_norm(self):
        print("normalized_random_3d_vector: ", self.v)
        print("norm: ", np.linalg.norm(self.v))
        self.assertAlmostEqual(np.linalg.norm(self.v), 1.0) 


class RandomVectorTest(unittest.TestCase):

    def setUp(self): 
        self.distance = 3.0
        normalized_v = normalized_random_3d_vector()
        self.v = random_vector(normalized_v, self.distance)

    def test_norm(self):
        print("random_3d_vector: ", self.v)
        print("distance: ", self.distance)
        print("norm: ", np.linalg.norm(self.v))
        self.assertLessEqual(np.linalg.norm(self.v), self.distance)


class PerturbAroundAPointTest(unittest.TestCase):

    def test(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))
        center = [0.0, 0.0, 0.0]
        cutoff = 3.0
        distance = 0.2
        # TODO: test the displacement distances
        perturbed_defect_structure, perturbed_sites = \
            perturb_around_a_point(structure, center, cutoff, distance)
        true_perturbed_sites = [0, 40, 44, 48, 50, 56, 57]
        self.assertEqual(perturbed_sites, true_perturbed_sites)


class GetIntFromStringTest(unittest.TestCase):

    def test(self):
        actual_1 = get_int_from_string("Mg1")
        actual_2 = get_int_from_string("Mg1.Al2")
        expected_1 = 1
        expected_2 = 12
        self.assertEqual(expected_1, actual_1)
        self.assertEqual(expected_2, actual_2)


class ParseDefectNameTest(unittest.TestCase):

    def test(self):
        actual_1 = parse_defect_name("Va_Mg1_0")
        expected_1 = ("Va", "Mg1", 0)
        self.assertEqual(actual_1, expected_1)


class PrintAlreadyExistTest(unittest.TestCase):

    def setUp(self):
        self.org_stdout, sys.stdout = sys.stdout, StringIO()

    def tearDown(self):
        sys.stdout = self.org_stdout

    def test(self):
        print_already_exist(name="Va_O1_0")
        expected = "   Va_O1_0 already exists, so nothing is done.\n"
        self.assertEqual(sys.stdout.getvalue(), expected)


class PrintIsBeingConstructedTest(unittest.TestCase):

    def setUp(self):
        self.org_stdout, sys.stdout = sys.stdout, StringIO()

    def tearDown(self):
        sys.stdout = self.org_stdout

    def test(self):
        print_is_being_constructed(name="Va_O1_0")
        expected = "   Va_O1_0 is being constructed.\n"
        self.assertEqual(sys.stdout.getvalue(), expected)


class FilterDefectNameSetTest(unittest.TestCase):

    def test(self):
        name_set = \
            ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_0', 'Va_O1_-1',
             'Va_O1_-2', 'Va_O2_0', 'Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2',
             'O_i1_-1', 'O_i1_0', 'O_Mg1_-4', 'O_Mg1_-3', 'O_Mg1_-2',
             'O_Mg1_-1', 'O_Mg1_0', 'Al_Mg1_0', 'Al_Mg1_1']
        expected_Va = ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_0',
                       'Va_O1_-1', 'Va_O1_-2', 'Va_O2_0']
        expected__i = ['Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2', 'O_i1_-1',
                       'O_i1_0']
        expected_Va_and__i = expected_Va + expected__i
        expected_Va_O = ['Va_O1_0', 'Va_O1_-1', 'Va_O1_-2', 'Va_O2_0']
        expected_Va_O1 = ['Va_O1_0', 'Va_O1_-1', 'Va_O1_-2']

        actual_Va = filter_name_set(name_set, ["Va"])
        actual__i = filter_name_set(name_set, ["_i"])
        actual_Va_and__i = filter_name_set(name_set, ["Va", "_i"])
        actual_Va_O = filter_name_set(name_set, ["Va_O"])
        actual_Va_O1 = filter_name_set(name_set, ["Va_O1"])

        self.assertEqual(sorted(actual_Va), sorted(expected_Va))
        self.assertEqual(sorted(actual__i), sorted(expected__i))
        self.assertEqual(sorted(actual_Va_and__i), sorted(expected_Va_and__i))
        self.assertEqual(sorted(actual_Va_O), sorted(expected_Va_O))
        self.assertEqual(sorted(actual_Va_O1), sorted(expected_Va_O1))


class DefectMakerTest(unittest.TestCase):

    def setUp(self):
        self.structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))

        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              representative_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             representative_coords=[0.25, 0.25, 0.25])
        self.irreducible_sites = [Mg1, O1]
        self.interstitial_coords = [[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]

        _test_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO64atoms-Va_Mg1"))
        _removed_atoms = {0: [0, 0, 0]}
        _inserted_atoms = []
        _element_diff = {"Mg": -1}
        _charge = -2

        self.test_d = DefectEntry(_test_structure, _removed_atoms,
                                  _inserted_atoms, _element_diff, _charge)

    def test_vacancies(self):
        vac1_d = \
            DefectMaker("Va_Mg1_-2", self.structure, self.irreducible_sites,
                        self.interstitial_coords)
        self.assertEqual(vac1_d.defect.as_dict(), self.test_d.as_dict())


if __name__ == "__main__":
    unittest.main()
