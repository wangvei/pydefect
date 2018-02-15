import unittest
from input_maker import *
from pydefect.input_generator.defect import Defect, IrreducibleSite
from defect_in import DefectSetting
import numpy as np
from pymatgen.core.structure import Structure

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

DEFAULT_POTCAR_DIR="/home/common/default_POTCAR"
FILENAME_PerturbAroundAPointTest_POSCAR = "examples/POSCAR-MgO64atoms"
FILENAME_POSCAR_VA1 = "examples/POSCAR-MgO64atoms-Va_Mg1"
FILENAME_POSCAR_IN1 = "examples/POSCAR-MgO64atoms-O_i1"
FILENAME_POSCAR_IN3 = "examples/POSCAR-MgO64atoms-Al_i1"
FILENAME_POSCAR_AS1 = "examples/POSCAR-MgO64atoms-O_Mg"
FILENAME_POSCAR_SS1 = "examples/POSCAR-MgO64atoms-N_O1"
#FILENAME_JSON_VA1 = "examples/MgO64atoms-Va_Mg1.json"

class NormedRandom3dVectorTest(unittest.TestCase):

    def setUp(self): 
        self.v = normed_random_3d_vector()

    def test_norm(self):
        print("normed_random_3d_vector: ", self.v)
        print("norm: ", np.linalg.norm(self.v))
        self.assertAlmostEqual(np.linalg.norm(self.v), 1.0) 


class RandomVectorTest(unittest.TestCase):

    def setUp(self): 
        self.distance = 3.0
        normed_v = normed_random_3d_vector()
        self.v = random_vector(normed_v, self.distance)

    def test_norm(self):
        print("random_3d_vector: ", self.v)
        print("distance: ", self.distance)
        print("norm: ", np.linalg.norm(self.v))
        self.assertLessEqual(np.linalg.norm(self.v), self.distance)


class PerturbAroundAPointTest(unittest.TestCase):

    def setUp(self):
        self.structure = \
            Structure.from_file(FILENAME_PerturbAroundAPointTest_POSCAR)
        self.center = [0.0, 0.0, 0.0]
        self.cutoff = 3.0
        self.distance = 0.2

    def test(self):
        # TODO: test the displacement distances
        perturbed_defect_structure, perturbed_sites = \
        perturb_around_a_point(
                        self.structure, self.center, self.cutoff, self.distance)
        true_perturbed_sites = [0, 40, 44, 48, 50, 56, 57]
        self.assertEqual(perturbed_sites, true_perturbed_sites)


class DefectMakerTest(unittest.TestCase):

    def setUp(self):
        self.structure = \
            Structure.from_file(FILENAME_PerturbAroundAPointTest_POSCAR)
        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              repr_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             repr_coords=[0.25, 0.25, 0.25])
        self.irreducible_sites = [Mg1, O1]
        self.interstitial_coords = [[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]

        _test_structure = Structure.from_file(FILENAME_POSCAR_VA1)
        _removed_atom_index = 0
        _inserted_atom_index = None
        _defect_coords = [0, 0, 0]
        _in_name = "Va"
        _out_name = "Mg1"
        _charge = -2

        self.test_d = Defect(_test_structure, _removed_atom_index,
                             _inserted_atom_index, _defect_coords,
                             _in_name, _out_name, _charge)

    def test_vacancies(self):
        vac1_d = \
            DefectMaker("Va_Mg1_-2", self.structure, self.irreducible_sites,
                        self.interstitial_coords)
        self.assertEqual(vac1_d.defect.as_dict(), self.test_d.as_dict())


#class DefectInputMakerTest(unittest.TestCase):
#
#    def setUp(self):
#        structure = Structure.from_file(FILENAME_PerturbAroundAPointTest_POSCAR)
#        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
#                              first_index=1, last_index=32,
#                              repr_coords=[0, 0, 0])
#        O1 = IrreducibleSite(irreducible_name="O1", element="O",
#                             first_index=33, last_index=64,
#                             repr_coords=[0.25, 0.25, 0.25])
#        irrep_elements = [Mg1, O1]
#        dopant_configs = [["Al", "Mg"], ["Al", "O"], ["N","Mg"], ["N", "O"]]
#        antisite_configs = [["Mg", "O"], ["O","Mg"]]
#        interstitial_coords = [[0.1, 0.1, 0.1]]
#        self.interstitial_coords = interstitial_coords
#        included = ["Va_O1_-1", "Va_O1_-2"]
#        excluded = ["Va_O1_1", "Va_O1_2"]
#        distance = 0.15
#        cutoff = 2.0
#        symprec = 0.001
#        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
#        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}
#
#        self._mgo = DefectSetting(
#                   structure, irrep_elements, dopant_configs, antisite_configs,
#                   interstitial_coords, included, excluded, distance, cutoff,
#                   symprec, oxidation_states, electronegativity)
#
#    def test_vacancies(self):
#        vacancy1 = "Va_Mg1_-2"
#        vacancy2 = "Va_Mg2_-2"
#        va1 = DefectInputMaker(vacancy1, self._mgo)
#        va1_check_structure = Structure.from_file(FILENAME_POSCAR_VA1)
#        self.assertTrue(va1.defect_structure, va1_check_structure)
#        with self.assertRaises(ValueError):
#            va2 = DefectInputMaker(vacancy2, self._mgo)
#
#        self.assertEqual(va1.defect_index, 0)
#        self.assertEqual(va1.defect_coords, [0.0, 0.0, 0.0])
#
#    def test_interstitials(self):
#        interstital1 = "O_i1_2"
#        interstital2 = "O_i2_-2"
#        interstital3 = "Al_i1_3"
#        in1 = DefectInputMaker(interstital1, self._mgo)
#        in1_check_structure = Structure.from_file(FILENAME_POSCAR_IN1)
#        self.assertTrue(in1.defect_structure, in1_check_structure)
#        with self.assertRaises(ValueError):
#            DefectInputMaker(interstital2, self._mgo)
#        in3 = DefectInputMaker(interstital3, self._mgo)
#        in3_check_structure = Structure.from_file(FILENAME_POSCAR_IN3)
#        self.assertTrue(in3.defect_structure, in3_check_structure)
#
#        self.assertEqual(in1.defect_index, 32)
#        self.assertEqual(in3.defect_index, 0)
#
#        self.assertEqual(in1.defect_coords, self.interstitial_coords[0])
#        self.assertEqual(in3.defect_coords, self.interstitial_coords[0])
#
#    def test_antisites(self):
#        antisite1 = "O_Mg1_0"
#        as1 = DefectInputMaker(antisite1, self._mgo)
#        as1_check_structure = Structure.from_file(FILENAME_POSCAR_AS1)
#        self.assertTrue(as1.defect_structure, as1_check_structure)
#
#    def test_substituted(self):
#        substituted1 = "N_O1_-2"
#        substituted2 = "N_O2_-2"
#        ss1 = DefectInputMaker(substituted1, self._mgo)
#        ss1_check_structure = Structure.from_file(FILENAME_POSCAR_SS1)
#        self.assertTrue(ss1.defect_structure, ss1_check_structure)
#        with self.assertRaises(ValueError):
#            ss2 = DefectInputMaker(substituted2, self._mgo)
#
#
#class DefectInputSetMakerTest(unittest.TestCase):
#
#    def setUp(self):
#        structure = Structure.from_file(FILENAME_PerturbAroundAPointTest_POSCAR)
#        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
#                              first_index=1, last_index=32,
#                              repr_coords=[0, 0, 0])
#        O1 = IrreducibleSite(irreducible_name="O1", element="O",
#                             first_index=33, last_index=64,
#                             repr_coords=[0.25, 0.25, 0.25])
#        irrep_elements = [Mg1, O1]
#        dopant_configs = [["Al", "Mg"]]
#        antisite_configs = [["O","Mg"]]
#        interstitial_coords = [[0.1, 0.1, 0.1]]
#        self.interstitial_coords = interstitial_coords
#        included = ["Va_O1_-1", "Va_O1_-2"]
#        excluded = ["Va_O1_1", "Va_O1_2"]
#        distance = 0.15
#        cutoff = 2.0
#        symprec = 0.001
#        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
#        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}
#
#        self._mgo = DefectSetting(
#            structure, irrep_elements, dopant_configs, antisite_configs,
#            interstitial_coords, included, excluded, distance, cutoff,
#            symprec, oxidation_states, electronegativity)
#
#    def test(self):
#        mgo = DefectInputSetMaker(self._mgo)
#        check_defect_set = \
#            ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0',
#             'Va_O1_0', 'Va_O1_-1', 'Va_O1_-2',
#             'Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2',
#             'O_i1_-2', 'O_i1_-1', 'O_i1_0',
#             'O_Mg1_-4', 'O_Mg1_-3', 'O_Mg1_-2', 'O_Mg1_-1', 'O_Mg1_0',
#             'Al_Mg1_0', 'Al_Mg1_1']
#        self.assertTrue(mgo.defect_set, check_defect_set)

