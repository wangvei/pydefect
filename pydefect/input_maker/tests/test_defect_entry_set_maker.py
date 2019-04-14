# -*- coding: utf-8 -*-

from io import StringIO
import os
import sys
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.defect_entry_set_maker \
    import get_int_from_string, parse_defect_name, log_is_being_removed, \
    log_already_exist, log_is_being_constructed, is_name_matched, \
    select_defect_names

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.irreducible_site import IrreducibleSite

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


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


class FilterNameTest(unittest.TestCase):
    def test(self):
        self.assertTrue(is_name_matched("Va_O11_-2",
                                        keywords=["Va_O[0-9]+_-[0-9]+"]))
        self.assertFalse(is_name_matched("Mg_i1_0", keywords=["Va", "O_i1"]))
        self.assertTrue(is_name_matched("Mg_i1_0", keywords=("Va", "Mg_i1")))
        self.assertFalse(is_name_matched("Mg_i1_0", keywords="Va"))


class FilterNameSetTest(unittest.TestCase):

    def test(self):
        name_set = \
            ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_0', 'Va_O1_-1',
             'Va_O1_-2', 'Va_O2_0', 'Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2',
             'O_i1_-1', 'O_i1_0', 'O_Mg1_-4', 'O_Mg1_-3', 'O_Mg1_-2',
             'O_Mg1_-1', 'O_Mg1_0', 'Al_Mg1_0', 'Al_Mg1_1']

        actual_Va = select_defect_names(name_set, ["Va"])
        expected_Va = ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_0',
                       'Va_O1_-1', 'Va_O1_-2', 'Va_O2_0']

        actual__i = select_defect_names(name_set, ["_i"])
        expected__i = ['Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2', 'O_i1_-1',
                       'O_i1_0']

        actual_Va_and__i = select_defect_names(name_set, ["Va", "_i"])
        expected_Va_and__i = expected_Va + expected__i

        actual_Va_O = select_defect_names(name_set, ["Va_O"])
        expected_Va_O = ['Va_O1_0', 'Va_O1_-1', 'Va_O1_-2', 'Va_O2_0']

        actual_Va_O_0 = select_defect_names(name_set, ["Va_O[0-9]_0"])
        expected_Va_O_0 = ['Va_O1_0', 'Va_O2_0']

        actual_Va_O1 = select_defect_names(name_set, ["Va_O1"])
        expected_Va_O1 = ['Va_O1_0', 'Va_O1_-1', 'Va_O1_-2']

        self.assertEqual(sorted(actual_Va), sorted(expected_Va))
        self.assertEqual(sorted(actual__i), sorted(expected__i))
        self.assertEqual(sorted(actual_Va_and__i), sorted(expected_Va_and__i))
        self.assertEqual(sorted(actual_Va_O), sorted(expected_Va_O))
        self.assertEqual(sorted(actual_Va_O_0), sorted(expected_Va_O_0))
        self.assertEqual(sorted(actual_Va_O1), sorted(expected_Va_O1))

import os
import shutil
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.defect_entry_set_maker import DefectEntrySetMaker
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting

class DefectEntrySetMakerTest(unittest.TestCase):

    def setUp(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))
        space_group_symbol = "Fm-3m"
        Mg1 = IrreducibleSite(irreducible_name="Mg1",
                              element="Mg",
                              first_index=1,
                              last_index=32,
                              representative_coords=[0, 0, 0],
                              wyckoff="a",
                              site_symmetry="m3m",
                              coordination_distances=None)
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             representative_coords=[0.25, 0.25, 0.25],
                             wyckoff="a",
                             site_symmetry="m3m",
                             coordination_distances=None)
        irreducible_elements = [Mg1, O1]
        dopant_configs = [["Al", "Mg"]]
        antisite_configs = [["Mg", "O"], ["O", "Mg"]]
        interstitial_coords = [[0.1, 0.1, 0.1]]
        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        distance = 0.15
        cutoff = 2.3
        symprec = 0.001
        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}

        self._mgo = DefectInitialSetting(
            structure, space_group_symbol, irreducible_elements, dopant_configs,
            antisite_configs, interstitial_coords, included, excluded, distance,
            cutoff, symprec, oxidation_states, electronegativity)

    def test_mgo(self):
        desm = DefectEntrySetMaker(defect_initial_setting=self._mgo,
                                   particular_defects=["Sc_Mg1_1", "Sc_Mg1_0"])
        des = desm.defect_entries

        for de in des:
            print(de)

# class DefectEntryMakerTest(unittest.TestCase):

    # def setUp(self):
    #     self.structure = \
    #         Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))

        # Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
        #                       first_index=1, last_index=32,
        #                       representative_coords=[0, 0, 0])
        # O1 = IrreducibleSite(irreducible_name="O1", element="O",
        #                      first_index=33, last_index=64,
        #                      representative_coords=[0.25, 0.25, 0.25])
        # self.irreducible_sites = [Mg1, O1]
        # self.interstitial_coords = [[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]

        # name = "Va_Mg1"
        # test_structure = Structure.from_file(
        #    os.path.join(test_dir, "POSCAR-MgO64atoms-Va_Mg1"))
        # removed_atoms = {0: [0, 0, 0]}
        # inserted_atoms = []
        # element_diff = {"Mg": -1}
        # charge = -2

        # self.test_d = DefectEntry(name, test_structure, removed_atoms,
        #                           inserted_atoms, element_diff, charge)

    # def test_vacancies(self):
    #     vac1_d = \
    #         DefectEntryMaker("Va_Mg1_-2", self.structure, self.irreducible_sites,
    #                          self.interstitial_coords)
    #     self.assertEqual(vac1_d.defect.as_dict(), self.test_d.as_dict())



class DefectStructureSetMakerTest(unittest.TestCase):
    def setUp(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))

        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              representative_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             representative_coords=[0.25, 0.25, 0.25])
        irreducible_elements = [Mg1, O1]
        dopant_configs = [["Al", "Mg"]]
        antisite_configs = [["Mg", "O"], ["O", "Mg"]]
        interstitial_coords = [[0.1, 0.1, 0.1]]
        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        distance = 0.15
        cutoff = 2.0
        symprec = 0.001
        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}

        self._mgo = DefectInitialSetting(
            structure, irreducible_elements, dopant_configs, antisite_configs,
            interstitial_coords, included, excluded, distance, cutoff,
            symprec, oxidation_states, electronegativity)

    def test_mgo(self):
        dssm = DefectStructureSetMaker(defect_initial_setting=self._mgo,
                                       particular_defects=["Sc_Mg1_0"])
        print(dssm)

    # def test_vo_mgo(self):
    #     if os.path.exists(test_vo_mgo_dir):
    #         shutil.rmtree(test_vo_mgo_dir)
    #     # Note that the type of filtering_words is list.
    #     VaspDefectInputSetMaker(defect_initial_setting=self.MgO,
    #                             keywords=["perfect", "Va_O"])

if __name__ == "__main__":
    unittest.main()
