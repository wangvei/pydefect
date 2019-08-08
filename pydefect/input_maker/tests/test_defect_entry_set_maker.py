# -*- coding: utf-8 -*-

import os
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.defect_entry_set_maker \
    import get_int_from_string, select_defect_names
from pydefect.input_maker.defect_entry_set_maker import DefectEntrySetMaker
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.core.defect_name import SimpleDefectName

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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


class SelectDefectNamesTest(unittest.TestCase):

    def test(self):
        name_set = \
            ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_0', 'Va_O1_-1',
             'Va_O1_-2', 'Va_O2_0', 'Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2',
             'O_i1_-1', 'O_i1_0', 'O_Mg1_-4', 'O_Mg1_-3', 'O_Mg1_-2',
             'O_Mg1_-1', 'O_Mg1_0', 'Al_Mg1_0', 'Al_Mg1_1']

        name_set = [SimpleDefectName.from_str(name) for name in name_set]

        actual_va = select_defect_names(name_set, ["Va"], return_str=True)
        expected_va = ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_0',
                       'Va_O1_-1', 'Va_O1_-2', 'Va_O2_0']

        actual__i = select_defect_names(name_set, ["_i"], return_str=True)
        expected__i = ['Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2', 'O_i1_-1',
                       'O_i1_0']

        actual_va_and__i = \
            select_defect_names(name_set, ["Va", "_i"], return_str=True)
        expected_va_and__i = expected_va + expected__i

        actual_va_o = select_defect_names(name_set, ["Va_O"], return_str=True)
        expected_va_o = ['Va_O1_0', 'Va_O1_-1', 'Va_O1_-2', 'Va_O2_0']

        actual_va_o_0 = \
            select_defect_names(name_set, ["Va_O[0-9]_0"], return_str=True)
        expected_va_o_0 = ['Va_O1_0', 'Va_O2_0']

        actual_va_o1 = select_defect_names(name_set, ["Va_O1"], return_str=True)
        expected_va_o1 = ['Va_O1_0', 'Va_O1_-1', 'Va_O1_-2']

        self.assertEqual(sorted(actual_va), sorted(expected_va))
        self.assertEqual(sorted(actual__i), sorted(expected__i))
        self.assertEqual(sorted(actual_va_and__i), sorted(expected_va_and__i))
        self.assertEqual(sorted(actual_va_o), sorted(expected_va_o))
        self.assertEqual(sorted(actual_va_o_0), sorted(expected_va_o_0))
        self.assertEqual(sorted(actual_va_o1), sorted(expected_va_o1))


class DefectEntrySetMakerTest(unittest.TestCase):

    def setUp(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))
        space_group_symbol = "Fm-3m"
        mg1 = IrreducibleSite(irreducible_name="Mg1",
                              element="Mg",
                              first_index=1,
                              last_index=32,
                              representative_coords=[0, 0, 0],
                              wyckoff="a",
                              site_symmetry="m3m",
                              coordination_distances=None)
        o1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             representative_coords=[0.25, 0.25, 0.25],
                             wyckoff="a",
                             site_symmetry="m3m",
                             coordination_distances=None)
        irreducible_elements = [mg1, o1]
        dopant_configs = [["Al", "Mg"]]
        antisite_configs = [["Mg", "O"], ["O", "Mg"]]

        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        distance = 0.15
        cutoff = 2.3
        symprec = 0.001
        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}

        self._mgo = DefectInitialSetting(
            structure=structure,
            space_group_symbol=space_group_symbol,
            transformation_matrix= [2, 2, 2],
            cell_multiplicity=32,
            irreducible_sites=irreducible_elements,
            dopant_configs=dopant_configs,
            antisite_configs=antisite_configs,
            interstitial_sites=["i1"],
            included=included,
            excluded=excluded,
            displacement_distance=distance,
            cutoff=cutoff,
            symprec=symprec,
            angle_tolerance=5,
            oxidation_states=oxidation_states,
            electronegativity=electronegativity)

    def test_mgo(self):
        desm = DefectEntrySetMaker(defect_initial_setting=self._mgo,
                                   particular_defects=["Sc_Mg1_1", "Sc_Mg1_0"])
        des = desm.defect_entries

        for de in des:
            print(de)


if __name__ == "__main__":
    unittest.main()
