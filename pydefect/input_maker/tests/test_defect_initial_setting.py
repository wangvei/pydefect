# -*- coding: utf-8 -*-
import os
import tempfile

from pymatgen.core.structure import Structure

from pydefect.input_maker.defect_initial_setting import (
    candidate_charge_set, get_electronegativity, get_oxidation_state,
    dopant_info, get_distances_from_string, insert_atoms, select_defects,
    DefectInitialSetting)
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.util.testing import PydefectTest
from pydefect.core.defect_entry import DefectType

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class CandidateChargeSetTest(PydefectTest):
    def test_range1(self):
        actual_positive = [i for i in candidate_charge_set(2)]
        expected_positive = [0, 1, 2]
        self.assertEqual(actual_positive, expected_positive)

    def test_range2(self):
        actual_negative = [i for i in candidate_charge_set(-2)]
        expected_negative = [-2, -1, 0]
        self.assertEqual(actual_negative, expected_negative)

    def test_range3(self):
        actual_positive = [i for i in candidate_charge_set(3)]
        expected_positive = [-1, 0, 1, 2, 3]
        self.assertEqual(actual_positive, expected_positive)

    def test_range4(self):
        actual_negative = [i for i in candidate_charge_set(-3)]
        expected_negative = [-3, -2, -1, 0, 1]
        self.assertEqual(actual_negative, expected_negative)


class GetElectronegativityTest(PydefectTest):
    def test_success(self):
        true_element = "Mg"
        self.assertEqual(get_electronegativity(true_element), 1.31)

    def test_fail(self):
        fake_element = "Yk"
        self.assertEqual(get_electronegativity(fake_element), 0)


class GetOxidationStateTest(PydefectTest):
    def test_success(self):
        true_element = "Mg"
        self.assertEqual(get_oxidation_state(true_element), 2)

    def test_fail(self):
        fake_element = "Yk"
        self.assertEqual(get_oxidation_state(fake_element), 0)


class DopantInfoTest(PydefectTest):
    def test_success(self):
        true_element = "Mg"
        expected = """   Dopant element: Mg
Electronegativity: 1.31
  Oxidation state: 2"""
        self.assertEqual(dopant_info(true_element), expected)

    def test_fail(self):
        fake_element = "Ak"
        self.assertEqual(dopant_info(fake_element), None)


class GetDistancesFromStringTest(PydefectTest):
    def test_success(self):
        fine_string_list = "Mg: 2.1 2.2 O: 2.3 2.4".split()
        expected = {"Mg": [2.1, 2.2], "O": [2.3, 2.4]}
        self.assertEqual(get_distances_from_string(fine_string_list), expected)

    def test_fail1(self):
        bad_string_list = "Mg 2.1 2.2 O 2.3 2.4".split()
        with self.assertRaises(KeyError):
            get_distances_from_string(bad_string_list)

    def test_fail2(self):
        bad_string_list = "Mg: 2.1 a O: 2.3 2.4".split()
        with self.assertRaises(ValueError):
            get_distances_from_string(bad_string_list)


class InsertAtoms(PydefectTest):
    def setUp(self) -> None:
        structure = self.get_structure_by_name("MgO64atoms")
        atoms = [{"element": "Al", "coords": [0, 0, 0.01]},
                 {"element": "Na", "coords": [0, 0, 0.02]},
                 {"element": "Mg", "coords": [0, 0, 0.03]},
                 {"element": "O",  "coords": [0, 0, 0.04]}]
        self.inserted_structure, self.inserted_atoms \
            = insert_atoms(structure, atoms)

    def test(self):
        expected = [{'element': 'Al', 'index': 1, 'coords': [0, 0, 0.01]},
                    {'element': 'Na', 'index': 0, 'coords': [0, 0, 0.02]},
                    {'element': 'Mg', 'index': 2, 'coords': [0, 0, 0.03]},
                    {'element': 'O', 'index': 35, 'coords': [0, 0, 0.04]}]
        self.assertEqual(expected, self.inserted_atoms)
        self.assertEqual("Al", str(self.inserted_structure[1].specie))
        self.assertArrayAlmostEqual(
            [0, 0, 0.04], list(self.inserted_structure[35].frac_coords), 5)


class SelectDefectsTest(PydefectTest):

    def setUp(self) -> None:
        self.name_set = [{"name": 'Va_Mg1', "charge": -1},
                         {"name": 'Va_Mg1', "charge": 0},
                         {"name": 'Va_O1', "charge": 0},
                         {"name": 'O_i1', "charge": 0},
                         {"name": 'O_i1', "charge": 1}]

    def test_keywords1(self):
        actual_va = select_defects(self.name_set, keywords=["Va"])
        expected_va = [{'name': 'Va_Mg1', 'charge': -1},
                       {'name': 'Va_Mg1', 'charge': 0},
                       {'name': 'Va_O1', 'charge': 0}]
        self.assertEqual(expected_va, actual_va)

    def test_keywords2(self):
        actual_va = select_defects(self.name_set, keywords=["Mg_i1"])
        self.assertEqual([], actual_va)

    def test_specified(self):
        actual_specified = \
            select_defects(self.name_set, specified_defects=["Va_Mg1_-1"])
        expected_specified = [{'name': 'Va_Mg1', 'charge': -1}]
        self.assertEqual(expected_specified, actual_specified)

    def test_specified_not_exit(self):
        with self.assertRaises(ValueError):
            select_defects(self.name_set, specified_defects=["Va_Mg2_-1"])

    def test_set_both_fail(self):
        with self.assertRaises(ValueError):
            select_defects(self.name_set, keywords=["Va"],
                           specified_defects=["Va_Mg1_-1"])


class DefectInitialSettingTest(PydefectTest):

    def setUp(self):
        """
        """
        self.structure = self.get_structure_by_name("MgO64atoms")
        space_group_symbol = "Fm-3m"
        transformation_matrix = [2, 0, 0, 0, 2, 0, 0, 0, 2]
        cell_multiplicity = 32
        coordination_distances_mg = {"O": [2.12, 2.12, 2.12, 2.12, 2.12, 2.12]}
        coordination_distances_o = {"Mg": [2.12, 2.12, 2.12, 2.12, 2.12, 2.12]}
        mg1 = IrreducibleSite(irreducible_name="Mg1",
                              element="Mg",
                              first_index=1,
                              last_index=32,
                              representative_coords=[0.0, 0.0, 0.0],
                              wyckoff="a",
                              site_symmetry="m-3m",
                              coordination_distances=coordination_distances_mg)
        o1 = IrreducibleSite(irreducible_name="O1",
                             element="O",
                             first_index=33,
                             last_index=64,
                             representative_coords=[0.25, 0.25, 0.25],
                             wyckoff="b",
                             site_symmetry="m-3m",
                             coordination_distances=coordination_distances_o)
        irreducible_sites = [mg1, o1]
        dopant_configs = [["Al", "Mg"], ["Al", "O"], ["N", "Mg"], ["N", "O"]]
        antisite_configs = [["Mg", "O"], ["O", "Mg"]]
        interstitial_sites = ["i1"]
        complex_defect_names = []
        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        displacement_distance = 0.15
        cutoff = 2.0
        symprec = 0.001
        angle_tolerance = 5
        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}

        self.MgO = DefectInitialSetting(
            structure=self.structure,
            space_group_symbol=space_group_symbol,
            transformation_matrix=transformation_matrix,
            cell_multiplicity=cell_multiplicity,
            irreducible_sites=irreducible_sites,
            dopant_configs=dopant_configs,
            antisite_configs=antisite_configs,
            interstitial_site_names=interstitial_sites,
            complex_defect_names=complex_defect_names,
            included=included,
            excluded=excluded,
            displacement_distance=displacement_distance,
            cutoff=cutoff,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            oxidation_states=oxidation_states,
            electronegativity=electronegativity)

        d1 = {"name": "Va_O1",
              "defect_type": DefectType.vacancy,
              "initial_structure": self.structure,
              "removed_atoms": [{"element": "O",
                                 "index": 8,
                                 "coords": [0.25, 0.25, 0.25]}],
              "inserted_atoms": list(),
              "changes_of_num_elements": {"O": -1},
              "initial_site_symmetry": "m-3m",
              "charges": 2,
              "num_equiv_sites": 64,
              "center": [0, 0, 0]}
        self.defect_set = [d1]

    def test_dict(self):
        # roundtrip: object -> dict -> object
        mgo_dict = self.MgO.as_dict()
        mgo_from_dict = DefectInitialSetting.from_dict(mgo_dict)
        # Note: irreducible_sites usually return pointers, so __eq__ is
        # overloaded in DefectInitialSetting.
        print(mgo_from_dict.as_dict())
        print(self.MgO.as_dict())
        self.assertTrue(mgo_from_dict.as_dict() == self.MgO.as_dict())
        self.MgO.to()

    def test_to_json_file(self):
        # round trip test
        with tempfile.NamedTemporaryFile() as fp:
            tmp_json = fp.name
            self.MgO.to_json_file(tmp_json)
            mgo_from_json = DefectInitialSetting.load_json(tmp_json)
            self.assertTrue(mgo_from_json.as_dict() == self.MgO.as_dict())

    def test_msonalbe(self):
        self.assertMSONable(self.MgO)

    # def test_from_defect_in(self):
    #     mgo_from_defect_in = \
    #         DefectInitialSetting.from_defect_in(
    #             poscar=os.path.join(test_dir, "POSCAR-MgO64atoms"),
    #             defect_in_file=os.path.join(test_dir, "defect.in.example"))

        # self.assertTrue(mgo_from_defect_in == self.MgO)

    def test_from_basic_settings(self):
        mgo_from_basic_settings = \
            DefectInitialSetting.from_basic_settings(
                structure=self.structure,
                transformation_matrix=[2, 2, 2],
                cell_multiplicity=32,
                dopants=["Al", "N"],
                is_antisite=True,
                interstitial_sites=["i1"],
                complex_defect_names=None,
                en_diff=4.0,
                included=["Va_O1_-1", "Va_O1_-2"],
                excluded=["Va_O1_1", "Va_O1_2"],
                displacement_distance=0.15,
                cutoff=2.0,
                symprec=0.001)

        print(mgo_from_basic_settings.interstitials["i1"])
#        self.assertTrue(mgo_from_basic_settings == self.MgO)

    def test_make_defect_set(self):
        # Sequence of expected is changed for easy view. Thus, sort is needed
        # # for comparison.
        self.MgO.make_defect_set()
        print(self.MgO.defect_entries)
        # expected = \
        #     ['Va_Mg1_-2', 'Va_Mg1_-1', 'Va_Mg1_0', 'Va_O1_2', 'Mg_i1_0',
        #      'Mg_i1_1', 'Mg_i1_2', 'O_i1_-2', 'O_i1_-1', 'O_i1_0', 'N_i1_-3',
        #      'N_i1_-2', 'N_i1_-1', 'N_i1_0', 'N_i1_1', 'Al_i1_-1', 'Al_i1_0',
        #      'Al_i1_1', 'Al_i1_2', 'Al_i1_3', 'Mg_O1_0', 'Mg_O1_1', 'Mg_O1_2',
        #      'Mg_O1_3', 'Mg_O1_4', 'O_Mg1_-4', 'O_Mg1_-3', 'O_Mg1_-2',
        #      'O_Mg1_-1', 'O_Mg1_0', 'Al_Mg1_-1', 'Al_Mg1_0', 'Al_Mg1_1',
        #      'Al_O1_-1', 'Al_O1_0', 'Al_O1_1', 'Al_O1_2', 'Al_O1_3', 'Al_O1_4',
        #      'Al_O1_5', 'N_Mg1_-5', 'N_Mg1_-4', 'N_Mg1_-3', 'N_Mg1_-2',
        #      'N_Mg1_-1', 'N_Mg1_0', 'N_Mg1_1', 'N_O1_-1', 'N_O1_0', 'N_O1_1']
        #
#        actual = [str(i) for i in self.MgO.make_defect_set()]
        # self.assertEqual(sorted(actual), sorted(expected))

