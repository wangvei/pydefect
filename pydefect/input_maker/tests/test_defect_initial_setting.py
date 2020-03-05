# -*- coding: utf-8 -*-
import tempfile
from copy import deepcopy
from pathlib import Path

from pydefect.core.defect_entry import DefectType
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.input_maker.defect_initial_setting import (
    default_charge_set, get_electronegativity, get_oxidation_state,
    dopant_info, get_distances_from_string, insert_atoms, select_defects,
    DefectInitialSetting)
from pydefect.util.testing import PydefectTest
from pydefect.core.defect_name import DefectName

from pymatgen.core.structure import Structure

parent = Path(__file__).parent


class CandidateChargeSetTest(PydefectTest):
    def test_range1(self):
        actual_positive = default_charge_set(2)
        expected_positive = {0, 1, 2}
        self.assertEqual(expected_positive, actual_positive)

    def test_range2(self):
        actual_negative = default_charge_set(-2)
        expected_negative = {-2, -1, 0}
        self.assertEqual(expected_negative, actual_negative)

    def test_range3(self):
        actual_positive = default_charge_set(3)
        expected_positive = {-1, 0, 1, 2, 3}
        self.assertEqual(expected_positive, actual_positive)

    def test_range4(self):
        actual_negative = default_charge_set(-3)
        expected_negative = {-3, -2, -1, 0, 1}
        self.assertEqual(expected_negative, actual_negative)


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
        self.name_set = {"Va_Mg1": {"charges": {-1, 0, 1, 2}},
                         "Va_O1": {"charges": {0}},
                         "O_i1": {"charges": {0, 1}}}

    def test_keywords1(self):
        actual_va = select_defects(self.name_set, keywords=["Va"])
        expected_va = {"Va_Mg1": {"charges": {-1, 0, 1, 2}},
                       "Va_O1": {"charges": {0}},
                       "O_i1": {"charges": set()}}
        self.assertEqual(expected_va, actual_va)

    def test_keywords2(self):
        actual_va = select_defects(self.name_set, keywords=["Mg_i1"])
        expected_va = {"Va_Mg1": {"charges": set()},
                       "Va_O1": {"charges": set()},
                       "O_i1": {"charges": set()}}
        self.assertEqual(expected_va, actual_va)

    def test_specified(self):
        actual_specified = \
            select_defects(self.name_set,
                           specified_defects=["Va_Mg1_1", "Va_Mg1_3"])
        print(actual_specified)
        expected_specified = {"Va_Mg1": {"charges": {1, 3}},
                              "Va_O1": {"charges": set()},
                              "O_i1": {"charges": set()}}
        self.assertEqual(expected_specified, actual_specified)

    def test_included_excluded(self):
        actual = select_defects(self.name_set,
                                included=["Va_Mg1_3"], excluded=["Va_O1_0"])
        print(actual)
        expected = {"Va_Mg1": {"charges": {-1, 0, 1, 2, 3}},
                    "Va_O1": {"charges": set()},
                    "O_i1": {"charges": {0, 1}}}
        self.assertEqual(expected, actual)


class DefectInitialSettingTest(PydefectTest):

    def setUp(self):
        """ MgO 64atoms"""
        self.structure = self.get_structure_by_name("defect_initial_setting")
        space_group_symbol = "Fm-3m"
        transformation_matrix = [-2, 2, 2, 2, -2, 2, 2, 2, -2]
        cell_multiplicity = 32
        coordination_distances_mg = {"O": [2.1] * 6}
        coordination_distances_o = {"Mg": [2.1] * 6}
        mg1 = IrreducibleSite(irreducible_name="Mg1",
                              element="Mg",
                              first_index=0,
                              last_index=31,
                              representative_coords=[0.5, 0.0, 0.0],
                              wyckoff="a",
                              site_symmetry="m-3m",
                              cutoff=2.74,
                              coordination_distances=coordination_distances_mg)
        o1 = IrreducibleSite(irreducible_name="O1",
                             element="O",
                             first_index=32,
                             last_index=63,
                             representative_coords=[0.75, 0.25, 0.25],
                             wyckoff="b",
                             site_symmetry="m-3m",
                             cutoff=2.74,
                             coordination_distances=coordination_distances_o)
        irreducible_sites = [mg1, o1]
        dopant_configs = [["Al", "Mg"], ["N", "O"]]
        antisite_configs = []
        interstitial_sites = ["i1"]
        complex_defect_names = ["divacancy"]
        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        displacement_distance = 0.2
        symprec = 0.01
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
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            oxidation_states=oxidation_states,
            electronegativity=electronegativity,
            interstitials_yaml=str(parent / "interstitials.yaml"),
            complex_defect_yaml=str(parent / "complex_defects.yaml"))

    def test_dict(self):
        expected = self.MgO.as_dict()
        actual = DefectInitialSetting.from_dict(expected).as_dict()
        # Note: irreducible_sites usually return pointers, so __eq__ is
        # overloaded in DefectInitialSetting.
        self.assertEqual(expected, actual)
#        self.MgO.to()

    def test_to_json_file(self):
        expected = self.MgO.as_dict()
        with tempfile.NamedTemporaryFile() as fp:
            tmp_json = fp.name
            self.MgO.to_json_file(tmp_json)
            actual = DefectInitialSetting.load_json(tmp_json).as_dict()
            self.assertEqual(expected, actual)

    def test_msonalbe(self):
        self.assertMSONable(self.MgO)

    def test_from_defect_in(self):
        actual = DefectInitialSetting.from_defect_in(
            poscar=str(parent / "DPOSCAR-defect_initial_setting"),
            defect_in_file=str(parent / "defect_unittest.in"),
            interstitials_yaml=str(parent / "interstitials.yaml"),
            complex_defect_yaml=str(parent / "complex_defects.yaml"))

        self.assertEqual(self.MgO.structure.lattice, actual.structure.lattice)

        actual = actual.as_dict()

        for key, expected in self.MgO.as_dict().items():
            if key not in \
                    ["structure", "interstitials_yaml", "complex_defect_yaml"]:
                self.assertEqual(expected, actual[key])

    def test_from_basic_settings(self):
        actual = DefectInitialSetting.from_basic_settings(
            structure=self.structure,
            transformation_matrix=[-2, 2, 2, 2, -2, 2, 2, 2, -2],
            cell_multiplicity=32,
            dopants=["Al", "N"],
            is_antisite=True,
            interstitial_sites=["i1"],
            complex_defect_names=["divacancy"],
            en_diff=1.0,
            included=["Va_O1_-1", "Va_O1_-2"],
            excluded=["Va_O1_1", "Va_O1_2"],
            displacement_distance=0.2,
            symprec=0.01,
            interstitials_yaml=str(parent / "interstitials.yaml"),
            complex_defect_yaml=str(parent / "complex_defects.yaml"))

        actual.to(defect_in_file=str(parent / "defect_unittest.in"),
                  poscar_file=str(parent / "DPOSCAR-defect_initial_setting"))
        self.assertEqual(self.MgO.structure, actual.structure)

        actual = actual.as_dict()

        for key, expected in self.MgO.as_dict().items():
            if key not in \
                    ["structure", "interstitials_yaml", "complex_defect_yaml"]:
                self.assertEqual(expected, actual[key])

    def test_make_all_defect_set(self):
        mgo_all = deepcopy(self.MgO)
        mgo_all.make_defect_set()

        # vacancies: 3 + 3 = 6
        # dopants: 3 + 3 = 6
        # interstitials: 3(Mg) + 3(O) + 5(Al) + 5(N) = 16
        # divacancies: 1
        self.assertEqual(29, len(mgo_all.defect_entries))

    def test_make_defect_set_keywords(self):
        mgo_interstitial = deepcopy(self.MgO)
        mgo_interstitial.make_defect_set(keywords=["i1"])
        actual = {str(DefectName(n.name, n.charge))
                  for n in mgo_interstitial.defect_entries}
        expected = {'Mg_i1_0', 'Mg_i1_1', 'Mg_i1_2', 'O_i1_0', 'O_i1_-2',
                    'O_i1_-1', 'Al_i1_0', 'Al_i1_1', 'Al_i1_2', 'Al_i1_3',
                    'Al_i1_-1', 'N_i1_0', 'N_i1_1', 'N_i1_-1', 'N_i1_-3',
                    'N_i1_-2'}
        self.assertEqual(expected, actual)

    def test_make_specified_defect_set(self):
        mgo_specified = deepcopy(self.MgO)
        mgo_specified.make_defect_set(
            specified_defects=["Al_Mg1_5", "N_O1_2"])
        actual = {str(DefectName(n.name, n.charge))
                  for n in mgo_specified.defect_entries}
        expected = {"Al_Mg1_5", "N_O1_2"}
        self.assertEqual(expected, actual)

        actual = mgo_specified.defect_entries[0].as_dict()

        structure = Structure.from_file(
            parent / "DPOSCAR-defect_initial_setting")
        structure.replace(0, "Al")
        expected = {"name": "Al_Mg1",
                    "defect_type": str(DefectType.substituted),
                    "initial_structure": structure,
                    "perturbed_initial_structure": structure,
                    "removed_atoms": [{"element": "Mg",
                                       "index": 0,
                                       "coords": [0.5, 0.0, 0.0]}],
                    "inserted_atoms": [{"element": "Al",
                                        "index": 0,
                                        "coords": [0.5, 0.0, 0.0]}],
                    "changes_of_num_elements": {"Al": 1, "Mg": -1},
                    "charge": 5,
                    "initial_site_symmetry": "m-3m",
                    "cutoff": 2.74,
                    "neighboring_sites": [46, 54, 59, 61, 62, 63],
                    "annotation": None,
                    "multiplicity": 32}

#        self.assertEqual(expected, actual)
        for key, value in expected.items():
            if key == "perturbed_initial_structure":
                continue
            if key == "neighboring_sites":
                self.assertEqual(set(value), set(actual[key]))
            else:
                self.assertEqual(value, actual[key])


# class DefectInitialSettingTest2(PydefectTest):

    # def setUp(self):
    #     """ Al2O3 """
    #     self.structure = PydefectTest.get_structure_by_name("Al2O3")

