# -*- coding: utf-8 -*-
import tempfile

import numpy as np
from pydefect.core.defect_entry import (
    DefectType, DefectEntry, determine_defect_type, anchor_atom_index,
    divide_defect_name)
from pydefect.util.testing import PydefectTest
from pydefect.core.config import CUTOFF_FACTOR


class DefectTypeTest(PydefectTest):

    def setUp(self) -> None:
        self.vacancy = DefectType.from_string("vacancy")
        self.substituted = DefectType.substituted

    def tests_str(self):
        self.assertEqual("vacancy", str(self.vacancy))

    def test__raise_error(self):
        with self.assertRaises(AttributeError):
            DefectType.from_string("antisite")

    def test_is_defect_center_atom(self):
        self.assertTrue(self.substituted.is_defect_center_atom)


class DetermineDefectTypeTest(PydefectTest):

    def setUp(self) -> None:
        self.vacancy = determine_defect_type(
            inserted_atoms=[], removed_atoms=[{"coords":[0, 0, 0]}])
        self.interstitial = determine_defect_type(
            inserted_atoms=[{"coords":[0, 0, 0]}], removed_atoms=[])
        self.substituted = determine_defect_type(
            inserted_atoms=[{"coords":[0, 0, 0]}],
            removed_atoms=[{"coords":[0, 0, 0]}])
        self.complex = determine_defect_type(
            inserted_atoms=[{"coords":[0, 0, 0.5]}],
            removed_atoms=[{"coords":[0, 0, 0]}])
        self.complex2 = determine_defect_type(
            inserted_atoms=[{"coords":[0, 0, 0]}],
            removed_atoms=[{"coords":[0, 0, 0]}, {"coords": [0, 0, 0.5]}])

    def test(self):
        self.assertEqual(DefectType.vacancy, self.vacancy)
        self.assertEqual(DefectType.interstitial, self.interstitial)
        self.assertEqual(DefectType.substituted, self.substituted)
        self.assertEqual(DefectType.complex, self.complex)
        self.assertEqual(DefectType.complex, self.complex2)

    def test_raise_error(self):
        with self.assertRaises(ValueError):
            determine_defect_type(inserted_atoms=[], removed_atoms=[])


class DefectEntryTest(PydefectTest):

    def setUp(self):
        """ """
        # DefectEntry class object for a single vacancy
        name = "Va_O1"
        defect_type = DefectType.from_string("vacancy")
        self.initial_structure_vacancy = \
            self.get_structure_by_name(name="MgO64atoms-Va_O_0-unrelaxed")
        perturbed_initial_structure = self.initial_structure_vacancy.copy()
        removed_atoms = [{"element": "O",
                          "index": 32,
                          "coords": [0.25, 0, 0]}]
        inserted_atoms = []
        changes_of_num_elements = {"O": -1}
        charge = 2
        initial_site_symmetry = "m-3m"
        multiplicity = 32
        neighboring_sites = [0, 4, 16, 17, 24, 26]
        cutoff = round(8.419456 / 4 * CUTOFF_FACTOR, 2)
        self.MgO_Va_O1_2 = \
            DefectEntry(name=name,
                        defect_type=defect_type,
                        initial_structure=self.initial_structure_vacancy,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        cutoff=cutoff,
                        neighboring_sites=neighboring_sites,
                        multiplicity=multiplicity)

        name = "complex"
        defect_type = DefectType.from_string("complex")
        self.initial_structure_complex = \
            self.get_structure_by_name(name="MgO64atoms-Va_Mg+Va_O+Ca_i")
        perturbed_initial_structure = self.initial_structure_complex.copy()
        removed_atoms = [{"element": "Mg",
                          "index": 0,
                          "coords": [0, 0, 0]},
                         {"element": "O",
                          "index": 32,
                          "coords": [0.25, 0, 0]}]
        inserted_atoms = [{"element": "Ca",
                           "index": 0,
                           "coords": [0.125, 0, 0]}]
        changes_of_num_elements = {"Mg": -1, "Ca": 1, "O": -1}
        charge = 2
        initial_site_symmetry = "4mm"
        multiplicity = 192
        neighboring_sites = [16, 17, 24, 26, 47, 48, 55, 57]
        cutoff = round(8.419456 / 4 * CUTOFF_FACTOR, 2)
        self.MgO_complex = \
            DefectEntry(name=name,
                        defect_type=defect_type,
                        initial_structure=self.initial_structure_complex,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        cutoff=cutoff,
                        neighboring_sites=neighboring_sites,
                        multiplicity=multiplicity)

    def test_from_defect_structure(self):
        perfect_structure = self.get_structure_by_name(name="MgO64atoms")
        expected = self.MgO_Va_O1_2.as_dict()
        actual = DefectEntry.from_defect_structure(
            defect_structure=self.initial_structure_vacancy,
            perfect_structure=perfect_structure,
            defect_name="Va_O1_2").as_dict()
        for d in expected:
            self.assertEqual(expected[d], actual[d])

    def test_from_defect_structure_complex(self):
        perfect_structure = self.get_structure_by_name(name="MgO64atoms")
        expected = self.MgO_complex.as_dict()
        actual = DefectEntry.from_defect_structure(
            defect_structure=self.initial_structure_complex,
            perfect_structure=perfect_structure,
            defect_name="complex_2").as_dict()
        for d in expected:
            self.assertEqual(expected[d], actual[d])

    def test_msonable(self):
        self.assertMSONable(self.MgO_Va_O1_2)
        self.assertMSONable(self.MgO_complex)

    def test_dict_round_trip(self):
        """ round trip test of as_dict and from_dict """
        expected = self.MgO_Va_O1_2.as_dict()
        actual = DefectEntry.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)

    def test_json_round_trip(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.MgO_Va_O1_2.to_json_file(tmp_file.name)
        expected = self.MgO_Va_O1_2.as_dict()
        actual = DefectEntry.load_json(tmp_file.name).as_dict()
        self.assertTrue(expected, actual)

    def test_atom_mapping_to_perfect(self):
        expected = list(range(64))
        expected.pop(32)
        actual = self.MgO_Va_O1_2.atom_mapping_to_perfect
        self.assertEqual(expected, actual)

    def test_atom_mapping_to_perfect_complex(self):
        expected = list(range(64))
        expected.pop(32)
        expected.pop(0)
        expected = [None] + expected
        actual = self.MgO_complex.atom_mapping_to_perfect
        self.assertEqual(expected, actual)

    def test_defect_center(self):
        expected = [0.25, 0, 0]
        actual = self.MgO_Va_O1_2.defect_center_coords
        self.assertArrayEqual(expected, actual)

    def test_defect_center_complex(self):
        expected = [0.125, 0, 0]
        actual = self.MgO_complex.defect_center_coords
        self.assertArrayEqual(expected, actual)

    def test_anchor_atom_index(self):
        expected = 38    # Line 47 [0.75, 0.5, 0.5]
        actual = self.MgO_Va_O1_2.anchor_atom_index
        self.assertEqual(actual, expected)

    def test_anchor_atom_index_complex(self):
        expected = [8, 38] # 8 or 38
        actual = self.MgO_complex.anchor_atom_index
        self.assertTrue(actual in expected)


class AnchorAtomIndexTest(PydefectTest):
    def test(self):
        structure = self.get_structure_by_name(name="KZn4P3")
        actual = anchor_atom_index(structure, [0.5, 0.5, 0.5])
        self.assertEqual(7, actual)


class DivideDirnameTest(PydefectTest):
    def setUp(self):
        # DefectEntry class object for a single vacancy
        self.dirname1 = "Va_Mg1_2"
        self.dirname2 = "Va_O1_2_inward"
        self.dirname3 = "Mg_i+Va_O1*2_-2_coord1"

    def test_dirname1(self):
        self.assertEqual(("Va_Mg1", 2, None), divide_defect_name(self.dirname1))

    def test_dirname2(self):
        self.assertEqual(("Va_O1", 2, "inward"),
                         divide_defect_name(self.dirname2))

    def test_dirname3(self):
        self.assertEqual(("Mg_i+Va_O1*2", -2, "coord1"),
                         divide_defect_name(self.dirname3))


