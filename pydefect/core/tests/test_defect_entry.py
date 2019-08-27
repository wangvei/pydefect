# -*- coding: utf-8 -*-
import tempfile

import numpy as np
from pydefect.core.defect_entry import (
    DefectType, DefectEntry, determine_defect_type, divide_dirname)
from pydefect.util.testing import PydefectTest
from pydefect.core.config import CUTOFF_FACTOR

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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
        self.file_dir = self.TEST_FILES_DIR / "files"

        name = "Va_O1"
        defect_type = DefectType.from_string("vacancy")
        initial_structure = \
            self.get_structure_by_name(name="MgO64atoms-Va_O_0-unrelaxed")
        perturbed_initial_structure = initial_structure.copy()
        removed_atoms = [{"element": "O",
                          "index": 32,
                          "coords": [0.25, 0, 0]}]
        inserted_atoms = []
        changes_of_num_elements = {"O": -1}
        charge = 2
        initial_site_symmetry = "Oh"
        multiplicity = 4
        neighboring_sites = [0, 1]
        cutoff = round(8.419456 / 4 * CUTOFF_FACTOR, 2)
        self.MgO_Va_O1_2 = \
            DefectEntry(name=name,
                        defect_type=defect_type,
                        initial_structure=initial_structure,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        cutoff=cutoff,
                        neighboring_sites=neighboring_sites,
                        multiplicity=multiplicity)

    def test_msonable(self):
        self.assertMSONable(self.MgO_Va_O1_2)

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
#        expected = [0, 1, 2, 3, 4, 5, 6, 7, None, 10, 11, 12, 13, 14, 15]
#        actual = self.MgO_complex.atom_mapping_to_perfect
#        self.assertEqual(expected, actual)

    def test_defect_center(self):
        expected = [0.25, 0, 0]
        actual = self.MgO_Va_O1_2.defect_center_coords
        self.assertArrayEqual(actual, expected)

        # pos = [[0.25, 0.25, 0.25], [0.25, 0.25, -0.25], [0.25, 0.25, 0.25]]
        # expected = list(np.average(np.array(pos), axis=0))
        # actual = self.MgO_complex.defect_center_coords
        # self.assertArrayAlmostEqual(actual, expected)

    def test_anchor_atom_index(self):
        expected = 38    # Line 47 [0.75, 0.5, 0.5]
        actual = self.MgO_Va_O1_2.anchor_atom_index
        self.assertEqual(actual, expected)


class DivideDirnameTest(PydefectTest):
    def setUp(self):
        """ """
        # DefectEntry class object for a single vacancy
        self.dirname1 = "Va_Mg1_2"
        self.dirname2 = "Va_O1_2_inward"
        self.dirname3 = "Mg_i+Va_O1*2_-2_coord1"

    def test_dirname1(self):
        self.assertEqual(("Va_Mg1", 2, None), divide_dirname(self.dirname1))

    def test_dirname2(self):
        self.assertEqual(("Va_O1", 2, "inward"), divide_dirname(self.dirname2))

    def test_dirname3(self):
        self.assertEqual(("Mg_i+Va_O1*2", -2, "coord1"),
                         divide_dirname(self.dirname3))



