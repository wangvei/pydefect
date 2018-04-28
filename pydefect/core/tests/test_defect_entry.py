# -*- coding: utf-8 -*-

import tempfile
import os
import unittest

from pymatgen.core.structure import Structure

from pydefect.core.defect_entry import get_num_atoms_for_elements, \
    element_diff_from_poscar_files, DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class FunctionTest(unittest.TestCase):

    def test_get_num_atoms_for_elements(self):
        structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO64atoms-Al_doped"))

        expected = [1, 31, 32]
        actual = get_num_atoms_for_elements(structure)
        self.assertEqual(actual, expected)

    def test_element_diff_from_poscar_files(self):
        poscar1 = os.path.join(test_dir, "POSCAR-MgO8atoms")
        poscar2 = os.path.join(test_dir, "POSCAR-MgO8atoms-Va_O1+N_O")

        actual = element_diff_from_poscar_files(poscar1, poscar2)
        expected = {"O": 2, "N": -1}
        self.assertEqual(actual, expected)


class DefectEntryTest(unittest.TestCase):

    def setUp(self):
        """ """
        # DefectEntry class object for a single vacancy
        name = "Va_O1"
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO8atoms-Va_O1"))
        removed_atoms = {8: [0.25, 0.25, 0.25]}
        inserted_atoms = []
        element_diff = {"O": -1}
        charge = 2
        self._MgO_Va_O1_2 = \
            DefectEntry(name, initial_structure, removed_atoms, inserted_atoms,
                        element_diff, charge)

        # DefectEntry class object for a complex defect
        name = "2Va_O1+Mg_i1_2"
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO8atoms-2Va_O1-Mg_i1_2"))
        removed_atoms = {8: [0.25, 0.25, 0.25], 9: [0.25, 0.25, 0.75]}
        inserted_atoms = [8]
        element_diff = {"O": -2, "Mg": 1}
        charge = 2

        self._MgO_complex = \
            DefectEntry(name, initial_structure, removed_atoms, inserted_atoms,
                        element_diff, charge)

    def test_dict(self):
        # object -> dict -> object
        d = self._MgO_Va_O1_2.as_dict()
        defect_entry_from_dict = DefectEntry.from_dict(d)
        self.assertTrue(defect_entry_from_dict == self._MgO_Va_O1_2)

    def test_from_yaml(self):
        defect_entry_from_yaml = DefectEntry.from_yaml(
            os.path.join(test_dir, "defect_entry-2Va_O1-Mg_i1_2.yaml"))
        self.assertTrue(defect_entry_from_yaml == self._MgO_complex)

    def test_from_simpler_yaml(self):
        simpler_dir = os.path.join(test_dir, "MgO/defects/2Va_O1-Mg_i1_2")
        os.chdir(simpler_dir)
        defect_entry_from_simpler_yaml = \
            DefectEntry.from_yaml("defect_entry.yaml")
        print(defect_entry_from_simpler_yaml)
#        self.assertTrue(defect_entry_from_simpler_yaml == self._MgO_complex)

    def test_from_yaml_fail(self):
        with self.assertRaises(Exception) as context:
            DefectEntry.from_yaml(
                os.path.join(test_dir, "defect_entry-2Va_O1-Mg_i1_2_fail.yaml"))

            self.assertTrue('This is broken' in context.exception)

    def test_json(self):
        # object -> json file -> object
        tmp_file = tempfile.NamedTemporaryFile()
        self._MgO_Va_O1_2.to_json_file(tmp_file.name)
        defect_entry_from_json = DefectEntry.load_json(tmp_file.name)
        self.assertTrue(defect_entry_from_json == self._MgO_Va_O1_2)

    def test_atom_mapping_to_perfect(self):
        expected = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
        actual = self._MgO_Va_O1_2.atom_mapping_to_perfect
        self.assertTrue(actual == expected)

#        expected = [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, None]
#        actual = self._MgO_complex.atom_mapping_to_perfect
#        self.assertTrue(actual == expected)


if __name__ == "__main__":
    unittest.main()

