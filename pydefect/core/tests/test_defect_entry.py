# -*- coding: utf-8 -*-

import tempfile
import os
import unittest
from copy import deepcopy

from pymatgen.core.structure import Structure

from pydefect.core.defect_entry import get_num_atoms_for_elements, \
    DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class GetNumAtomsForElementsTest(unittest.TestCase):

    def test(self):
        structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO64atoms-Al_doped"))

        expected = [1, 31, 32]
        actual = get_num_atoms_for_elements(structure)
        self.assertEqual(actual, expected)


class DefectEntryTest(unittest.TestCase):

    def setUp(self):
        """ """
        # DefectEntry class object for a single vacancy
        name = "Va_O1"
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO8atoms-Va_O1"))
        perturbed_initial_structure = initial_structure.copy()
        removed_atoms = {8: [0.25, 0.25, 0.25]}
        inserted_atoms = []
        changes_of_num_elements = {"O": -1}
        charge = 2
        initial_site_symmetry = "Oh"
        num_equiv_sites = 4
        perturbed_sites = []
        self._MgO_Va_O1_2 = \
            DefectEntry(name=name,
                        initial_structure=initial_structure,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        perturbed_sites=perturbed_sites,
                        num_equiv_sites=num_equiv_sites)

        # DefectEntry class object for a complex defect
        name = "2Va_O1+Mg_i1_2"
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO8atoms-2Va_O1-Mg_i1_2"))
        perturbed_initial_structure = initial_structure.copy()
        removed_atoms = {8: [0.25, 0.25, 0.25], 9: [0.25, 0.25, 0.75]}
        inserted_atoms = [8]
        changes_of_num_elements = {"O": -2, "Mg": 1}
        charge = 2
        initial_site_symmetry = "mmm"
        perturbed_sites = []
        num_equiv_sites = 0

        self._MgO_complex = \
            DefectEntry(name=name,
                        initial_structure=initial_structure,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        perturbed_sites=perturbed_sites,
                        num_equiv_sites=num_equiv_sites)

    def test_from_yaml(self):
        defect_entry_from_yaml = DefectEntry.from_yaml(
            os.path.join(test_dir, "defect_entry-2Va_O1-Mg_i1_2.yaml"))
        print(defect_entry_from_yaml)
        print(self._MgO_complex)
        self.assertTrue(
            defect_entry_from_yaml.as_dict() == self._MgO_complex.as_dict())

    # def test_from_simpler_yaml(self):
    #     simpler_dir = os.path.join(test_dir, "MgO/defects/2Va_O1-Mg_i1_2")
    #     os.chdir(simpler_dir)
    #     defect_entry_from_simpler_yaml = \
    #         DefectEntry.from_yaml("defect_entry.yaml")
    #     print(defect_entry_from_simpler_yaml)
#        self.assertTrue(defect_entry_from_simpler_yaml == self._MgO_complex)

    # def test_from_yaml_fail(self):
    #     with self.assertRaises(Exception) as context:
    #         DefectEntry.from_yaml(
    #             os.path.join(test_dir, "defect_entry-2Va_O1-Mg_i1_2_fail.yaml"))

            # self.assertTrue('This is broken' in context.exception)

    def test_dict_roundtrip(self):
        """ round trip test of as_dict and from_dict
        """
        Va_O1_2_from_dict = DefectEntry.from_dict(self._MgO_Va_O1_2.as_dict())
        self.assertTrue(Va_O1_2_from_dict.as_dict() ==
                        self._MgO_Va_O1_2.as_dict())

    def test_json(self):
        """ round trip test of to_json and from_json
        """
        tmp_file = tempfile.NamedTemporaryFile()
        self._MgO_Va_O1_2.to_json_file(tmp_file.name)
        defect_entry_from_json = DefectEntry.load_json(tmp_file.name)
        self.assertTrue(defect_entry_from_json.as_dict() ==
                        self._MgO_Va_O1_2.as_dict())

    def test_atom_mapping_to_perfect(self):
        expected = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
        actual = self._MgO_Va_O1_2.atom_mapping_to_perfect
        self.assertTrue(actual == expected)

        expected = [0, 1, 2, 3, 4, 5, 6, 7, None, 10, 11, 12, 13, 14, 15]
        actual = self._MgO_complex.atom_mapping_to_perfect
        self.assertTrue(actual == expected)


if __name__ == "__main__":
    unittest.main()

