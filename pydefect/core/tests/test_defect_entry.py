# -*- coding: utf-8 -*-
import tempfile
import os
import unittest
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

from pydefect.core.defect_entry import DefectType, DefectEntry, divide_dirname
from pydefect.core.config import CUTOFF_RADIUS

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectEntryTest(PymatgenTest):

    def setUp(self):
        """ """
        # DefectEntry class object for a single vacancy
        name = "Va_O1"
        defect_type = DefectType.from_string("vacancy")
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO8atoms-Va_O1"))
        perturbed_initial_structure = initial_structure.copy()
        removed_atoms = {8: [0.25, 0.25, 0.25]}
        inserted_atoms = {}
        changes_of_num_elements = {"O": -1}
        charge = 2
        initial_site_symmetry = "Oh"
        num_equiv_sites = 4
        neighboring_sites = [0, 1]
        self._MgO_Va_O1_2 = \
            DefectEntry(name=name,
                        defect_type=defect_type,
                        initial_structure=initial_structure,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        cutoff=CUTOFF_RADIUS,
                        neighboring_sites=neighboring_sites,
                        num_equiv_sites=num_equiv_sites)

        # DefectEntry class object for a complex defect
        name = "2Va_O1+Mg_i1"
        defect_type = DefectType.from_string("complex")
        initial_structure = Structure.from_file(
            os.path.join(test_dir, "POSCAR-MgO8atoms-2Va_O1-Mg_i1_2"))
        perturbed_initial_structure = initial_structure.copy()
        removed_atoms = {8: [0.25, 0.25, 0.25], 9: [0.25, 0.25, 0.75]}
        inserted_atoms = {8: [0.25, 0.25, 0.25]}
        changes_of_num_elements = {"O": -2, "Mg": 1}
        charge = 2
        initial_site_symmetry = "mmm"
        neighboring_sites = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
        num_equiv_sites = 24

        self._MgO_complex = \
            DefectEntry(name=name,
                        defect_type=defect_type,
                        initial_structure=initial_structure,
                        perturbed_initial_structure=perturbed_initial_structure,
                        removed_atoms=removed_atoms,
                        inserted_atoms=inserted_atoms,
                        changes_of_num_elements=changes_of_num_elements,
                        charge=charge,
                        initial_site_symmetry=initial_site_symmetry,
                        cutoff=CUTOFF_RADIUS,
                        neighboring_sites=neighboring_sites,
                        num_equiv_sites=num_equiv_sites)

    def test_from_yaml(self):
        defect_entry_from_yaml = DefectEntry.from_yaml(
            filename=os.path.join(test_dir, "defect_entry-2Va_O1-Mg_i1_2.yaml"),
            defect_name="2Va_O1+Mg_i1_2")
        print(defect_entry_from_yaml)
        print(self._MgO_complex)
        self.assertTrue(defect_entry_from_yaml.as_dict() == self._MgO_complex.as_dict())

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

    def test_dict_round_trip(self):
        """ round trip test of as_dict and from_dict """
        dict = self._MgO_Va_O1_2.as_dict()
        Va_O1_2_from_dict = DefectEntry.from_dict(dict)
        for i in dict.keys():
            self.assertTrue(Va_O1_2_from_dict.as_dict()[i] ==
                        self._MgO_Va_O1_2.as_dict()[i])

    def test_json_round_trip(self):
        """ round trip test of to_json and from_json """
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

    def test_defect_center(self):
        pos = [[0.25, 0.25, 0.25], [0.25, 0.25, -0.25], [0.25, 0.25, 0.25]]
        expected = list(np.average(np.array(pos), axis=0))

        actual = self._MgO_complex.defect_center_coords

        self.assertArrayAlmostEqual(actual, expected)

    def test_anchor_atom_index(self):
        expected = 14    # [0.75, 0.75, 0.75]
        actual = self._MgO_Va_O1_2.anchor_atom_index

        self.assertEqual(actual, expected)


class DivideDirnameTest(PymatgenTest):
    def setUp(self):
        """ """
        # DefectEntry class object for a single vacancy
        self.dirname1 = "Va_Mg1_2"
        self.dirname2 = "Va_O1_2_inward"
        self.dirname3 = "Mg_i+Va_O1*2_-2_coord1"

    def test_dirname1(self):
        self.assertEqual(divide_dirname(self.dirname1),
                         ("Va_Mg1", 2, None))

    def test_dirname2(self):
        self.assertEqual(divide_dirname(self.dirname2),
                         ("Va_O1", 2, "inward"))

    def test_dirname3(self):
        self.assertEqual(divide_dirname(self.dirname3),
                         ("Mg_i+Va_O1*2", -2, "coord1"))


if __name__ == "__main__":
    unittest.main()

