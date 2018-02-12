import unittest
from defect import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

FILENAME_TO_JSON_FILE_VAC = "examples/Va_Mg1_-2.json"
FILENAME_TO_JSON_FILE_INT = "examples/Mg_i1_1.json"

class DefectTest(unittest.TestCase):

    def setUp(self):
        self._vac = Defect(removed_atom_index=1, inserted_atom_index=None,
                           defect_coords=[0, 0, 0],
                           in_name="Va", out_name="Mg1", charge=-2)
        self._int = Defect(removed_atom_index=None, inserted_atom_index=1,
                           defect_coords=[0.5, 0.5, 0.5],
                           in_name="Mg", out_name="i1", charge=1)

    def test_dict(self):
        d_vac = Defect.from_dict(self._vac.as_dict())
        d_int = Defect.from_dict(self._int.as_dict())
        self.assertTrue(vars(self._vac) == vars(d_vac))
        self.assertTrue(vars(self._int) == vars(d_int))

    def test_to_json_file(self):
        self._vac.to_json_file(FILENAME_TO_JSON_FILE_VAC)
        jf_vac = Defect.json_load(FILENAME_TO_JSON_FILE_VAC)
        self.assertTrue(vars(self._vac) == vars(jf_vac))

        self._int.to_json_file(FILENAME_TO_JSON_FILE_INT)
        jf_int = Defect.json_load(FILENAME_TO_JSON_FILE_INT)
        self.assertTrue(vars(self._int) == vars(jf_int))


class IrreducibleSiteTest(unittest.TestCase):
    
    def setUp(self):
        self._mg = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                                   first_index=1, last_index=32,
                                   repr_coords=[0, 0, 0])

    def test_dict(self):
        d = IrreducibleSite.from_dict(self._mg.as_dict())
        self.assertTrue(self._mg == d)

    def test_natoms(self):
        self.assertEqual(self._mg.natoms, 32)

