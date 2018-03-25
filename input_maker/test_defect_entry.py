import unittest
from defect_entry import *
from pymatgen.io.vasp.inputs import Poscar

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

DEFAULT_POTCAR_DIR="/home/common/default_POTCAR"
TEST_DIRECTORY = "../examples/MgO"
DIRNAME_VAC = TEST_DIRECTORY + "/defects/Va_O1_2"


class DefectEntryTest(unittest.TestCase):

    def setUp(self):
        """ """
        initial_structure = Poscar.from_file(DIRNAME_VAC + "/POSCAR").structure
        removed_atom_index = 9
        inserted_atom_index = None
        defect_coords = [0.25, 0.25, 0.25]
        in_name = None
        out_name = "O1"
        charge = 2

        self._MgO_Va_O1_2 = DefectEntry(
            initial_structure, removed_atom_index, inserted_atom_index,
            defect_coords, in_name, out_name, charge)

    def test_Dict(self):
        self.dict = self._MgO_Va_O1_2.as_dict()
        self._MgO_Va_O1_2_from_dict = DefectEntry.from_dict(self.dict)


if __name__ == "__main__":
    unittest.main()

