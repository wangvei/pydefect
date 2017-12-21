import unittest
from pydefect.input_generator.make_vasp_input import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

class PerturbAroundDefectTest(unittest.TestCase):

    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgO64atoms")

    def test_atom_index(self):
#        defect_position = 1
        defect_position = [0.1, 0.1, 0.1]
        cutoff = 3.0
        distance = 0.1
        new_structure = perturb_around_defect(self.structure, defect_position, cutoff, distance)

class VaspInputMakerTest(unittest.TestCase):

    def test_
        

