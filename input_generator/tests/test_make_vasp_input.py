import unittest
from pydefect.input_generator.make_vasp_input import *
from pydefect.input_generator.defect_input import IrrepElement, DefectSetting
import numpy as np

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

DEFAULT_POTCAR_DIR="/home/common/default_POTCAR"

class NormedRandom3DVectorTest(unittest.TestCase): 

    def setUp(self): 
        self.v = normed_random_3D_vector()

    def test_norm(self):
        print("normed_random_3D_vector: ", self.v)
        print("norm: ", np.linalg.norm(self.v))
        self.assertAlmostEqual(np.linalg.norm(self.v), 1.0) 


class RandomVectorTest(unittest.TestCase):

    def setUp(self): 
        self.distance = 3.0
        normed_v = normed_random_3D_vector()
        self.v = random_vector(normed_v, self.distance)

    def test_norm(self):
        print("random_3D_vector: ", self.v)
        print("distance: ", self.distance)
        print("norm: ", np.linalg.norm(self.v))
        self.assertLessEqual(np.linalg.norm(self.v), self.distance) 


class PotcarDirTest(unittest.TestCase):

    def test_potcar_dir(self):
        directory = potcar_dir()
        print("POTCAR directory: ", directory)


class MakePOTCARTest(unittest.TestCase):

    def setUp(self):
        self.dirname="."
        self.elements=["Mg", "O"]

    def test_make_POTCAR(self):
        make_POTCAR(self.dirname, self.elements, DEFAULT_POTCAR_DIR)


class DefectNameTest(unittest.TestCase):

    def setUp(self):
        self.a = "Va_Mg1_-2"
        self.b = "Mg_i2_2"
        self.c = "Mg_i_2"

    def test_DefectName(self):
        d_a = DefectName(self.a)
        d_b = DefectName(self.b)
        # return error message
#        d_c = DefectName(self.c)
        print(self.a, "-->", d_a.__dict__)
        print(self.b, "-->", d_b.__dict__)
#        print(self.c, "-->", d_c.__dict__)


class VaspInputMakerTest(unittest.TestCase):

    def setUp(self):
        self.vacancy = "Va_Mg1_-2"
        self.interstital1 = "Mg_i1_2"
        self.interstital2 = "O_i2_-2"
        self.interstital3 = "Al_i1_3"
        self.antisite1 = "Mg_O1_2"
        self.antisite1 = "O_Mg1_-2"

        structure = Structure.from_file("POSCAR-MgO64atoms")                       
        Mg1 = IrrepElement(irrepname="Mg1", element="Mg", first_index=1, last_index=32, repr_coord=[0, 0, 0]) 
        O1 = IrrepElement(irrepname="O1", element="O", first_index=33, last_index=64, repr_coord=[0.25, 0.25, 0.25])
        irrep_elements = [Mg1, O1]                                                 
        dopant_configs = "Al_Mg1"                                                  
        antisite_configs = ["Mg_O1", "O_Mg1"]                                      
        interstitial_coords = [[0.1, 0.1, 0.1], [0.3, 0.3, 0.3]]
        include = ["Va_O1_-1", "Va_O1_-2"]                                         
        exclude = ["Va_O1_1", "Va_O1_2"]                                           
        symbreak = True                                                            
        displace = 0.15                                                            
        cutoff = 3.0                                                               
        symprec = 0.001                                                            
        oxidation_states = {"Mg": 2, "O": -2}                                      
        electronegativity = {"Mg": 1.31, "O": 3.44}                                
                                                                                   
        self.defect_setting = \
            DefectSetting(structure, irrep_elements, dopant_configs,         
                          antisite_configs, interstitial_coords, include, 
                          exclude, symbreak, displace, cutoff, symprec, 
                          oxidation_states, electronegativity)  

    def test_vacancy(self):
        va = VaspInputMaker(self.defect_setting, self.vacancy)
        print(va.defect_structure)
        print(va.defect_position)

    def test_interstitial(self):
        in1 = VaspInputMaker(self.defect_setting, self.interstital1)
        in2 = VaspInputMaker(self.defect_setting, self.interstital2)
        in3 = VaspInputMaker(self.defect_setting, self.interstital3)
        print(in1.defect_structure)
        print(in2.defect_structure)
        print(in3.defect_structure)

#class PerturbAroundDefectTest(unittest.TestCase):
#
#    def setUp(self):
#        self.structure = Structure.from_file("POSCAR-MgO64atoms")
#
#    def test_atom_index(self):
##        defect_position = 1
#        defect_position = [0.1, 0.1, 0.1]
#        cutoff = 3.0
#        distance = 0.1
#        new_structure = perturb_around_defect(self.structure, defect_position, cutoff, distance)
#
#
#    def test_
#        
#
