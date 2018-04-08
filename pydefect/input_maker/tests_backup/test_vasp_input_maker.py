import unittest
from pydefect.input_maker.make_vasp_input import *
from pydefect.input_maker.defect_input import IrrepElement, DefectSetting
import numpy as np

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

DEFAULT_POTCAR_DIR="/home/common/default_POTCAR"


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

class VaspInputMakerTest(unittest.TestCase):

    def setUp(self):
        self.vacancy = "Va_Mg1_-2"
        self.interstital1 = "Mg_i1_2"
        self.interstital2 = "O_i2_-2"
        self.interstital3 = "Al_i1_3"
        self.antisite1 = "Mg_O1_2"
        self.antisite2 = "O_Mg1_-2"

        structure = Structure.from_file("POSCAR-MgO64atoms")                       
        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg", first_index=1, last_index=32, repr_coords=[0, 0, 0]) 
        O1 = IrreducibleSite(irreducible_name="O1", element="O", first_index=33, last_index=64, repr_coords=[0.25, 0.25, 0.25])
        irreducible_sites = [Mg1, O1]                                                 
        dopant_configs = "Al_Mg"                                                  
        antisite_configs = ["Mg_O", "O_Mg"]                                      
        interstitial_coords = [[0.1, 0.1, 0.1], [0.3, 0.3, 0.3]]
        included = ["Va_O1_-1", "Va_O1_-2"]                                         
        excluded = ["Va_O1_1", "Va_O1_2"]                                           
        symbreak = True                                                            
        displace = 0.15                                                            
        cutoff = 3.0                                                               
        symprec = 0.001                                                            
        oxidation_states = {"Mg": 2, "O": -2}                                      
        electronegativity = {"Mg": 1.31, "O": 3.44}                                
                                                                                   
        self.defect_setting = \
            DefectSetting(structure, irreducible_sites, dopant_configs,         
                          antisite_configs, interstitial_coords, included, 
                          excluded, displace, cutoff, symprec, 
                          oxidation_states, electronegativity)  

    def test_vacancy(self):
        va = InputMaker(self.vacancy, self.defect_setting)
        va.analyze_name()
        print(va.defect_structure)
        print(va.defect_index)
        print(va.defect_coords)

    def test_interstitial(self):
        in1 = InputMaker(self.interstital1, self.defect_setting)
        in2 = InputMaker(self.interstital2, self.defect_setting)
        in3 = InputMaker(self.interstital3, self.defect_setting)

        in1.analyze_name()
        in2.analyze_name()
        in3.analyze_name()

#class VaspInputSetMakerTest(unittest.TestCase):
#
#    def setUp(self):
#
#        structure = Structure.from_file("POSCAR-MgO64atoms")                       
#        Mg1 = IrrepElement(irrepname="Mg1", element="Mg", first_index=1, last_index=32, repr_coord=[0, 0, 0]) 
#        O1 = IrrepElement(irrepname="O1", element="O", first_index=33, last_index=64, repr_coord=[0.25, 0.25, 0.25])
#        irrep_elements = [Mg1, O1]                                                 
#        dopant_configs = ["Al_Mg"]                                                 
#        antisite_configs = ["Mg_O", "O_Mg"]
#        interstitial_coords = [[0.1, 0.1, 0.1], [0.3, 0.3, 0.3]]
#        include = ["Va_O1_-1", "Va_O1_-2"]                                         
#        exclude = ["Va_O1_1", "Va_O1_2"]                                           
#        symbreak = True                                                            
#        displace = 0.15                                                            
#        cutoff = 3.0                                                               
#        symprec = 0.001                                                            
#        oxidation_states = {"Al": 3, "Mg": 2, "O": -2}                                      
#        electronegativity = {"Al": 1.3, "Mg": 1.31, "O": 3.44}                                
#                                                                                   
#        self.defect_setting = \
#            DefectInitialSetting(structure, irrep_elements, dopant_configs,
#                          antisite_configs, interstitial_coords, include, 
#                          exclude, symbreak, displace, cutoff, symprec, 
#                          oxidation_states, electronegativity)  
#
#    def test(self):        
#        a = VaspInputSetMaker(self.defect_setting)
#        a.constructor()
