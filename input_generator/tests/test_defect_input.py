import unittest
from pydefect.input_generator.defect_input import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

class IrrepElementTest(unittest.TestCase):
    
    def setUp(self):
        self._mg = IrrepElement(irrepname="Mg1", element="Mg", first_index=1, 
                              last_index=32, repr_coord=[0, 0, 0])

    def test_dict(self):
        d = IrrepElement.from_dict(self._mg.as_dict())
        self.assertTrue(self._mg == d)

    def test_multiplicity(self):
        self.assertEqual(self._mg.multiplicity, 32)

class DefectSettingTest(unittest.TestCase):

    def setUp(self):
        structure = Structure.from_file("POSCAR-MgO64atoms")
        Mg1 = IrrepElement(irrepname="Mg1", element="Mg", first_index=1, last_index=32, repr_coord=[0, 0, 0])
        O1 = IrrepElement(irrepname="O1", element="O", first_index=33, last_index=64, repr_coord=[0.25, 0.25, 0.25])
        irrep_elements = [Mg1, O1] 
        dopant_configs = "Al_Mg1"
        antisite_configs = ["Mg_O1", "O_Mg1"]
        interstitial_coords = [[0.1, 0.1, 0.1]]
        include = ["Va_O1_-1", "Va_O1_-2"]
        exclude = ["Va_O1_1", "Va_O1_2"]
        symbreak = True
        displace = 0.15
        symprec = 0.001
        oxidation_states = {"Mg": 2, "O": -2}
        electron_negativity = {"Mg": 1.31, "O": 3.44}

        self._a = DefectSetting(structure, irrep_elements, dopant_configs, 
                 antisite_configs, interstitial_coords, include, exclude, 
                 symbreak, displace, symprec, oxidation_states, 
                 electron_negativity)

    def test_as_dict(self):
        d = self._a.as_dict()
#        print(d)

    def test_to_json(self):
        self.assertTrue(self._a.to_json())

class DefectInTest(unittest.TestCase):

    def setUp(self):
        structure = Structure.from_file("POSCAR-MgO64atoms")
        dopants = ["Al"]
        interstitial_coords = "0.1 0.1 0.1"
        is_antisite = True
        ElNeg_diff = 1.0
        include = ["Va_O1_-1", "Va_O1_-2"]                                      
        exclude = ["Va_O1_1", "Va_O1_2"]                                        
        symbreak = True                                                         
        displace = 0.15                                                         
        symprec = 0.001         

        self._b = DefectIn(structure, dopants, interstitial_coords, 
                           is_antisite, ElNeg_diff, include="", exclude="", 
                           symbreak=False, displace=0.2, symprec=0.01)

    def test_to(self):
        self._b.to()

    def test_from_str_file(self):
        self._c = DefectIn.from_str_file("POSCAR-MgO64atoms")

    

#    def test_dict(self):
#        print(self._a.as_dict())
#        print("------------------")
#        d = DefectSetting.from_dict(self._a.as_dict())
#        print(d.as_dict())
#        for k in d.as_dict():
            

#        print(self._a.as_dict())
#        self.assertTrue(self._a == d)

#    def test_from_defect_in(self):
#        print(self._a.__dict__)
#        f = DefectSetting.from_defect_in("POSCAR-MgO64atoms", "defect.in.test_MgO")
#        print(f.__dict__)
#        self.assertTrue(self._a == f)

