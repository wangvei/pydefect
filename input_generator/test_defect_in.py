import unittest
from defect_in import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

FILENAME_TO_JSON_FILE = "examples/defect_setting_test.json"

class DefectSettingTest(unittest.TestCase):

    def setUp(self):
        structure = Structure.from_file("POSCAR-MgO64atoms")
        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              repr_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             repr_coords=[0.25, 0.25, 0.25])
        irrep_elements = [Mg1, O1] 
        dopant_configs = "Al_Mg1"
        antisite_configs = ["Mg_O1", "O_Mg1"]
        interstitial_coords = [[0.1, 0.1, 0.1]]
        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        displace = 0.15
        cutoff = 3.0
        symprec = 0.001
        oxidation_states = {"Mg": 2, "O": -2}
        electronegativity = {"Mg": 1.31, "O": 3.44}

        self._a = DefectSetting(structure, irrep_elements, dopant_configs, 
                antisite_configs, interstitial_coords, included, excluded, 
                displace, cutoff, symprec, oxidation_states, electronegativity)

    def test_dict(self):
        self._dict = self._a.as_dict()
        self.object_from_dict = DefectSetting.from_dict(self._dict)
        self.assertTrue(vars(self.object_from_dict)) == vars(self._a)

    def test_to_json_file(self):
        self._a.to_json_file(FILENAME_TO_JSON_FILE)
        jf = DefectSetting.json_load(FILENAME_TO_JSON_FILE_VAC)
        self.assertTrue(vars(self._a) == vars(jf))


     def setUp(self):
         structure = Structure.from_file("POSCAR-MgO64atoms-strange_sequence")
         dopants = ["Al"]
         interstitial_coords = "0.1 0.1 0.1"
         is_antisite = True
         EN_diff = 100.0
         included = ["Va_O1_-1", "Va_O1_-2"]
         excluded = ["Va_O1_1", "Va_O1_2"]
         symbreak = True
         displace = 0.15
         symprec = 0.001

         self._b = DefectInMaker(structure, dopants, interstitial_coords,
                         is_antisite, EN_diff, included="", excluded="",
                         displace=0.2, cutoff=3.0, symprec=0.01)

     def a(self):
         print(self._b)

    def test_to(self):
        self._b.to()

    def test_from_str_file(self):
        self._c = DefectInMaker.from_str_file("POSCAR-MgO64atoms")



    def test_dict(self):
        print(self._a.as_dict())
        print("------------------")
        d = DefectSetting.from_dict(self._a.as_dict())
        print(d.as_dict())
        for k in d.as_dict():


        print(self._a.as_dict())
        self.assertTrue(self._a == d)

    def test_from_defect_in(self):
        print(self._a.__dict__)
        f = DefectSetting.from_defect_in("POSCAR-MgO64atoms", "defect.in.test_MgO")
        print(f.__dict__)
        self.assertTrue(self._a == f)

