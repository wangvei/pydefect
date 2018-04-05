import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.vasp_input_maker import *
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

FILENAME_POSCAR = "examples/POSCAR-MgO64atoms"
FILENAME_POSCAR_VA1 = "examples/POSCAR-MgO64atoms-Va_Mg1"
FILENAME_POSCAR_IN1 = "examples/POSCAR-MgO64atoms-O_i1"
FILENAME_POSCAR_IN3 = "examples/POSCAR-MgO64atoms-Al_i1"
FILENAME_POSCAR_AS1 = "examples/POSCAR-MgO64atoms-O_Mg"
FILENAME_POSCAR_SS1 = "examples/POSCAR-MgO64atoms-N_O1"
#FILENAME_JSON_VA1 = "examples/MgO64atoms-Va_Mg1.json"


class VaspDefectInputSetMakerTest(unittest.TestCase):

    def setUp(self):
        structure = Structure.from_file(FILENAME_POSCAR)
        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              representative_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             representative_coords=[0.25, 0.25, 0.25])
        irreducible_elements = [Mg1, O1]
        dopant_configs = [["Al", "Mg"]]
        antisite_configs = []
        interstitial_coords = [[0.1, 0.1, 0.1]]
#        self.interstitial_coords = interstitial_coords
        included = ["Va_O1_-1", "Va_O1_-2"]
        excluded = ["Va_O1_1", "Va_O1_2"]
        distance = 0.15
        cutoff = 2.0
        symprec = 0.001
        oxidation_states = {"Mg": 2, "O": -2, "Al": 3, "N": -3}
        electronegativity = {"Mg": 1.31, "O": 3.44, "Al": 1.61, "N": 3.04}

        self._mgo = DefectInitialSetting(
            structure, irreducible_elements, dopant_configs, antisite_configs,
            interstitial_coords, included, excluded, distance, cutoff,
            symprec, oxidation_states, electronegativity)

    def test(self):
        os.chdir("examples")
        VaspDefectInputSetMaker(defect_initial_setting=self._mgo)
        VaspDefectInputSetMaker(defect_initial_setting=self._mgo,
                                particular_defects=["Sc_Mg1_0"])


if __name__ == "__main__":
    unittest.main()

