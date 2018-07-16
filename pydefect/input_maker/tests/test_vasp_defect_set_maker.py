# -*- coding: utf-8 -*-

import os
import shutil
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.vasp_defect_set_maker import VaspDefectInputSetMaker
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")

class VaspDefectInputSetMakerTest(unittest.TestCase):

    def setUp(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))

        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              representative_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             representative_coords=[0.25, 0.25, 0.25])
        irreducible_elements = [Mg1, O1]
        dopant_configs = [["Al", "Mg"]]
        antisite_configs = [["Mg", "O"], ["O", "Mg"]]
        interstitial_coords = [[0.1, 0.1, 0.1]]
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

    def test_mgo(self):
        test_mgo_dir = os.path.join(test_dir, "MgO")
        if os.path.exists(test_mgo_dir):
            shutil.rmtree(test_mgo_dir)
        os.mkdir(test_mgo_dir)
        os.chdir(test_mgo_dir)
        shutil.copyfile("../INCAR-MgO64atoms", "INCAR")
        shutil.copyfile("../KPOINTS-MgO64atoms", "KPOINTS")
        VaspDefectInputSetMaker(defect_initial_setting=self._mgo)
        VaspDefectInputSetMaker(defect_initial_setting=self._mgo,
                                particular_defects=["Sc_Mg1_0"])

    def test_vo_mgo(self):
        test_vo_mgo_dir = os.path.join(test_dir, "MgO_Va_O")
        if os.path.exists(test_vo_mgo_dir):
            shutil.rmtree(test_vo_mgo_dir)
        os.mkdir(test_vo_mgo_dir)
        os.chdir(test_vo_mgo_dir)
        shutil.copyfile("../INCAR-MgO64atoms", "INCAR")
        shutil.copyfile("../KPOINTS-MgO64atoms", "KPOINTS")
        # Note that the type of filtering_words is list.
        VaspDefectInputSetMaker(defect_initial_setting=self._mgo,
                                keywords=["perfect", "Va_O"])


if __name__ == "__main__":
    unittest.main()

