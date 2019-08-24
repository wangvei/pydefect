# -*- coding: utf-8 -*-
import os
import unittest

from pydefect.util.vasp_util import element_diff_from_structures, \
    calc_participation_ratio, calc_orbital_character, calc_orbital_difference
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Procar
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ElementDiffFromPoscarFilesTest(PydefectTest):
    def test(self):
        structure = self.get_structure_by_name("MgO64atoms-Va_Mg+Va_O-unrelax")
        ref_structure = self.get_structure_by_name("MgO64atoms")

        actual = element_diff_from_structures(structure, ref_structure)
        expected = {"Mg": -1, "O": -1}
        self.assertEqual(expected, actual)


class CalcParticipationRatioTest(PydefectTest):
    def test(self):
        procar = Procar(self.DEFECTS_MGO_DIR / "Va_O1_2" / "PROCAR")

        actual = calc_participation_ratio(procar=procar,
                                          spin=Spin.up,
                                          band_index=0,
                                          atom_indices=[32, 33])
        #    33  0.051  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.051
        #    34  0.051  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.051
        # tot    0.732  0.001  0.001  0.001  0.000  0.000  0.000  0.000  0.000  0.736
        expected = 0.051 * 2 / 0.736
        self.assertAlmostEqual(expected, actual, 2)


class CalcOrbitalCharacterTest(PydefectTest):
    def test(self):
        procar = Procar(self.DEFECTS_MGO_DIR / "Va_O1_2" / "PROCAR")
        structure = self.get_structure_by_name("MgO64atoms-Va_O1_2")
        actual = calc_orbital_character(procar=procar,
                                        structure=structure,
                                        spin=Spin.up,
                                        band_index=0,
                                        kpoint_index=0)

        mg_sum = (0.002 + 0.001 + 0.001 + 0.001 + 0.002 + 0.001 + 0.001 + 0.001
                  + 0.001 + 0.001 + 0.001 + 0.001 + 0.001 + 0.001 + 0.001
                  + 0.001 + 0.002 + 0.002 + 0.001 + 0.001 + 0.001 + 0.001
                  + 0.001 + 0.001 + 0.002 + 0.001 + 0.002 + 0.001 + 0.001
                  + 0.001 + 0.001 + 0.001)
        o_sum = (0.051 + 0.051 + 0.017 + 0.051 + 0.017 + 0.017 + 0.008 + 0.028
                 + 0.024 + 0.024 + 0.028 + 0.014 + 0.014 + 0.014 + 0.014 + 0.024
                 + 0.028 + 0.014 + 0.014 + 0.028 + 0.024 + 0.014 + 0.014 + 0.024
                 + 0.014 + 0.028 + 0.014 + 0.028 + 0.014 + 0.024 + 0.014)

        expected = {'Mg': {'s': mg_sum, 'p': 0.0, 'd': 0.0, 'f': 0.0},
                    'O': {'s': o_sum, 'p': 0.0, 'd': 0.0, 'f': 0.0}}
        self.assertAlmostEqual(expected["Mg"]["s"], actual["Mg"]["s"], 4)
        self.assertAlmostEqual(expected["O"]["s"], actual["O"]["s"], 4)


class CalcOrbitalDifferenceTest(PydefectTest):
    def test(self):
        test_dir3 = os.path.join(test_dir, "MgO", "defects", "perfect")
        procar_perfect = Procar(os.path.join(test_dir3, "PROCAR"))
        structure_perfect = \
            Structure.from_file(os.path.join(test_dir3, "CONTCAR"))

        # MgO vbm
        orbital_1 = calc_orbital_character(procar=procar_perfect,
                                           structure=structure_perfect,
                                           spin=Spin.up,
                                           band_index=31,
                                           kpoint_index=0)

        test_dir4 = os.path.join(test_dir, "MgO", "defects", "Va_O1_2")
        procar_vo2 = Procar(os.path.join(test_dir4, "PROCAR"))
        structure_vo2 = \
            Structure.from_file(os.path.join(test_dir4, "CONTCAR"))

        # MgO vbm
        orbital_2 = calc_orbital_character(procar=procar_vo2,
                                           structure=structure_vo2,
                                           spin=Spin.up,
                                           band_index=27,
                                           kpoint_index=0)

        print(calc_orbital_difference(orbital_1, orbital_2))

