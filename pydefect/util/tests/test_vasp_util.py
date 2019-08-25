# -*- coding: utf-8 -*-

from pydefect.util.testing import PydefectTest
from pydefect.util.vasp_util import element_diff_from_structures, \
    calc_participation_ratio, calc_orbital_character, calc_orbital_difference
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Procar

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ElementDiffFromPoscarFilesTest(PydefectTest):
    def test_element_diff_from_poscar_files(self):
        structure = self.get_structure_by_name("MgO64atoms-Va_Mg+Va_O-unrelax")
        ref_structure = self.get_structure_by_name("MgO64atoms")

        actual = element_diff_from_structures(structure, ref_structure)
        expected = {"Mg": -1, "O": -1}
        self.assertEqual(expected, actual)


class CalcParticipationRatioTest(PydefectTest):
    def test_calc_participation_ratio(self):
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
    def test_calc_orbital_character(self):
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
    def setUp(self) -> None:
        procar_perfect = Procar(self.DEFECTS_MGO_DIR / "perfect" / "PROCAR")
        structure_perfect = self.get_structure_by_name("MgO64atoms")

        procar_va_o_2 = Procar(self.DEFECTS_MGO_DIR / "Va_O1_2" / "PROCAR")
        structure_va_o_2 = self.get_structure_by_name("MgO64atoms-Va_O1_2")

        # MgO vbm
        self.orbital_1 = calc_orbital_character(procar=procar_perfect,
                                                structure=structure_perfect,
                                                spin=Spin.up,
                                                band_index=127,
                                                kpoint_index=0)

        # MgO vbm
        self.orbital_2 = calc_orbital_character(procar=procar_va_o_2,
                                                structure=structure_va_o_2,
                                                spin=Spin.up,
                                                band_index=123,
                                                kpoint_index=0)

    def test_calc_orbital_difference(self):
        # {'Mg': {'s': 0.0, 'p': 0.0, 'd': 0.0, 'f': 0.0},
        #   'O': {'s': 0.0, 'p': 0.736, 'd': 0.0, 'f': 0.0}}
        # {'Mg': {'s': 0.0, 'p': 0.004, 'd': 0.0, 'f': 0.0},
        #   'O': {'s': 0.0, 'p': 0.73, 'd': 0.0, 'f': 0.0}}

        actual = calc_orbital_difference(self.orbital_1, self.orbital_2)
        self.assertAlmostEqual(0.01, actual, 3)

