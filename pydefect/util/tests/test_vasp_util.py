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
        ref_structure =  self.get_structure_by_name("MgO64atoms")

        actual = element_diff_from_structures(structure, ref_structure)
        expected = {"Mg": -1, "O": -1}
        self.assertEqual(expected, actual)


class CalcParticipationRatioTest(PydefectTest):
    def test(self):
        procar = Procar(os.path.join(test_dir2, "CO", "PROCAR"))

        actual = calc_participation_ratio(procar=procar, spin=Spin.up,
                                          band_index=6, kpoint_index=0,
                                          atom_indices=[1])
        expected = 0.39768339768339767
        self.assertEqual(actual, expected)

        actual = calc_participation_ratio(procar=procar, spin=Spin.down,
                                          band_index=7, kpoint_index=0,
                                          atom_indices=[0, 1])
        self.assertEqual(actual, 1.0)


class CalcOrbitalCharacterTest(PydefectTest):
    def test(self):
        procar = Procar(os.path.join(test_dir2, "CO", "PROCAR"))
        structure = Structure.from_file(os.path.join(test_dir2, "CO", "CONTCAR"))

        actual = calc_orbital_character(procar=procar, structure=structure,
                                        spin=Spin.up, band_index=6,
                                        kpoint_index=0)
        print(actual)


class CalcOrbitalSimilarityTest(PydefectTest):
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

