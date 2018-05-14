#!/usr/bin/env python
# -*- coding: utf-8 -*-
import filecmp
import unittest

from pydefect.input_maker.vasp_input_maker import ModIncar, make_band_kpoints, \
    make_kpoints, make_incar

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "April 4, 2018"


class ModIncarTest(unittest.TestCase):

    def test_ModIncar(self):
        i = ModIncar.from_file("INCAR-ModIncar_before")
        i.pretty_write_file("INCAR-ModIncar_actual")
        self.assertTrue(
            filecmp.cmp("INCAR-ModIncar_actual", "INCAR-ModIncar_expected"))


class MakeBandKpointsTest(unittest.TestCase):

    def test_make_band_kpoints_mgo(self):
        make_band_kpoints(ibzkpt="IBZKPT", poscar="PPOSCAR-MgO",
                          pposcar="PPOSCAR-MgO", num_split_kpoints=2)
        self.assertTrue(
            filecmp.cmp("KPOINTS-0", "KPOINTS-make_band_kpoints_expected-0"))
        self.assertTrue(
            filecmp.cmp("KPOINTS-1", "KPOINTS-make_band_kpoints_expected-1"))

        make_band_kpoints(ibzkpt="IBZKPT", poscar="PPOSCAR-YMnO3-bc_exchanged",
                          pposcar="PPOSCAR-YMnO3", num_split_kpoints=1)


class MakeKpointsTest(unittest.TestCase):

    def test_make_kpoints(self):
        # make_kpoints(task="structure_opt", poscar="PPOSCAR-MgO")
        # self.assertTrue(filecmp.cmp("KPOINTS",
        #                             "KPOINTS-so_expected_MgO"))
        # make_kpoints(task="structure_opt", poscar="PPOSCAR-YMnO3")
        # self.assertTrue(filecmp.cmp("KPOINTS",
        #                             "KPOINTS-so_expected_YMnO3"))
        # make_kpoints(task="dos", poscar="PPOSCAR-YMnO3")
        # self.assertTrue(filecmp.cmp("KPOINTS",
        #                             "KPOINTS-dos_expected_YMnO3"))
        #   make_kpoints(task="competing_phase", poscar="PPOSCAR-Ru",
        #                is_metal=True)
        #   self.assertTrue(filecmp.cmp("KPOINTS",
        #                               "KPOINTS-cp_expected_Ru"))
        make_kpoints(task="competing_phase_molecule", poscar="POSCAR-F2")

    #    def test_make_band_kpoints(self):

    #       make_kpoints(task="band", ibzkpt="IBZKPT", poscar="BPOSCAR-MgO",
    #                    pposcar="PPOSCAR-MgO", num_split_kpoints=1)
    #       self.assertTrue(
    #           filecmp.cmp("KPOINTS", "KPOINTS-make_band_kpoints_expected"))

    #        make_kpoints(task="defect")

    #        make_incar(task="defect", functional="hse06")


class MakeIncarTest(unittest.TestCase):

    def test_make_incar(self):
        make_incar(task="structure_opt", functional="hse06",
                   defect_in="defect.in", hfscreen=0.2, aexx=0.3,
                   is_magnetization=True,
                   my_incar_setting="my_INCAR_setting.yaml")

