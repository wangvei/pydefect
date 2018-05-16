# -*- coding: utf-8 -*-
import os
import unittest

from pydefect.vasp_util.util import element_diff_from_poscar_files

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class ElementDiffFromPoscarFilesTest(unittest.TestCase):
    def test(self):
        poscar1 = os.path.join(test_dir, "POSCAR-MgO8atoms")
        poscar2 = os.path.join(test_dir, "POSCAR-MgO8atoms-Va_O1+N_O")

        actual = element_diff_from_poscar_files(poscar1, poscar2)
        expected = {"O": 2, "N": -1}
        self.assertEqual(actual, expected)