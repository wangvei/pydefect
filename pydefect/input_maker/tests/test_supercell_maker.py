# -*- coding: utf-8 -*-
import os
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.supercell_maker import Supercell, Supercells

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class SupercellTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("PPOSCAR-MgO")
#       from_file(os.path.join(test_dir, "PPOSCAR-MgO"))

    def test_init1(self):
        multi = [2, 1, 1]
        s1 = Supercell(structure=self.structure, multi=multi, comment='')
        s1.supercell_structure.to(filename="PPOSCAR-MgO-2x1x1")

    def test_init2(self):
        multi = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        s2 = Supercell(structure=self.structure, multi=multi, comment='')
        s2.supercell_structure.to(filename="PPOSCAR-MgO-conv")

    def test_recommended_supercell(self):
        s3 = Supercells(self.structure,
                        primitive=False,
                        max_num_atoms=200,
                        min_num_atoms=50,
                        isotropy_criterion=1.1)
        s = s3.sorted_by_isotropy()[0]

        s.supercell_structure.to(filename="POSCAR-MgO-recommended_supercell")
        print(s.comment)

