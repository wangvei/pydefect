# -*- coding: utf-8 -*-
import os
import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.supercell_maker import Supercell, Supercells

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class SupercellTest(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file("PPOSCAR-MgO")
#       from_file(os.path.join(test_dir, "PPOSCAR-MgO"))

    def test_init1(self):
        multi = [2, 1, 1]
        s1 = Supercell(structure=self.structure, trans_mat=multi, comment='')
        s1.structure.to(filename="PPOSCAR-MgO-2x1x1")

    def test_init2(self):
        multi = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        s2 = Supercell(structure=self.structure, trans_mat=multi, comment='')
        s2.structure.to(filename="PPOSCAR-MgO-conv")

    def test_recommended_supercell(self):
        s3 = Supercells(self.structure,
                        is_conventional=False,
                        max_num_atoms=500,
                        min_num_atoms=50,
                        isotropy_criterion=1.1)
        s = s3.create_sorted_supercells_by_num_atoms()[-1]

        s.structure.to(filename="POSCAR-MgO-recommended_supercell")
        print(s.comment)

