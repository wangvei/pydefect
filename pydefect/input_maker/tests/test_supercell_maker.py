# -*- coding: utf-8 -*-
import os

import numpy as np
from pydefect.input_maker.supercell_maker import (
    Supercell, Supercells, calc_isotropy, sanitize_matrix)
from pydefect.util.testing import PydefectTest
from pymatgen.core import IStructure, Lattice


class CalcIsotropyTest(PydefectTest):
    def setUp(self):
        coords = [[0, 0, 0]]
        lattice = Lattice([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        self.struct = IStructure(lattice, ["H"], coords)

    def test(self):
        trans_mat = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
        actual = calc_isotropy(self.struct, trans_mat)
        # (1 + 1) / 3 / 3
        expected = (0.2222, 90)
        self.assertEqual(expected, actual)


class SanitizeMatrixTest(PydefectTest):
    def setUp(self) -> None:
        self.matrix1 = [2]
        self.matrix2 = [2, 3, 4]
        self.matrix3 = [[2, 3, 4], [3, 6, 7], [4, 7, 5]]
        self.matrix4 = [2, 3, 4, 3, 6, 7, 4, 7, 5]

    def test_len1(self):
        expected = np.array([[2, 0, 0],
                             [0, 2, 0],
                             [0, 0, 2]])
        actual = sanitize_matrix(self.matrix1)
        self.assertArrayAlmostEqual(expected, actual, 7)

    def test_len2(self):
        expected = np.array([[2, 0, 0],
                             [0, 3, 0],
                             [0, 0, 4]])
        actual = sanitize_matrix(self.matrix2)
        self.assertArrayAlmostEqual(expected, actual, 7)

    def test_len3(self):
        expected = np.array([[2, 3, 4],
                             [3, 6, 7],
                             [4, 7, 5]])
        actual = sanitize_matrix(self.matrix3)
        self.assertArrayAlmostEqual(expected, actual, 7)

    def test_len4(self):
        expected = np.array([[2, 3, 4],
                             [3, 6, 7],
                             [4, 7, 5]])
        actual = sanitize_matrix(self.matrix4)
        self.assertArrayAlmostEqual(expected, actual, 7)


class SupercellTest(PydefectTest):
    def setUp(self):
        self.mgo_struct = self.get_structure_by_name("MgO")
        self.mgo64_struct = self.get_structure_by_name("MgO64atoms")
        self.kzn4p3_struct = self.get_structure_by_name("KZn4P3")

    def test_init1(self):
        multi = [2, 1, 1]
        actual = Supercell(structure=self.mgo_struct,
                           trans_mat=multi,
                           multiplicity=2).structure
        expected = IStructure.from_str(
            """Mg2 O2
             1.0
             0.000000 4.246894 4.246894
             2.123447 0.000000 2.123447
             2.123447 2.123447 0.000000
             Mg O
             2 2
             direct
             0.000000 0.000000 0.000000 Mg
             0.500000 0.000000 0.000000 Mg
             0.250000 0.500000 0.500000 O
             0.750000 0.500000 0.500000 O""", fmt="poscar")
        self.assertEqual(expected, actual)

    def test_init2(self):
        multi = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        actual = Supercell(structure=self.mgo_struct,
                           trans_mat=multi,
                           multiplicity=4).structure
        expected = IStructure.from_str(
            """Mg4 O4
            1.0
            4.246894 0.000000 0.000000
            0.000000 4.246894 0.000000
            0.000000 0.000000 4.246894
            Mg O
            4 4
            direct
            0.000000 0.000000 0.000000 Mg
            0.500000 0.500000 0.000000 Mg
            0.500000 0.000000 0.500000 Mg
            0.000000 0.500000 0.500000 Mg
            0.500000 0.500000 0.500000 O
            1.000000 1.000000 0.500000 O
            1.000000 0.500000 1.000000 O
            0.500000 1.000000 1.000000 O""", fmt="poscar")
        self.assertEqual(expected, actual)

    def test_supercell_input(self):
        multi = [2, 1, 1]
        supercell = Supercell(structure=self.mgo64_struct,
                              trans_mat=multi,
                              check_unitcell=True)
        supercell.to("POSCAR")
        os.remove("POSCAR")
        os.remove("UPOSCAR")
        actual = supercell.structure
        expected = IStructure.from_str(
        """trans_mat: [3.9999999999999996, 0.0, 0.0, 0.0, 1.9999999999999998, 0.0, 0.0, 0.0, 1.9999999999999998], multi: 16, isotropy: (0.5333, 59.99999999999999)
        1.0
        0.000000 6.314592 6.314592
        2.104864 0.000000 2.104864
        2.104864 2.104864 0.000000
        Mg O
        3 3
        direct
        0.000000 0.000000 0.000000 Mg
        0.333333 0.000000 0.000000 Mg
        0.666667 0.000000 -0.000000 Mg
        0.166667 0.500000 0.500000 O
        0.500000 0.500000 0.500000 O
        0.833333 0.500000 0.500000 O""", fmt="poscar")
        self.assertEqual(expected, actual)


class SupercellsTest(PydefectTest):
    def setUp(self):
        self.mgo_struct = self.get_structure_by_name("MgO")
        self.kzn4p3_struct = self.get_structure_by_name("KZn4P3")
        self.tetra_struct = self.get_structure_by_name("Zn3P2")

    def test_recommended_supercell_conv(self):
        supercells_mgo = Supercells(self.mgo_struct,
                                    conventional_base=True,
                                    max_num_atoms=500,
                                    min_num_atoms=1,
                                    criterion=1.1)
        supercell_mgo = supercells_mgo.sorted_supercells_by_num_atoms[0]
        actual = supercell_mgo.structure
        supercell_mgo.to("POSCAR-MgO")
        os.remove("POSCAR-MgO")
        expected = IStructure.from_str(
            """trans_mat: -1 1 1 1 -1 1 1 1 -1, multi: 4.0, isotropy: 0.0
            1.0
            4.246894 0.000000 0.000000
            0.000000 4.246894 0.000000
            0.000000 0.000000 4.246894
            Mg O
            4 4
            direct
            0.000000 0.000000 0.000000 Mg
            0.500000 0.500000 0.000000 Mg
            0.500000 0.000000 0.500000 Mg
            0.000000 0.500000 0.500000 Mg
            0.500000 0.500000 0.500000 O
            1.000000 1.000000 0.500000 O
            1.000000 0.500000 1.000000 O
            0.500000 1.000000 1.000000 O""", fmt="poscar")
        self.assertEqual(expected, actual)

    def test_recommended_supercell(self):
        supercells_mgo = Supercells(self.mgo_struct,
                                    conventional_base=False,
                                    max_num_atoms=500,
                                    min_num_atoms=16,
                                    criterion=1.1)
        supercell_mgo = supercells_mgo.sorted_supercells_by_num_atoms[0]
        actual = supercell_mgo.structure
        expected = IStructure.from_str(
            """Mg8 O8
            1.0
            0.000000 4.246894 4.246894
            4.246894 0.000000 4.246894
            4.246894 4.246894 0.000000
            Mg O
            8 8
            direct
            0.000000 0.000000 0.000000 Mg
            0.000000 0.000000 0.500000 Mg
            0.000000 0.500000 0.000000 Mg
            0.000000 0.500000 0.500000 Mg
            0.500000 0.000000 0.000000 Mg
            0.500000 0.000000 0.500000 Mg
            0.500000 0.500000 0.000000 Mg
            0.500000 0.500000 0.500000 Mg
            0.250000 0.250000 0.250000 O
            0.250000 0.250000 0.750000 O
            0.250000 0.750000 0.250000 O
            0.250000 0.750000 0.750000 O
            0.750000 0.250000 0.250000 O
            0.750000 0.250000 0.750000 O
            0.750000 0.750000 0.250000 O
            0.750000 0.750000 0.750000 O""", fmt="poscar")
        self.assertEqual(expected, actual)
        actual_comment = supercell_mgo.comment
        expected_comment = "trans_mat: [2, 0, 0, 0, 2, 0, 0, 0, 2], multi: 8, isotropy: (0.0, 59.99999999999999)"
        self.assertEqual(expected_comment, actual_comment)

    def test_recommended_supercell_rhombo(self):
        supercells_kzn4p3 = Supercells(self.kzn4p3_struct,
                                       conventional_base=False,
                                       max_num_atoms=500,
                                       criterion=1.1)
        supercell = supercells_kzn4p3.sorted_supercells_by_num_atoms[0]
        actual = supercell.structure
        actual.to(fmt="poscar", filename="POSCAR")
        os.remove("POSCAR")
        expected = IStructure.from_str(
            """K16 Zn64 P48
            1.0
            -7.933452 4.580381 11.216710
            0.000000 -9.160762 11.216710
            7.933452 4.580381 11.216710
            K Zn P
            16 64 48
            direct
            0.750000 1.250000 0.500000 K
            1.000000 1.000000 0.500000 K
            1.250000 0.750000 0.500000 K
            0.500000 1.250000 0.750000 K
            0.500000 0.500000 0.500000 K
            0.750000 1.000000 0.750000 K
            1.000000 0.750000 0.750000 K
            1.250000 1.250000 1.000000 K
            1.250000 0.500000 0.750000 K
            0.500000 1.000000 1.000000 K
            0.750000 0.750000 1.000000 K
            1.000000 1.250000 1.250000 K
            1.000000 0.500000 1.000000 K
            1.250000 1.000000 1.250000 K
            0.500000 0.750000 1.250000 K
            0.750000 0.500000 1.250000 K
            0.946131 1.446131 0.696131 Zn
            1.196131 1.196131 0.696131 Zn
            1.446131 0.946131 0.696131 Zn
            0.696131 1.446131 0.946131 Zn
            0.696131 0.696131 0.696131 Zn
            0.946131 1.196131 0.946131 Zn
            1.196131 0.946131 0.946131 Zn            
            1.446131 1.446131 1.196131 Zn
            1.446131 0.696131 0.946131 Zn
            0.696131 1.196131 1.196131 Zn
            0.946131 0.946131 1.196131 Zn
            1.196131 1.446131 1.446131 Zn
            1.196131 0.696131 1.196131 Zn
            1.446131 1.196131 1.446131 Zn
            0.696131 0.946131 1.446131 Zn
            0.946131 0.696131 1.446131 Zn
            0.553869 1.053869 0.303869 Zn
            0.803869 0.803869 0.303869 Zn
            1.053869 0.553869 0.303869 Zn
            0.303869 1.053869 0.553869 Zn
            0.303869 0.303869 0.303869 Zn
            0.553869 0.803869 0.553869 Zn
            0.803869 0.553869 0.553869 Zn
            1.053869 1.053869 0.803869 Zn
            1.053869 0.303869 0.553869 Zn
            0.303869 0.803869 0.803869 Zn
            0.553869 0.553869 0.803869 Zn
            0.803869 1.053869 1.053869 Zn
            0.803869 0.303869 0.803869 Zn
            1.053869 0.803869 1.053869 Zn
            0.303869 0.553869 1.053869 Zn
            0.553869 0.303869 1.053869 Zn
            1.166211 1.666211 0.916211 Zn
            1.416211 1.416211 0.916211 Zn
            1.666211 1.166211 0.916211 Zn
            0.916211 1.666211 1.166211 Zn
            0.916211 0.916211 0.916211 Zn
            1.166211 1.416211 1.166211 Zn
            1.416211 1.166211 1.166211 Zn
            1.666211 1.666211 1.416211 Zn
            1.666211 0.916211 1.166211 Zn
            0.916211 1.416211 1.416211 Zn
            1.166211 1.166211 1.416211 Zn
            1.416211 1.666211 1.666211 Zn
            1.416211 0.916211 1.416211 Zn
            1.666211 1.416211 1.666211 Zn
            0.916211 1.166211 1.666211 Zn
            1.166211 0.916211 1.666211 Zn
            0.333789 0.833789 0.083789 Zn
            0.583789 0.583789 0.083789 Zn
            0.833789 0.333789 0.083789 Zn
            0.083789 0.833789 0.333789 Zn
            0.083789 0.083789 0.083789 Zn
            0.333789 0.583789 0.333789 Zn
            0.583789 0.333789 0.333789 Zn
            0.833789 0.833789 0.583789 Zn
            0.833789 0.083789 0.333789 Zn
            0.083789 0.583789 0.583789 Zn
            0.333789 0.333789 0.583789 Zn
            0.583789 0.833789 0.833789 Zn
            0.583789 0.083789 0.583789 Zn
            0.833789 0.583789 0.833789 Zn
            0.083789 0.333789 0.833789 Zn
            0.333789 0.083789 0.833789 Zn
            1.014789 1.514789 0.764789 P
            1.264789 1.264789 0.764789 P
            1.514789 1.014789 0.764789 P
            0.764789 1.514789 1.014789 P
            0.764789 0.764789 0.764789 P
            1.014789 1.264789 1.014789 P
            1.264789 1.014789 1.014789 P
            1.514789 1.514789 1.264789 P
            1.514789 0.764789 1.014789 P
            0.764789 1.264789 1.264789 P
            1.014789 1.014789 1.264789 P
            1.264789 1.514789 1.514789 P
            1.264789 0.764789 1.264789 P
            1.514789 1.264789 1.514789 P
            0.764789 1.014789 1.514789 P
            1.014789 0.764789 1.514789 P
            0.485211 0.985211 0.235211 P
            0.735211 0.735211 0.235211 P
            0.985211 0.485211 0.235211 P
            0.235211 0.985211 0.485211 P
            0.235211 0.235211 0.235211 P
            0.485211 0.735211 0.485211 P
            0.735211 0.485211 0.485211 P
            0.985211 0.985211 0.735211 P
            0.985211 0.235211 0.485211 P
            0.235211 0.735211 0.735211 P
            0.485211 0.485211 0.735211 P
            0.735211 0.985211 0.985211 P
            0.735211 0.235211 0.735211 P
            0.985211 0.735211 0.985211 P
            0.235211 0.485211 0.985211 P
            0.485211 0.235211 0.985211 P
            0.250000 0.750000 0.000000 P
            0.500000 0.500000 -0.000000 P
            0.750000 0.250000 0.000000 P
            0.000000 0.750000 0.250000 P
            0.000000 0.000000 0.000000 P
            0.250000 0.500000 0.250000 P
            0.500000 0.250000 0.250000 P
            0.750000 0.750000 0.500000 P
            0.750000 -0.000000 0.250000 P
            -0.000000 0.500000 0.500000 P
            0.250000 0.250000 0.500000 P
            0.500000 0.750000 0.750000 P
            0.500000 -0.000000 0.500000 P
            0.750000 0.500000 0.750000 P
            0.000000 0.250000 0.750000 P
            0.250000 -0.000000 0.750000 P""", fmt="poscar")
        self.assertEqual(expected, actual)
        actual_comment = supercell.comment
        expected_comment = "trans_mat: [-1, 3, -1, -1, -1, 3, 3, -1, -1], multi: 16, isotropy: (0.0, 66.43340130273333)"
        self.assertEqual(expected_comment, actual_comment)

    def test_recommended_supercell_tetra(self):
        supercells_zn3p2 = Supercells(self.tetra_struct,
                                      conventional_base=True,
                                      criterion=1.1)
        supercell = supercells_zn3p2.sorted_supercells_by_isotropy[0]
        actual = supercell.structure
        supercell.to("POSCAR-Zn3P2")
        os.remove("POSCAR-Zn3P2")
        expected = IStructure.from_str(
            """trans_mat: [1, 1, 0, -1, 1, 0, 0, 0, 1], multi: 2, isotropy: (0.0007, 90.0)
            1.0
            8.120691 8.120691 0.000000
            -8.120691 8.120691 0.000000
            0.000000 0.000000 11.466590
            Zn P
            48 32
            direct
            0.608541 0.108541 0.606249 Zn
            1.108541 0.608541 0.606249 Zn
            0.391459 -0.391459 0.106249 Zn
            0.891459 0.108541 0.106249 Zn
            0.108541 -0.108541 0.106249 Zn
            0.608541 0.391459 0.106249 Zn
            0.391459 0.391459 0.893751 Zn
            0.891459 0.891459 0.893751 Zn
            0.108541 0.108541 0.893751 Zn
            0.608541 0.608541 0.893751 Zn
            0.391459 -0.108541 0.606249 Zn
            0.891459 0.391459 0.606249 Zn
            0.608541 -0.108541 0.393751 Zn
            1.108541 0.391459 0.393751 Zn
            0.391459 0.108541 0.393751 Zn
            0.891459 0.608541 0.393751 Zn
            0.641522 0.141522 0.885170 Zn
            1.141522 0.641521 0.885170 Zn
            0.358478 -0.358478 0.385170 Zn
            0.858478 0.141522 0.385170 Zn
            0.141522 -0.141522 0.385170 Zn
            0.641522 0.358478 0.385170 Zn
            0.358478 0.358478 0.614830 Zn
            0.858479 0.858479 0.614830 Zn
            0.141522 0.141522 0.614830 Zn
            0.641522 0.641522 0.614830 Zn
            0.358478 -0.141522 0.885170 Zn
            0.858479 0.358479 0.885170 Zn
            0.641522 -0.141522 0.114830 Zn
            1.141522 0.358478 0.114830 Zn
            0.358478 0.141522 0.114830 Zn
            0.858478 0.641522 0.114830 Zn
            0.623398 0.123397 0.148155 Zn
            1.123398 0.623398 0.148155 Zn
            0.376603 -0.376603 0.648155 Zn
            0.876602 0.123397 0.648155 Zn
            0.123398 -0.123398 0.648155 Zn
            0.623398 0.376602 0.648155 Zn
            0.376603 0.376603 0.351845 Zn
            0.876602 0.876602 0.351845 Zn
            0.123398 0.123398 0.351845 Zn
            0.623398 0.623398 0.351845 Zn
            0.376603 -0.123397 0.148155 Zn
            0.876602 0.376602 0.148155 Zn
            0.623398 -0.123397 0.851845 Zn
            1.123398 0.376603 0.851845 Zn
            0.376603 0.123397 0.851845 Zn
            0.876602 0.623398 0.851845 Zn
            0.500000 0.000000 0.746106 P
            1.000000 0.500000 0.746106 P
            0.000000 0.000000 0.246106 P
            0.500000 0.500000 0.246106 P
            0.000000 0.000000 0.753894 P
            0.500000 0.500000 0.753894 P
            0.500000 0.000000 0.253894 P
            1.000000 0.500000 0.253894 P
            0.250000 -0.250000 0.745801 P
            0.750000 0.250000 0.745801 P
            0.250000 -0.250000 0.245801 P
            0.750000 0.250000 0.245801 P
            0.250000 0.250000 0.754199 P
            0.750000 0.750000 0.754199 P
            0.250000 0.250000 0.254199 P
            0.750000 0.750000 0.254199 P
            0.754391 0.000000 0.500000 P
            1.254391 0.500000 0.500000 P
            0.500000 -0.245609 0.000000 P
            1.000000 0.254391 0.000000 P
            0.500000 0.245609 0.000000 P
            1.000000 0.745609 0.000000 P
            0.245609 0.000000 0.500000 P
            0.745609 0.500000 0.500000 P
            0.745609 0.000000 0.000000 P
            1.245609 0.500000 0.000000 P
            0.500000 -0.254391 0.500000 P
            1.000000 0.245609 0.500000 P
            0.500000 0.254391 0.500000 P
            1.000000 0.754391 0.500000 P
            0.254391 0.000000 0.000000 P
            0.754391 0.500000 0.000000 P""", fmt="poscar")
        self.assertEqual(expected, actual)
