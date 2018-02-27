#!/usr/bin/env python

import unittest
from defect import *
from pymatgen.core.structure import Structure

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

FILENAME_VAC_POSCAR = "examples/POSCAR-MgO64atoms-Va_Mg1"
FILENAME_INT_POSCAR = "examples/POSCAR-MgO64atoms-O_i1"
FILENAME_TO_JSON_FILE_VAC = "examples/Va_Mg1_-2.json"
FILENAME_TO_JSON_FILE_INT = "examples/Mg_i1_1.json"

TEST_DIRECTORY = "../analysis/test_files/SbLi2Na_dft_result/"


class PerfectTest(unittest.TestCase):

    def setUp(self):
        self._supercell = Supercell(TEST_DIRECTORY)

    def test_structure(self):
        # expected values are from POSCAR file
        s = self._supercell.structure
        lattice_val = 3.3989999999999898
        expected_lattice = [[0, lattice_val, lattice_val],
                            [lattice_val, 0, lattice_val],
                            [lattice_val, lattice_val, 0]]
        np.testing.assert_allclose(np.array(s.lattice.matrix), expected_lattice)
        self.assertEqual(s.sites[0].specie, Element("Sb"))
        self.assertEqual(s.sites[1].specie, Element("Li"))
        self.assertEqual(s.sites[2].specie, Element("Li"))
        self.assertEqual(s.sites[3].specie, Element("Na"))
        np.allclose(s.sites[0].coords, np.array([0.00, 0.00, 0.00]))
        np.allclose(s.sites[1].coords, np.array([0.00, 0.75, 0.75]))
        np.allclose(s.sites[2].coords, np.array([0.00, 0.25, 0.25]))
        np.allclose(s.sites[3].coords, np.array([0.50, 0.50, 0.50]))

    def test_energy(self):
        e = self._supercell.energy
        expected = -12.23262413
        self.assertAlmostEqual(e, expected)

    def test_electrostatic_potential(self):
        ep = self._supercell.electrostatic_potential
        # from OUTCAR file
        expected = [-87.6376,
                    -26.7608,
                    -26.7608,
                    -39.9980]
        np.testing.assert_allclose(np.array(ep), np.array(expected))

    def test_eigen_value(self):
        ev = self._supercell.eigenvalue
        # from EIGENVAL file
        expected_first_line = \
            [[-5.461771, 1.000000],
             [3.574613, 1.000000],
             [3.574613, 1.000000],
             [3.574614, 1.000000],
             [5.252315, 0.000000],
             [6.222399, 0.000000],
             [6.222399, 0.000000],
             [6.222402, 0.000000],
             [9.639593, 0.000000],
             [9.639593, 0.000000],
             [9.639599, 0.000000],
             [10.042995, 0.000000],
             [10.042995, 0.000000],
             [11.587756, 0.000000],
             [16.725048, 0.000000],
             [19.028426, 0.000000]]
        actual_first_line = ev[Spin.up][0]
        np.testing.assert_allclose(np.array(actual_first_line),
                                   np.array(expected_first_line), rtol = 1e-4)


class DefectTest(unittest.TestCase):

    def setUp(self):
        _vac_structure = Structure.from_file(FILENAME_VAC_POSCAR)
        _int_structure = Structure.from_file(FILENAME_INT_POSCAR)

        self._vac = \
            DefectInput(initial_structure=_vac_structure, removed_atom_index=0,
                        inserted_atom_index=None, defect_coords=[0, 0, 0],
                        in_name="Va", out_name="Mg1", charge=-2)
        self._int = \
            DefectInput(initial_structure=_int_structure, removed_atom_index=None,
                        inserted_atom_index=32, defect_coords=[0.1, 0.1, 0.1],
                        in_name="O", out_name="i1", charge=-2)

    def test_dict(self):
        d_vac = DefectInput.from_dict(self._vac.as_dict())
        d_int = DefectInput.from_dict(self._int.as_dict())
        self.assertTrue(vars(self._vac) == vars(d_vac))
        self.assertTrue(vars(self._int) == vars(d_int))

    def test_to_json_file(self):
        self._vac.to_json_file(FILENAME_TO_JSON_FILE_VAC)
        jf_vac = DefectInput.json_load(FILENAME_TO_JSON_FILE_VAC)
        self.assertTrue(vars(self._vac) == vars(jf_vac))

        self._int.to_json_file(FILENAME_TO_JSON_FILE_INT)
        jf_int = DefectInput.json_load(FILENAME_TO_JSON_FILE_INT)
        self.assertTrue(vars(self._int) == vars(jf_int))

    def test_atom_mapping_to_perfect(self):
        atom_mapping_vac_test = [i for i in range(1, 64)]
#        print(atom_mapping_vac_test)
#        print(self._vac.atom_mapping_to_perfect())
        self.assertEqual(
            self._vac.atom_mapping_to_perfect(), atom_mapping_vac_test)

        atom_mapping_int_test = \
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
             19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, None, 32, 33,
             34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
             51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

#        print(atom_mapping_int_test)
#        print(self._int.atom_mapping_to_perfect())
        self.assertEqual(
            self._int.atom_mapping_to_perfect(), atom_mapping_int_test)

    def test_anchor_atom_index(self):
        correct_farthest_atom_site = 6
        self.assertEqual(self._vac.anchor_atom_index()[0],
                         correct_farthest_atom_site)


class IrreducibleSiteTest(unittest.TestCase):
    
    def setUp(self):
        self._mg = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                                   first_index=1, last_index=32,
                                   repr_coords=[0, 0, 0])

    def test_dict(self):
        d = IrreducibleSite.from_dict(self._mg.as_dict())
        self.assertTrue(self._mg == d)

    def test_natoms(self):
        self.assertEqual(self._mg.natoms, 32)

if __name__ == "__main__":
    unittest.main()