import unittest
import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.core.periodic_table import Element
from pydefect.taka_my_classes.perfect import Unitcell, Supercell

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 19, 2018"

TEST_DIRECTORY = "./test_files/SbLi2Na_dft_result/"


class UnitcellTest(unittest.TestCase):

    def setUp(self):
        self._unitcell = Unitcell(TEST_DIRECTORY)

    def test_structure(self):
        # expected values are from POSCAR file
        s = self._unitcell.structure
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
        e = self._unitcell.energy
        expected = -12.23262413
        self.assertAlmostEqual(e, expected)

    def test_dielectric_tensor(self):
        dielectric_electron = [[9.042893, 0.000000, 0.000000],
                               [-0.000000, 9.042893, 0.000000],
                               [0.000000, 0.000000, 9.042893]]

        dielectric_ionic = [[7.561335, -0.000000, 0.000000],
                            [-0.000000, 7.561323, 0.000000],
                            [0.000000, 0.000000, 7.561328]]
        expected = np.array(dielectric_electron) + np.array(dielectric_ionic)
        d = self._unitcell.dielectric_tensor
        np.testing.assert_allclose(d, expected)

    def test_eigen_value(self):
        ev = self._unitcell.eigen_value
        # from EIGENVAL file
        expected_first_line =\
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

    def test_optimize_ewald(self):
        ewald = self._unitcell.get_ewald_param()
        # Following expected value is printed by written code.
        # Numbers of real lattices and reciprocal lattices
        # are 11051 and 10692, respectively.
        expected_ewald = 0.4140771427121976
        self.assertAlmostEqual(ewald, expected_ewald)

    def test_generate_neighbor_lattices(self):
        expected = [[0.0, -3.39899999999999, -3.39899999999999],
                    [3.39899999999999, 0.0, -3.39899999999999],
                    [3.39899999999999, -3.39899999999999, 0.0],
                    [-3.39899999999999, 0.0, -3.39899999999999],
                    [0.0, 3.39899999999999, -3.39899999999999],
                    [-3.39899999999999, -3.39899999999999, 0.0],
                    [0.0, 0.0, 0.0],
                    [3.39899999999999, 3.39899999999999, 0.0],
                    [0.0, -3.39899999999999, 3.39899999999999],
                    [3.39899999999999, 0.0, 3.39899999999999],
                    [-3.39899999999999, 3.39899999999999, 0.0],
                    [-3.39899999999999, 0.0, 3.39899999999999],
                    [0.0, 3.39899999999999, 3.39899999999999]]
        count = 0
        for v in self._unitcell.generate_neighbor_lattices(
                3.4 * np.sqrt(2), include_self = True):
            count += 1
            is_v_included = any([abs(np.linalg.norm(e - v)) < 1e-3
                                 for e in expected])
            self.assertTrue(is_v_included)
            # TODO: Better to implement test of reciprocal lattice version
        self.assertEqual(count, len(expected))


class SupercellTest(unittest.TestCase):

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
        ev = self._supercell.eigen_value
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


if __name__ == "__main__":
    unittest.main()
