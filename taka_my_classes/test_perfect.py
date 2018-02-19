import unittest
import pydefect.taka_my_classes.perfect.py

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 19, 2018"

#TODO: make test directory(including poscar, outcar, vasprun.xml)
#currently only print type and values.
TEST_DIRECTORY = "TODO!"

class PerfectTest(unittest.TestCase):

    def setUp(self):
        self._perfect = Perfect(TEST_DIRECTORY)

    def test_structure(self):
        s = self._perfect.structure
        print(type(s))
        print(s)

    def test_energy(self):
        e = self._perfect.energy
        print(type(e))
        print(e)

    def test_dielectric_tensor(self):
        d = self._perfect.dielectric_tensor
        print(type(d))
        print(d)

    def test_electrostatic_potential(self):
        ep = self._perfect.electrostatic_potential
        print(type(ep))
        print(ep)

    def test_eigen_value(self):
        ev = self._perfect.eigen_value
        print(type(ev))
        print(ev)

