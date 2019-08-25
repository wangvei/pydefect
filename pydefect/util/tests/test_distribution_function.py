# -*- coding: utf-8 -*-

import numpy as np
from pydefect.util.distribution_function import fermi_dirac_distribution, \
    bose_einstein_distribution, maxwell_boltzmann_distribution
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class DistributionTest(PydefectTest):

    def setUp(self) -> None:
        self.energy = 20
        self.fermi_level = 3
        self.temperature = 10000

    def test_fermi_dirac_distribution(self):
        actual = fermi_dirac_distribution(self.energy, self.fermi_level,
                                          self.temperature)
        expected = np.reciprocal(np.exp((20 - 3) / (10000 / 11604.505)) + 1)

        self.assertAlmostEqual(expected, actual)

    def test_bose_einstein_distribution(self):
        actual = bose_einstein_distribution(self.energy, self.fermi_level,
                                            self.temperature)
        expected = np.reciprocal(np.exp((20 - 3) / (10000 / 11604.505)) - 1)

        self.assertAlmostEqual(expected, actual)

    def test_maxwell_boltzmann_distribution(self):
        actual = maxwell_boltzmann_distribution(self.energy, self.temperature)
        expected = np.reciprocal(np.exp(20 / (10000 / 11604.505)))

        self.assertAlmostEqual(actual, expected)


