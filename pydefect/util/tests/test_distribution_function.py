# -*- coding: utf-8 -*-

import numpy as np
import os
from pymatgen.util.testing import PymatgenTest

from pydefect.util.distribution_function import fermi_dirac_distribution, \
    bose_einstein_distribution, maxwell_boltzmann_distribution

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class DistributionTest(PymatgenTest):

    def test(self):
        energy = 20 # eV
        fermi_level = 3 # eV
        temperature = 10000

        actual_fd = fermi_dirac_distribution(energy, fermi_level, temperature)
        expected_fd = \
            np.reciprocal(np.exp((20 - 3) / (10000 / 11604.505)) + 1)

        self.assertAlmostEqual(actual_fd, expected_fd)

        actual_be = bose_einstein_distribution(energy, fermi_level, temperature)
        expected_be = \
            np.reciprocal(np.exp((20 - 3) / (10000 / 11604.505)) - 1)

        self.assertAlmostEqual(actual_be, expected_be)

        actual_mb = maxwell_boltzmann_distribution(energy, temperature)
        expected_mb = \
            np.reciprocal(np.exp(20 / (10000 / 11604.505)))

        self.assertAlmostEqual(actual_mb, expected_mb)


