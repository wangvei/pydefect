#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import physical_constants

EV = physical_constants['Boltzmann constant in eV/K'][0]


def fermi_dirac_dist(energy, fermi_level, temperature):

    if (energy - fermi_level) / (EV * temperature) < 100:
        return np.reciprocal(
            np.exp((energy - fermi_level) / (EV * temperature)) + 1)
    else:
        return 0.0


def bose_einstein_dist(energy, fermi_level, temperature):
    return np.reciprocal(
        np.exp((energy - fermi_level) / (EV * temperature)) - 1)


def maxwell_boltzmann_dist(energy, temperature):
    return np.exp(- energy / (EV * temperature))


