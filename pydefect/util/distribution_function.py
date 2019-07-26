# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import physical_constants

EV = physical_constants['Boltzmann constant in eV/K'][0]


def fermi_dirac_distribution(energy, fermi_level, temperature):
    """Calculate Fermi-Dirac distribution at a given energy"""
    if (energy - fermi_level) / (EV * temperature) < 50:
        return np.reciprocal(
            np.exp((energy - fermi_level) / (EV * temperature)) + 1)
    else:
        return 0.0


def bose_einstein_distribution(energy, fermi_level, temperature):
    return np.reciprocal(
        np.exp((energy - fermi_level) / (EV * temperature)) - 1)


def maxwell_boltzmann_distribution(energy, temperature):
    return np.exp(- energy / (EV * temperature))


