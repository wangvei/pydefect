# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import physical_constants

EV = physical_constants['Boltzmann constant in eV/K'][0]


def fermi_dirac_distribution(energy: float,
                             fermi_level: float,
                             temperature: float) -> float:
    """Calculate Fermi-Dirac distribution at a given energy
    Args:
        energy (float): Calculated energy in eV.
        fermi_level (float): Fermi level in eV.
        temperature (float): Temperature in K.

    Return:
        calculated Fermi-Dirac distribution. For very tiny value, return 0.0
    """
    if (energy - fermi_level) / (EV * temperature) < 50:
        return np.reciprocal(
            np.exp((energy - fermi_level) / (EV * temperature)) + 1)
    else:
        return 0.0


def bose_einstein_distribution(energy: float,
                               fermi_level: float,
                               temperature: float) -> float:
    """Calculate Bose-Einstein distribution at a given energy
    Args:
        energy (float): Calculated energy in eV.
        fermi_level (float): Fermi level in eV.
        temperature (float): Temperature in K.

    Return:
        calculated Bose-Einstein distribution.
    """
    return np.reciprocal(
        np.exp((energy - fermi_level) / (EV * temperature)) - 1)


def maxwell_boltzmann_distribution(energy: float, temperature: float) -> float:
    """Calculate Maxwell-Boltzmann distribution at a given energy
    Args:
        energy (float): Calculated energy in eV.
        temperature (float): Temperature in K.

    Return:
        calculated Maxwell-Boltzmann distribution.
    """
    return np.exp(- energy / (EV * temperature))


