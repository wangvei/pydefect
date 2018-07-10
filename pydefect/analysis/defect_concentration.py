# -*- coding: utf-8 -*-

import math
import numpy as np
import json
import warnings

from monty.json import MontyEncoder
from monty.serialization import loadfn

from collections import defaultdict, namedtuple
from itertools import combinations

from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.input_maker.defect_set_maker import is_name_selected
from pydefect.util.carrier_concentration import maxwell_boltzmann_dist, \
    CarrierConcentration


def calc_concentration(energies, temperature, e_f, vbm, cbm, total_dos,
                       num_sites, volume, magnetization,
                       original_concentration):
    """
    Calculates the carrier and defect concentration at a given temperature and
    Fermi level. When original_concentration is given, defect concentration at
    each defect type is fixed. For example, sum of Va_O1_0, Va_O1_1 and Va_O1_2
    is fixed.

    Args:
        energies (defaultdict):
            Defect formation energies. energies[name][charge]
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (dict):
            Total density of states.
        num_sites (dict):
            Number of sites. num_sites[name]
        volume (float):
            Volume in cm-3.
        magnetization (dict):
            Magnetization in \mu_B. magnetization[defect][charge]
        original_concentration (defaultdict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench. original_concentration[name][charge]
    Return:
        p (float):
            Hole concentration in cm-3.
        n (float):
            Electron concentration in cm-3.
        concentration (defaultdict):
            Defect concentration in cm-3. concentration[defect name][charge]
    """

    concentration = defaultdict(dict)
    for name in energies.keys():
        if name not in num_sites.keys():
            raise KeyError(
                "Number of sites for {} is not given.".format(name))

        c = {}
        for charge in energies[name].keys():

            m = magnetization[name][charge]
            if abs(m - round(m)) > 0.01:
                warnings.warn("The magnetization of {} in {} is {}, and not "
                              "almost integer".format(name, charge, m))
            num_mag_conf = round(m) + 1
            degree_of_freedom = num_sites[name] * num_mag_conf

            energy = energies[name][charge] + e_f * charge
            c[charge] = maxwell_boltzmann_dist(energy, temperature) / volume \
                        * degree_of_freedom

        if original_concentration is None:
            concentration[name] = c
        else:
            total = sum(original_concentration[name].values())
            c_sum = sum(c.values())

            if math.isnan(total) or math.isnan(c_sum):
                raise ValueError("Carrier/defect concentration is too small.")

            for charge in energies[name].keys():
                concentration[name][charge] = total * c[charge] / c_sum
    # Use absolute Fermi level for p and n.
    e_f_abs = e_f + vbm
    p = CarrierConcentration.p(temperature, e_f_abs, total_dos, vbm, volume)
    n = CarrierConcentration.n(temperature, e_f_abs, total_dos, cbm, volume)

    return p, n, concentration


def calc_equilibrium_concentration(energies, temperature, vbm, cbm, total_dos,
                                   num_sites, volume, magnetization,
                                   original_concentration, verbose=False,
                                   max_iteration=50, threshold=1e-5):
    """
    Calculates the equilibrium carrier and defect concentration at a given
    temperature.

            Defect formation energies. energies[name][charge]
        temperature (float):
            Temperature in K.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (dict):
            Total density of states.
        num_sites (dict):
            Number of sites. num_sites[name]
        volume (float):
            Volume in cm-3.
        magnetization (dict):
            Magnetization in \mu_B. magnetization[defect][charge]
        original_concentration (defaultdict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench. original_concentration[name][charge]
        max_iteration (int):
            Max iterations for the seek of the equilibrium defect concentration
        threshold (float):
            Threshold for the convergence of the selfconsistent calculation of
            the carrier and defect concentration, which is a ratio of the
            residual charge sum and the highest carrier or defect concentration.
        verbose (book):
            Whether print the carrier and defect concentration during the seek
            of the selfconsistent answer.
    """
    e_f = (cbm - vbm) / 2
    interval = e_f
    for iteration in range(max_iteration):
        p, n, defect_concentration = \
            calc_concentration(energies, temperature, e_f, vbm, cbm, total_dos,
                               num_sites, volume, magnetization,
                               original_concentration)
        total_defect_charge = \
            sum([sum([c * e for c, e in defect_concentration[d].items()])
                 for d in defect_concentration])
        charge_sum = p - n + total_defect_charge

        if verbose:
            print("- {} ----".format(iteration))
            print("Fermi level: {} eV.".format(e_f))
            print("p: {:.1e} cm-3.".format(p))
            print("n: {:.1e} cm-3.".format(n))
            for name, c_of_charge in defect_concentration.items():
                for charge, concent in c_of_charge.items():
                    print("{} {}: {:.1e} cm-3.".format(name, charge, concent))
            print("Charge sum: {:.1e} cm-3.".format(charge_sum))

        # In case the Fermi level locates in between vbm and cbm, the common
        # ratio 0.5 is okay. Otherwise, higher common ratio is needed.
        # Therefore, 0.7 is set here.
        interval *= 0.7
        e_f = e_f + np.sign(charge_sum) * interval
        # This line controls the accuracy.
        max_concentration = np.amax([n, p, charge_sum])
        if np.abs(charge_sum / max_concentration) < threshold:
            return e_f, p, n, defect_concentration

    raise NoConvergenceError("Equilibrium condition has not been reached.")


class DefectConcentration:
    """
    A class related to a set of carrier and defect concentration
    """

    def __init__(self, energies, temperature, e_f, p, n,
                 previous_concentration, concentration):
        """
        Args:
            energies (defaultdict):
                Defect formation energies. energies[name][charge]
            temperature (float):
                Temperature in K.
            e_f (float):
                Fermi level in the absolute scale.
            p (float):
                Hole concentration in cm-3.
            n (float):
                Electron concentration in cm-3.
            previous_concentration (DefectConcentration):
                Concentration with previous_concentration[name][charge] in cm-3
                that is used for calculating the concentration at the quenching
                temperature.
            concentration (defaultdict):
                Defect concentration in cm-3. concentration[defect name][charge]
        """

        self._energies = energies
        self._temperature = temperature
        self._e_f = e_f
        self._p = p
        self._n = n
        self._previous_concentration = previous_concentration
        self._concentration = concentration

    @classmethod
    def from_defect_energies(cls, defect_energies, temperature, unitcell,
                             num_sites_filename, e_f=None,
                             previous_concentration=None, verbose=False):
        """
        Calculates defect formation energies from some files.
        Args:
            defect_energies (DefectEnergies):
                DefectEnergies object used for calculating concentration.
            temperature (float):
                Temperature in K.
            unitcell (UnitcellDftResults):
                UnitcellDftResults object for volume and total_dos
            num_sites_filename (str):
                Yaml file name for the number of sites in a given volume like
                Va_Mg1 : 2
            e_f (float):
                Fermi level in the absolute scale. When it is given, the
                concentration is calculated at the Fermi level.
            previous_concentration (DefectConcentration):
                Concentration with previous_concentration[name][charge] form.
                When it is given, the total concentration at each defect is
                preserved.
            verbose (book):
                Whether print the carrier and defect concentration during the seek
                of the selfconsistent answer.
        """
        energies = defect_energies.energies
        volume = unitcell.volume * 10 ** -24  # [A^3] -> [cm^3]
        total_dos = unitcell.total_dos
        vbm = defect_energies.vbm
        cbm = defect_energies.cbm
        magnetization = defect_energies.magnetization
        num_sites = loadfn(num_sites_filename)

        if previous_concentration is None:
            original_concentration = None
        else:
            original_concentration = previous_concentration.concentration

        if e_f:
            p, n, c = calc_concentration(
                energies, temperature, e_f, vbm, cbm, total_dos, num_sites,
                volume, magnetization, original_concentration)
        else:
            e_f, p, n, c = calc_equilibrium_concentration(
                energies, temperature, vbm, cbm, total_dos, num_sites, volume,
                magnetization, original_concentration, verbose)

        return cls(energies, temperature, e_f, p, n, previous_concentration, c)

    def __str__(self):
        print("--------")
        print("Temperature: {} K.".format(self.temperature))
        print("Fermi level: {} eV.".format(self.e_f))
        print("p: {:.1e} cm-3.".format(self.p))
        print("n: {:.1e} cm-3.".format(self.n))
        print("p - n: {:.1e} cm-3.".format(self.p - self.n))
        for name, c_of_charge in self.concentration.items():
            print("---")
            for charge, concent in c_of_charge.items():
                print("{} {}: {:.1e} cm-3.".format(name, charge, concent))


    @property
    def energies(self):
        return self._energies

    @property
    def temperature(self):
        return self._temperature

    @property
    def e_f(self):
        return self._e_f

    @property
    def p(self):
        return self._p

    @property
    def n(self):
        return self._n

    @property
    def previous_concentration(self):
        return self._previous_concentration

    @property
    def concentration(self):
        return self._concentration


class NoConvergenceError(Exception):
    pass
