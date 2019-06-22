# -*- coding: utf-8 -*-

import math
import numpy as np
import warnings

from matplotlib import pyplot as plt
from monty.serialization import loadfn

from collections import defaultdict

from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.core.error_classes import NoConvergenceError
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.distribution_function import maxwell_boltzmann_distribution, \
    fermi_dirac_distribution
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def calc_defect_carrier_concentration(energies: dict,
                                      temperature: float,
                                      e_f: float,
                                      vbm: float,
                                      cbm: float,
                                      total_dos: dict,
                                      multiplicity: dict,
                                      magnetization: dict,
                                      volume: float,
                                      reference_concentration: dict = None):
    """ Calculate carrier & defect concentration at a temperature & Fermi level

    When the reference_concentration is provided, each defect specie
    concentration is fixed.

    Args:
        energies (dict):
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
        multiplicity (dict):
            Multiplicity in the supercell. It depends on the site symmetry.
        volume (float):
            Volume in cm-3.
        magnetization (dict):
            Magnetization in mu_B. magnetization[name][charge]
        reference_concentration (dict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench. reference_concentration[name][charge]
    Return:
        p (float):
            Hole concentration in cm-3.
        n (float):
            Electron concentration in cm-3.
        concentration (dict):
            Defect concentration in cm-3. concentration[name][charge]
    """
    concentration = defaultdict(dict)

    for name in energies.keys():

        c = {}
        for charge in energies[name].keys():

            m = magnetization[name][charge]

            if abs(m - round(m)) > 0.01:
                logger.warning(
                    "The total_magnetization of {} in {} is {}, and not "
                    "almost integer".format(name, charge, m))

            num_mag_conf = round(m) + 1
            degree_of_freedom = multiplicity[name][charge] * num_mag_conf

            energy = energies[name][charge] + e_f * charge
            c[charge] = maxwell_boltzmann_distribution(energy, temperature) \
                        / volume * degree_of_freedom

        if reference_concentration is None:
            concentration[name] = c
        else:
            total = sum(reference_concentration[name].values())
            c_sum = sum(c.values())

            if math.isnan(total) or math.isnan(c_sum):
                raise ValueError("Carrier/defect concentration is too small.")

            for charge in energies[name].keys():
                concentration[name][charge] = total * c[charge] / c_sum

    # Use absolute Fermi level for p and n.
    concentration["p"][1] = p(temperature, e_f, total_dos, vbm, volume)
    concentration["n"][-1] = n(temperature, e_f, total_dos, cbm, volume)

    return dict(concentration)


def calc_equilibrium_concentration(energies,
                                   temperature,
                                   vbm,
                                   cbm, total_dos,
                                   multiplicity, volume, magnetization,
                                   reference_concentration, verbose=False,
                                   max_iteration=50, threshold=1e-5):
    """ Calculates equilibrium carrier & defect concentration at a temperature

    Args:
        energies (dict):
            Defect formation energies. energies[name][charge]
        temperature (float):
            Temperature in K.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (dict):
            Total density of states.
        multiplicity (dict):
            Multiplicity in the supercell. It depends on the site symmetry.
        volume (float):
            Volume in cm-3.
        magnetization (dict):
            Magnetization in mu_B. total_magnetization[defect][charge]
        reference_concentration (dict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench. original_concentration[name][charge]
        max_iteration (int):
            Max iterations for the seek of the equilibrium defect concentration
        threshold (float):
            Threshold for the convergence of the selfconsistent calculation of
            the carrier and defect concentration, which is a ratio of the
            residual charge sum and the highest carrier or defect concentration.
        verbose (bool):
            Whether print the carrier and defect concentration during the seek
            of the selfconsistent answer.
    """
    e_f = (cbm - vbm) / 2
    interval = e_f

    for iteration in range(max_iteration):
        defect_concentration = calc_defect_carrier_concentration(
                energies, temperature, e_f, vbm, cbm, total_dos, multiplicity,
                magnetization, volume, reference_concentration)

        total_charge = \
            sum([sum([c * e for c, e in defect_concentration[d].items()])
                 for d in defect_concentration])

        if verbose:
            logger.info("- {} ----".format(iteration))
            logger.info("Fermi level: {} eV.".format(e_f))
            logger.info("p: {:.1e} cm-3.".format(p))
            logger.info("n: {:.1e} cm-3.".format(n))

            for name, c_of_charge in defect_concentration.items():
                for charge, concentration in c_of_charge.items():
                    logger.info("{} {}: {:.1e} cm-3.".
                                format(name, charge, concentration))

            logger.info("Charge sum: {:.1e} cm-3.".format(total_charge))

        # In case the Fermi level locates in between vbm and cbm, the common
        # ratio of 0.5 is okay. Otherwise, higher common ratio is needed.
        # Therefore, 0.7 is set here.
        interval *= 0.7
        e_f = e_f + np.sign(total_charge) * interval
        # This line controls the accuracy.
        max_concentration = np.amax([n, p, total_charge])
        if np.abs(total_charge / max_concentration) < threshold:
            return e_f, defect_concentration

    raise NoConvergenceError("Equilibrium condition has not been reached.")


# initial_num_symmops = \
#     num_symmetry_operation(d.defect_entry.initial_site_symmetry)
# final_num_symmops = \
#     num_symmetry_operation(d.dft_results.site_symmetry)

# mul = int(d.defect_entry.num_equiv_sites /
#           final_num_symmops * initial_num_symmops)

class DefectConcentration:
    """ A class related to a set of carrier and defect concentration """

    def __init__(self,
                 concentrations: list,
                 equilibrium_efs: list,
                 equilibrium_concentrations: list,
                 vbm: float,
                 cbm: float,
                 quenched_temperature: float = 298,
                 quenched_equilibrium_efs: list = None,
                 quenched_equilibrium_concentrations: list = None):
        """
        Args:
            concentrations (list):
                Carrier and defect concentration in cm-3.
                concentrations[temp][Fermi][name][charge]
                * Carrier electron: concentrations[temp][Fermi][e][-1]
                * Carrier hole: concentrations[temp][Fermi][p][1]
            equilibrium_efs (list):
                Equilibrium Fermi level at each temperature.
            equilibrium_concentrations (list):
                Equilibrium concentration at each temperature.
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            quenched_equilibrium_efs (list):
                Equilibrium Fermi level at quenched_temperature quenched from
                temperature.
            quenched_equilibrium_concentrations (list):
                Equilibrium concentration at quenched_temperature quenched from
                temperature.
            quenched_temperature (list):
                Temperature in K.
        """

        self.concentrations = concentrations
        self.equilibrium_efs = equilibrium_efs
        self.equilibrium_concentrations = equilibrium_concentrations
        self.vbm = vbm
        self.cbm = cbm

        self.quenched_temperature = quenched_temperature
        self.quenched_equilibrium_efs = quenched_equilibrium_efs
        self.quenched_equilibrium_concentrations = \
            quenched_equilibrium_concentrations

    @classmethod
    def from_calc_results(cls,
                          temperatures: list,
                          defect_energies: DefectEnergies,
                          unitcell: UnitcellCalcResults,
                          fermi_range: list = None,
                          fermi_mesh_num: int = 101,
                          verbose: bool = False):
        """ Calculates defect formation energies from some files.

        Args:
            temperatures (list):
                List of temperatures in K.
            defect_energies (DefectEnergies):
                DefectEnergies object used for calculating concentration.
            unitcell (UnitcellCalcResults):
                UnitcellDftResults object for volume and total_dos
            fermi_range (list):
                Range of Fermi level, [lower_limit, upper_limit]
            fermi_mesh_num (float):
                Number of mesh for Fermi level including boundary.
            verbose (book):
                Whether print the carrier and defect concentration during the
                seek of the selfconsistent results.
        """
        energies = defect_energies.energies
        multiplicity = defect_energies.multiplicity
        magnetization = defect_energies.magnetization

        volume = unitcell.volume * 10 ** -24  # [A^3] -> [cm^3]
        total_dos = unitcell.total_dos
        vbm = defect_energies.vbm
        cbm = defect_energies.cbm

        concentrations = defaultdict(dict)
        if fermi_range is None:
            fermi_range = [vbm, cbm]
        fermi_mesh = np.linspace(fermi_range[0], fermi_range[1], fermi_mesh_num)

        for t in temperatures:
            for f in fermi_mesh:
                concentrations[t][f] = \
                    calc_defect_carrier_concentration(
                        energies, t, f, vbm, cbm, total_dos, multiplicity,
                        magnetization, volume)


        # if e_f:
        #     p, n, c = calc_defect_carrier_concentration(
        #         energies, temperature, e_f, vbm, cbm, total_dos, multiplicity,
        #         volume, magnetization, original_concentration)
        # else:
        #     e_f, p, n, c = calc_equilibrium_concentration(
        #         energies, temperature, vbm, cbm, total_dos, multiplicity, volume,
        #         magnetization, original_concentration, verbose)

        # return cls(energies, temperature, e_f, p, n, previous_concentration, c)

    def __str__(self):
        outs = ["--------",
                "Temperature: {} K.".format(self.temperature),
                "Fermi level: {} eV.".format(self.e_f),
                "p: {:.1e} cm-3.".format(self.p),
                "n: {:.1e} cm-3.".format(self.n),
                "p - n: {:.1e} cm-3.".format(self.p - self.n)]

        for name, c_of_charge in self.concentration.items():
            outs.append("---")
            for charge, concent in c_of_charge.items():
                outs.append("{} {}: {:.1e} cm-3.".format(name, charge, concent))

        return "\n".join(outs)


class CarrierConcentration:
    """ Container class for carrier concentration """

    def __init__(self,
                 temperature: float,
                 vbm: float,
                 cbm: float,
                 fermi_levels: list,
                 ns: list,
                 ps: list):
        """
        Args:
            temperature (float): list of temperature considered.
            fermi_levels: Fermi levels considered in the absolute scale in eV.
            ns[temperature][Fermi level]: Carrier electron concentrations.
            ps[temperature][Fermi level]: Carrier hole concentrations.
        """
        self.temperature = temperature
        self.vbm = vbm
        self.cbm = cbm
        self.fermi_levels = fermi_levels
        self.ns = ns
        self.ps = ps

    @property
    def electron_concentration(self):
        return self.ns

    @property
    def hole_concentration(self):
        return self.ps

    # def __str__(self):
    #     outs = ["Temperature [K]: " + str(self.temperatures),
    #             "E_f [eV],  n [cm-3],  p [cm-3]"]
    #     for a, b, c in zip(self.fermi_levels, self.ns, self.ps):
    #         outs.append('%8.4f'%a + "   " + '%.2e'%b + "   " + '%.2e'%c)

        # return "\n".join(outs)

    @classmethod
    def from_unitcell(cls, temperature, unitcell, e_range=None):
        """ Calculates carrier concentration.

        Args:
            temperature (float):
            unitcell (UnitcellCalcResults):
            e_range (list):
        """
        volume = unitcell.volume * 10 ** -24  # [A^3] -> [cm^3]
        total_dos = unitcell.total_dos  # [1/eV] and [eV]
        vbm, cbm = unitcell.band_edge

        if e_range is None:
            e_range = [vbm - 1, cbm + 1]

        num_mesh = (cbm - vbm) / 0.01
        fermi_levels = list(np.linspace(e_range[0], e_range[1], num_mesh))

        ns = []
        ps = []

        for f in fermi_levels:

            ns.append(n(temperature, f, total_dos, cbm, volume))
            ps.append(p(temperature, f, total_dos, vbm, volume))

        return cls(temperature, vbm, cbm, fermi_levels, ns, ps)

    def get_plot(self, xlim=None, ylim=None, relative=True):
        """ Get a matplotlib plot.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            relative:
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.title("Temperature:" + str(self.temperature) + " K")

        ax.set_xlabel("Fermi level (eV)")
        ax.set_ylabel("Concentration (cm-3)")
        ax.set_yscale("log", nonposy='clip')

        max_y = max([max(self.ns), max(self.ps)])

        if xlim:
            plt.xlim(xlim[0], xlim[1])

        if ylim:
            plt.ylim(10 ** ylim[0], 10 ** ylim[1])
        else:
            plt.ylim(10**10, max_y)

        if relative is True:
            vbm = 0.0
            cbm = self.cbm - self.vbm
            fermi_levels = [f - self.vbm for f in self.fermi_levels]
        else:
            vbm = self.vbm
            cbm = self.cbm
            fermi_levels = self.fermi_levels

        plt.axvline(x=vbm, linewidth=1.0, linestyle='dashed')
        plt.axvline(x=cbm, linewidth=1.0, linestyle='dashed')

        ax.plot(fermi_levels, self.ns, '-', color="red", label="n")
        ax.plot(fermi_levels, self.ps, '-', color="blue", label="p")

        plt.show()


def p(temperature, fermi_level, total_dos, vbm, volume, threshold=0.05):
    """ hole carrier concentration at the given absolute fermi_level. """
    mesh_distance = total_dos[1][1] - total_dos[1][0]
    hole = sum(fermi_dirac_distribution(fermi_level, e, temperature) * td
            for td, e in zip(total_dos[0], total_dos[1])
            if e <= vbm + threshold)

    return hole * mesh_distance / volume


def n(temperature, fermi_level, total_dos, cbm, volume, threshold=0.05):
    """ electron carrier concentration at the given absolute fermi_level."""
    mesh_distance = total_dos[1][1] - total_dos[1][0]
    electron = sum(fermi_dirac_distribution(e, fermi_level, temperature) * td
            for td, e in zip(total_dos[0], total_dos[1])
            if e >= cbm - threshold)
    return electron * mesh_distance / volume