# -*- coding: utf-8 -*-

from collections import defaultdict
import math
from matplotlib import pyplot as plt
import numpy as np
from typing import List

from pydefect.analysis.defect import Defect
from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.core.error_classes import NoConvergenceError
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.distribution_function \
    import maxwell_boltzmann_distribution, fermi_dirac_distribution
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def hole_concentration(temperature: float,
                       e_f: float,
                       total_dos: np.array,
                       vbm: float,
                       volume: float,
                       threshold: float = 0.05):
    """ hole carrier concentration at the given absolute fermi_level.

    Args:
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        total_dos (np.array):
            Total density of states as a function of absolute energy.
                [[energy1, energy2, ...], [dos1, dos2, ...]]
            The units of energy and dos are [eV] and [1/eV], respectively.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        volume (float):
            Volume in A^3.
        threshold (float):
            Determines up to where electrons occupy the dos from the vbm.
            Note that the total dos near the band edges is not so accurate
            compared to band structure, so threshold is needed.
    """
    mesh_distance \
        = (total_dos[0][-1] - total_dos[0][0]) / (len(total_dos[0]) - 1)
    # Note that e_f and e are opposite for holes.
    hole = sum(fermi_dirac_distribution(e_f, e, temperature) * td
               for e, td in zip(total_dos[0], total_dos[1])
               if e <= vbm + threshold)

    return hole * mesh_distance / (volume / 10 ** 24)


def electron_concentration(temperature: float,
                           e_f: float,
                           total_dos: np.array,
                           cbm: float,
                           volume: float,
                           threshold: float = 0.05):
    """ electron carrier concentration at the given absolute fermi_level.

    Args:
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        total_dos (np.array):
            Total density of states as a function of absolute energy.
                [[energy1, energy2, ...], [dos1, dos2, ...]]
        cbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        volume (float):
            Volume in cm-3.
        threshold (float):
            Determines up to where electrons occupy the dos from the vbm.
            Note that the total dos near the band edges is not so accurate
            compared to band structure, so threshold is needed.
    """
    mesh_distance \
        = (total_dos[0][-1] - total_dos[0][0]) / (len(total_dos[0]) - 1)
    electron = sum(fermi_dirac_distribution(e, e_f, temperature) * td
                   for e, td in zip(total_dos[0], total_dos[1])
                   if e >= cbm - threshold)

    return electron * mesh_distance / (volume / 10 ** 24)


def calc_concentration(defect_energies: dict,
                       temperature: float,
                       e_f: float,
                       vbm: float,
                       cbm: float,
                       total_dos: np.array,
                       multiplicity: dict,
                       magnetization: dict,
                       volume: float,
                       reference_concentration: dict = None):
    """ Calculate carrier/defect concentration at given temp & Fermi level

    When the reference_concentration is provided, each defect specie
    concentration is fixed.

    Args:
        defect_energies (dict):
            Defect formation energies. energies[name][charge][annotation]
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (2xN numpy array):
                [[energy1, energy2, ...], [dos1, dos2, ...]]
        multiplicity (dict):
            Multiplicity in the supercell. It depends on the number of sites
            in the supercell and the site symmetry.
        volume (float):
            Volume in A-3.
        magnetization (dict):
            Magnetization in mu_B. magnetization[name][charge]
        reference_concentration (dict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench.
            reference_concentration[name][charge][annotation]

    Return:
        concentration (dict):
            Defect concentration in cm-3. concentration[name][charge]
        Note that "p" and "n" are special and mean the carrier hole and
        electron concentration in cm-3.
    """
    concentration = defaultdict(dict)

    for name in defect_energies.keys():
        for charge in defect_energies[name].keys():
            for annotation, de in defect_energies[name][charge].items():

                mul = multiplicity[name][charge][annotation]
                mag = magnetization[name][charge][annotation]

                if abs(mag - round(mag)) > 0.1:
                    logger.warning(
                        "The total_magnetization of {} in {} is {}, and not "
                        "integer".format(name, charge, mag))

                num_mag_conf = round(mag) + 1
                degree_of_freedom = mul * num_mag_conf

                # volume unit conversion from [A^3] to [cm^3]
                energy = de + e_f * charge
                c = (maxwell_boltzmann_distribution(energy, temperature)
                     * degree_of_freedom / (volume / 10 ** 24))

                concentration[name].setdefault(charge, {}).update({annotation: c})
#            if reference_concentration is None:
#            else:
#                dct_sum = {name: sum(v2.values()) for k, v2 in v1.items() for name, v1 in reference_concentration.items()}
            #     total = sum(reference_concentration[name].values())
            #     c_sum = sum(cc.values())

                # if math.isnan(total) or math.isnan(c_sum):
                #     raise ValueError("Carrier/defect concentration is too small.")

                # for charge in defect_energies[name].keys():
                #     concentration[name][charge] = total * charge[charge] / c_sum

    concentration["p"][1] = \
        {None: hole_concentration(temperature, e_f, total_dos, vbm, volume)}
    concentration["n"][-1] = \
        {None: electron_concentration(temperature, e_f, total_dos, cbm, volume)}

    return dict(concentration)


def calc_equilibrium_concentration(defect_energies: dict,
                                   temperature: float,
                                   vbm: float,
                                   cbm: float,
                                   total_dos: np.array,
                                   multiplicity: dict,
                                   magnetization: dict,
                                   volume: float,
                                   reference_concentration: dict = None,
                                   verbose: bool = False,
                                   max_iteration: int = 50,
                                   interval_decay_parameter: float = 0.7,
                                   threshold: float = 1e-5):
    """ Calculates equilibrium carrier & defect concentration at a temperature

    Args:
        defect_energies (dict):
            Defect formation energies. defect_energies[name][charge]
        temperature (float):
            Temperature in K.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (2xN numpy array):
            Total density of states.
            total_dos[[energy1, dos1], [energy2, dos2],...]
        multiplicity (dict):
            Multiplicity in the supercell. It depends on the site symmetry.
        magnetization (dict):
            Magnetization in mu_B. total_magnetization[defect][charge]
        volume (float):
            Volume in A-3.
        reference_concentration (dict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench. original_concentration[name][charge]
        verbose (bool):
            Whether print the carrier and defect concentration during the seek
            of the selfconsistent answer.
        max_iteration (int):
            Max iterations for the seek of the equilibrium defect concentration
        interval_decay_parameter (float):
            Determines how the interval decays.
            In case the Fermi level locates in between vbm and cbm, the common
            ratio of 0.5 is okay. Otherwise, higher common ratio is needed.
            Therefore, 0.7 is set here.
        threshold (float):
            Threshold for the convergence of the selfconsistent calculation of
            the carrier and defect concentration, which is a ratio of the
            residual charge sum and the highest carrier or defect concentration.
    """
    e_f = (cbm - vbm) / 2
    interval = e_f

    for iteration in range(max_iteration):
        defect_concentration = \
            calc_concentration(
                defect_energies=defect_energies, temperature=temperature,
                e_f=e_f, vbm=vbm, cbm=cbm, total_dos=total_dos,
                multiplicity=multiplicity, magnetization=magnetization,
                volume=volume, reference_concentration=reference_concentration)

        holes = defect_concentration["p"][1][None]
        electrons = defect_concentration["n"][-1][None]

        total_charge = \
            sum([e * c for d in defect_concentration
                 for c in defect_concentration[d]
                 for e in defect_concentration[d][c].values()])

        if verbose:
            logger.info("- {}th iteration ----".format(iteration))
            logger.info("Fermi level: {:.2f} eV.".format(e_f))

            for name, c_of_charge in defect_concentration.items():
                for charge, c_of_annotation in c_of_charge.items():
                    for annotation, concentration in c_of_annotation.items():
                        logger.info("{:>8} {:2d} {:>6}: {:.1e} cm-3.".
                                    format(name, charge, str(annotation),
                                           concentration))

            logger.info("Charge sum: {:.1e} cm-3.".format(total_charge))

        interval *= interval_decay_parameter
        e_f = e_f + np.sign(total_charge) * interval
        # This line controls the accuracy.
        max_concentration = np.amax([electrons, holes, total_charge])

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
                 defect_energies: dict,
                 multiplicity: dict,
                 magnetization: dict,
                 volume: float,
                 vbm: float,
                 cbm: float,
                 total_dos: np.array):

                 # concentrations: list,
                 # equilibrium_ef: list,
                 # equilibrium_concentrations: list,
                 # quenched_temperature: float = None,
                 # quenched_equilibrium_ef: list = None,
                 # quenched_equilibrium_concentrations: list = None):
        """
        Args:
            defect_energies (dict):
                DefectEnergy as a function of name, charge, and annotation.
                defect_energies[name][charge][annotation] = DefectEnergy object
            multiplicity (dict):
                Spatial multiplicity as a function of name, charge,
                and annotation.
                multiplicity[name][charge][annotation] = int
            magnetization (dict):
                Magnetization as a function of name, charge, and annotation.
                magnetization[name][charge][annotation] = float
            volume (float):
                Volume in A^3.
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            total_dos (np.array):
                Total density of states.

        Attributes:
            temp (float):
                Temperature in K.
            concentrations (list):
                Carrier and defect concentrations in cm-3.
                concentrations[Fermi][name][charge]
                * Carrier electron: concentrations[Fermi]["e"][-1]
                * Carrier hole: concentrations[Fermi]["p"][1]
            ef (float):
                Equilibrium Fermi level at each temperature.
            equilibrium_concentrations (list):
                Equilibrium concentration at each temperature as functions of
                name and charge.
                equilibrium_concentrations[name][charge] = float
            quenched_temp (float):
                Temperature in K.
            quenched_equilibrium_ef (float):
                Equilibrium Fermi level at quenched_temp quenched from temp.
            quenched_equilibrium_concentrations (list):
                Equilibrium concentration at quenched_temp quenched from temp.
                quenched_equilibrium_concentrations[name][charge] = float
        """

        self.defect_energies = defect_energies
        self.multiplicity = multiplicity
        self.magnetization = magnetization
        self.volume = volume
        self.vbm = vbm
        self.cbm = cbm
        self.total_dos = total_dos

        self.temp = None
        self.concentrations = None

        self.equilibrium_ef = None
        self.equilibrium_concentrations = None

        self.quenched_temp = None
        self.quenched_ef = None
        self.quenched_concentrations = None

    @classmethod
    def from_calc_results(cls,
                          defect_energies_obj: DefectEnergies,
                          unitcell: UnitcellCalcResults):
#                          verbose: bool = False):
        #                          temperatures: list,
        # fermi_range: list = None,
        # fermi_mesh_num: int = 101,
        """ Calculates defect formation energies from some files.

        Args:
            temperatures (list):
                List of temperatures in K.
            defect_energies_obj (DefectEnergies):
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
        defect_energies = defect_energies_obj.defect_energies
        multiplicity = defect_energies.multiplicity
        magnetization = defect_energies.magnetization

        volume = unitcell.volume  # [A^3]
        total_dos = unitcell.total_dos
        vbm = defect_energies_obj.vbm
        cbm = defect_energies_obj.cbm

        return cls(defect_energies=defect_energies,
                   multiplicity=multiplicity,
                   magnetization=magnetization,
                   volume=volume, vbm=vbm, cbm=cbm, total_dos=total_dos)


    # def _repr__(self):
    #     outs = ["--------",
    #             "Temperature: {} K.".format(self.temp),
    #             "Fermi level: {} eV.".format(self.e_f),
    #             "p: {:.1e} cm-3.".format(self.p),
    #             "n: {:.1e} cm-3.".format(self.n),
    #             "p - n: {:.1e} cm-3.".format(self.p - self.n)]

        # for name, c_of_charge in self.concentration.items():
        #     outs.append("---")
        #     for charge, concent in c_of_charge.items():
        #         outs.append("{} {}: {:.1e} cm-3.".format(name, charge, concent))

        # return "\n".join(outs)


    def calc_concentration(self, temp, fermi_range, num_mesh):
        concentrations = {}
        if fermi_range is None:
            fermi_range = [self.vbm, self.cbm]
        fermi_mesh = np.linspace(fermi_range[0], fermi_range[1], num_mesh)

        for e_f in fermi_mesh:
            concentrations[e_f] = \
                calc_concentration(
                    defect_energies=self.defect_energies,
                    temperature=temp, e_f=e_f, vbm=self.vbm, cbm=self.cbm,
                    total_dos=self.total_dos, multiplicity=self.multiplicity,
                    magnetization=self.magnetization, volume=self.volume)


        # if e_f:
        #     p, n, c = calc_defect_carrier_concentration(
        #         energies, temperature, e_f, vbm, cbm, total_dos, multiplicity,
        #         volume, magnetization, original_concentration)
        # else:
        #     e_f, p, n, c = calc_equilibrium_concentration(
        #         energies, temperature, vbm, cbm, total_dos, multiplicity, volume,
        #         magnetization, original_concentration, verbose)

        return


    def calc_quenched_concentration(self, temp):
        pass


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

    # def __repr__(self):
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
        volume = unitcell.volume
        total_dos = unitcell.total_dos  # [eV] and [1/eV]
        vbm, cbm = unitcell.band_edge

        if e_range is None:
            e_range = [vbm - 1, cbm + 1]

        num_mesh = (cbm - vbm) / 0.01
        fermi_levels = list(np.linspace(e_range[0], e_range[1], num_mesh))

        ns = []
        ps = []

        for f in fermi_levels:
            ns.append(electron_concentration(temperature, f, total_dos, cbm, volume))
            ps.append(hole_concentration(temperature, f, total_dos, vbm, volume))

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
