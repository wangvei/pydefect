# -*- coding: utf-8 -*-

from collections import defaultdict
from copy import deepcopy
from typing import Optional, Tuple

import numpy as np
from matplotlib import pyplot as plt
from monty.json import MSONable
from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.core.defect_name import DefectName
from pydefect.core.error_classes import NoConvergenceError
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.distribution_function \
    import maxwell_boltzmann_distribution, fermi_dirac_distribution
from pydefect.util.logger import get_logger
from pydefect.util.tools import (
    flatten_dict, defaultdict_to_dict, mod_defaultdict)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def hole_concentration(temperature: float,
                       e_f: float,
                       total_dos: np.array,
                       vbm: float,
                       volume: float,
                       threshold: float = 0.05) -> float:
    """ hole carrier concentration at the given absolute fermi_level.

    Args:
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        total_dos (np.array):
            Total density of states as a function of absolute energy.
                [[dos1, dos2, ...], [energy1, energy2, ...]]
            The units of energy and dos are [eV] and [1/eV], respectively.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        volume (float):
            Volume in A^3.
        threshold (float):
            Determines up to where electrons occupy the dos from the vbm.
            Note that the total dos near the band edges is not so accurate
            compared to band structure, so threshold is needed.

    Returns:
        Float of the hole concentration in cm-3.
    """
    doses, energies = total_dos
    energy_range = energies[-1] - energies[0]
    num_intervals = len(energies) - 1
    energy_interval = energy_range / num_intervals

    # Note that e_f and e are opposite for holes.
    hole = sum(fermi_dirac_distribution(e_f, e, temperature) * dos
               for dos, e in zip(doses, energies) if e <= vbm + threshold)

    return hole * energy_interval / (volume / 10 ** 24)


def electron_concentration(temperature: float,
                           e_f: float,
                           total_dos: np.array,
                           cbm: float,
                           volume: float,
                           threshold: float = 0.05) -> float:
    """ electron carrier concentration at the given absolute fermi_level.

    Args:
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        total_dos (np.array):
            Total density of states as a function of absolute energy.
                [[dos1, dos2, ...], [energy1, energy2, ...]]
        cbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        volume (float):
            Volume in A^3.
        threshold (float):
            Determines up to where electrons occupy the dos from the vbm.
            Note that the total dos near the band edges is not so accurate
            compared to band structure, so threshold is needed.

    Returns:
        Float of the electron concentration in cm-3.
    """
    doses, energies = total_dos

    energy_range = energies[-1] - energies[0]
    num_intervals = len(energies) - 1
    energy_interval = energy_range / num_intervals

    electron = sum(fermi_dirac_distribution(e, e_f, temperature) * dos
                   for dos, e in zip(doses, energies) if e >= cbm - threshold)

    return electron * energy_interval / (volume / 10 ** 24)


def calc_concentration(defect_energies: Optional[dict],
                       temperature: float,
                       e_f: float,
                       vbm: float,
                       cbm: float,
                       total_dos: np.array,
                       volume: float,
                       ref_concentration: Optional[dict] = None) -> dict:
    """ Calculate concentrations at given temperature & Fermi level

    When the reference_concentration is provided, each defect specie
    concentration is fixed.

    Args:
        defect_energies (dict):
            Defect formation energies. energies[name][charge]
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (2xN numpy array):
                [[dos1, dos2, ...], [energy1, energy2, ...]]
        volume (float):
            Volume in A-3.
        ref_concentration (dict):
            Original defect concentration in cm-3 used for calculating the
            concentration by quench.
            reference_concentration[name][charge][annotation]

    Return:
        concentration (dict):
            Defect concentration in cm-3. concentration[name][charge]
        Note that "p" and "n" have special meaning as the carrier hole and
        electron concentration in cm-3, respectively.
    """
    concentrations = mod_defaultdict(depth=2)

    concentrations["p"][1] = \
        {None: hole_concentration(temperature, e_f, total_dos, vbm, volume)}
    concentrations["n"][-1] = \
        {None: electron_concentration(temperature, e_f, total_dos, cbm, volume)}

    if defect_energies is None:
        return dict(concentrations)

    for name in defect_energies:
        concentration_by_name = defaultdict(dict)
        for charge, de in flatten_dict(defect_energies[name]):

            num_mag_conf = abs(de.magnetization) + 1
            degree_of_freedom = de.multiplicity * num_mag_conf

            energy = de.defect_energy + e_f * charge
            # volume unit conversion from [A^3] to [cm^3]
            c = (maxwell_boltzmann_distribution(energy, temperature)
                 * degree_of_freedom / (volume / 10 ** 24))

            concentration_by_name[charge] = {de.annotation: c}

        if ref_concentration:
            reference_total_concentration = \
                sum(sum(c.values()) for c in ref_concentration[name].values())
            total_concentration = \
                sum(sum(c.values()) for c in concentration_by_name.values())
            factor = (reference_total_concentration / total_concentration)

            for charge, annotation, _ \
                    in flatten_dict(concentration_by_name):
                concentration_by_name[charge][annotation] *= factor

        concentrations[name] = concentration_by_name

    return defaultdict_to_dict(concentrations)


def calc_equilibrium_concentration(defect_energies: dict,
                                   temperature: float,
                                   vbm: float,
                                   cbm: float,
                                   total_dos: np.array,
                                   volume: float,
                                   ref_concentration: Optional[dict] = None,
                                   verbose: bool = False,
                                   max_iteration: int = 200,
                                   interval_decay_parameter: float = 0.7,
                                   threshold: float = 1e-5
                                   ) -> Tuple[float, dict]:
    """ Calculates equilibrium carrier & defect concentration at a temperature

    When the ref_concentration is set, the total defect concentration with the
    same name is kept fixed.

    Args:
        defect_energies (dict):
            Defect formation energies. energies[name][charge]
        temperature (float):
            Temperature in K.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (2xN numpy array):
            Total density of states.
            total_dos[[energy1, dos1], [energy2, dos2],...]
        volume (float):
            Volume in A-3.
        ref_concentration (dict):
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
            residual charge sum and  highest carrier or defect concentration.
    """
    # Trial Fermi level begins at center of the band gap in the absolute scale.
    e_f = (cbm + vbm) / 2
    interval = (cbm - vbm) / 2
    for iteration in range(max_iteration):
        defect_concentration = \
            calc_concentration(defect_energies=defect_energies,
                               temperature=temperature,
                               e_f=e_f,
                               vbm=vbm,
                               cbm=cbm,
                               total_dos=total_dos,
                               volume=volume,
                               ref_concentration=ref_concentration)

        total_charge = \
            sum([energy * charge for name in defect_concentration
                 for charge in defect_concentration[name]
                 for energy in defect_concentration[name][charge].values()])

        if verbose:
            logger.info(f"- {iteration}th iteration ----")
            logger.info(f"Fermi level: {e_f:.2f} eV.")

            for name in defect_concentration:
                for charge in defect_concentration[name]:
                    for annotation in defect_concentration[name][charge]:
                        concentration = \
                            defect_concentration[name][charge][annotation]
                        logger.info(f"{name:>8}  {charge:2d}  {annotation:>6}:"
                                    f"   {concentration:.1e} cm-3.")

            logger.info(f"Charge sum: {total_charge:.1e} cm-3.")

        interval *= interval_decay_parameter
        e_f = e_f + np.sign(total_charge) * interval

        max_concentration = \
            np.amax([concentration for d in defect_concentration
                     for c in defect_concentration[d]
                     for concentration in defect_concentration[d][c].values()])

        # This part controls the accuracy.
        if np.abs(total_charge / max_concentration) < threshold:
            return e_f, defect_concentration

    raise NoConvergenceError("Equilibrium condition has not been reached.")


class DefectConcentration(MSONable):
    """ A class related to a set of carrier and defect concentration """

    def __init__(self,
                 defect_energies: dict,
                 volume: float,
                 vbm: float,
                 cbm: float,
                 total_dos: np.array,
                 temperature: float = None,
                 fermi_mesh: list = None,
                 concentrations: list = None,
                 equilibrium_ef: float = None,
                 equilibrium_concentration: dict = None,
                 quenched_temperature: float = None,
                 quenched_ef: float = None,
                 quenched_equilibrium_concentration: dict = None,
                 quenched_carrier_concentrations: list = None):
        """
        Args:
            defect_energies (dict):
                DefectEnergy as a function of name, charge, and annotation.
                energies[name][charge][annotation] = DefectEnergy object
            volume (float):
                Volume in A^3.
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            total_dos (np.array):
                Total density of states as a function of absolute energy.
                [[dos1, dos2, ...], [energy1, energy2, ...]]
            temperature (float):
                Temperature in K.
            fermi_mesh (list):
                List of calculated Fermi levels in the absolute scale.
            concentrations (list):
                A set of Carrier and defect concentrations in cm-3 as a
                function of the Fermi level
                concentrations[Fermi index][name][charge][annotation]
                * Carrier electron: concentrations[Fermi index]["e"][-1][None]
                * Carrier hole: concentrations[Fermi index]["p"][1][None]
            equilibrium_ef (float):
                Equilibrium Fermi level at temperature.
            equilibrium_concentration (dict):
                Equilibrium concentration at each temperature as functions of
                name and charge.
                equilibrium_concentrations[name][charge][annotation] = float
            quenched_temperature (float):
                Temperature in K.
            quenched_ef (float):
                Equilibrium Fermi level at quenched_temperature quenched from
                temperature.
            quenched_equilibrium_concentration (dict):
                Equilibrium concentration at quenched_temperature quenched from
                temperature.
                quenched_equilibrium_concentrations[name][charge][annotation]
            quenched_carrier_concentrations (list):
                A set of carrier concentrations as a function of the Fermi
                level at quenched_temperature
                * Carrier electron: concentrations[Fermi index]["e"][-1][None]
                * Carrier hole: concentrations[Fermi index]["p"][1][None]
        """
        self.defect_energies = defect_energies
        self.volume = volume
        self.vbm = vbm
        self.cbm = cbm
        self.total_dos = total_dos

        self.temperature = temperature
        self.fermi_mesh = fermi_mesh
        self.concentrations = concentrations

        self.equilibrium_ef = equilibrium_ef
        self.equilibrium_concentration = equilibrium_concentration

        self.quenched_temperature = quenched_temperature
        self.quenched_ef = quenched_ef
        self.quenched_equilibrium_concentration = \
            quenched_equilibrium_concentration
        self.quenched_carrier_concentrations = quenched_carrier_concentrations

    @classmethod
    def from_calc_results(cls,
                          defect_energies: DefectEnergies,
                          unitcell: UnitcellCalcResults,
                          fractional_magnetization_to_one: bool = False,
                          fractional_criterion: float = 0.1):
        """ Prepare object from DefectEnergies and UnitcellCalcResults objects

        Args:
            defect_energies (DefectEnergies):
                DefectEnergies object used for calculating concentration.
            unitcell (UnitcellCalcResults):
                UnitcellDftResults object for volume and total_dos
            fractional_magnetization_to_one (bool)
                Whether to set the fractional magnetization to 1.
            fractional_criterion (float):
                The criterion to determine if magnetization is fractional.
        """
        energies = deepcopy(defect_energies.defect_energies)

        for name, charge, defect_energy \
                in flatten_dict(defect_energies.defect_energies):
            n = DefectName(name, charge, defect_energy.annotation)
            mag = defect_energy.magnetization

            if abs(mag - round(mag)) > fractional_criterion:
                logger.warning(f"The total_magnetization of {n} "
                               f"is {mag}, and not integer")

                if fractional_magnetization_to_one:
                    logger.warning(f"The magnetization of {n} is set to 1.")
                    energies[name][charge].magnetization = 1.0

        return cls(defect_energies=energies,
                   volume=unitcell.volume,  # [A^3]
                   vbm=defect_energies.vbm,
                   cbm=defect_energies.cbm,
                   total_dos=unitcell.total_dos)

    # @classmethod
    # def from_dict(cls, d):
    #     """ Construct a class object from a dictionary. """
    #     defect_energies = sanitize_keys_in_dict(d["defect_energies"])
    #     equilibrium_concentration = \
    #         sanitize_keys_in_dict(d["equilibrium_concentration"])
    #     quenched_equilibrium_concentration = \
    #         sanitize_keys_in_dict(d["quenched_equilibrium_concentration"])

        # if d["concentrations"] is not None:
        #     concentrations = \
        #         [sanitize_keys_in_dict(d) for d in d["concentrations"]]
        # else:
        #     concentrations = None

        # if d["quenched_carrier_concentrations"] is not None:
        #     quenched_carrier_concentrations = \
        #         [sanitize_keys_in_dict(d)
        #          for d in d["quenched_carrier_concentrations"]]
        # else:
        #     quenched_carrier_concentrations = None

        # return cls(defect_energies=defect_energies,
        #            volume=d["volume"],
        #            vbm=d["vbm"],
        #            cbm=d["cbm"],
        #            total_dos=d["total_dos"],
        #            temperature=d["temperature"],
        #            fermi_mesh=d["fermi_mesh"],
        #            concentrations=concentrations,
        #            equilibrium_ef=d["equilibrium_ef"],
        #            equilibrium_concentration=equilibrium_concentration,
        #            quenched_temperature=d["quenched_temperature"],
        #            quenched_ef=d["quenched_ef"],
        #            quenched_equilibrium_concentration=
        #            quenched_equilibrium_concentration,
        #            quenched_carrier_concentrations=
        #            quenched_carrier_concentrations)

    # def to_json_file(self, filename="defect_concentrations.json"):
    #     with open(filename, 'w') as fw:
    #         json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    # @classmethod
    # def load_json(cls, filename):
    #     return loadfn(filename)

    def __repr__(self):

        outs = [f"volume: {round(self.volume, 2)}",
                f"vbm, cbm: {round(self.vbm, 2)}, {round(self.cbm, 2)}",
                f"band gap: {round(self.cbm - self.vbm, 2)}",
                ""]

        for i, c in enumerate([self.equilibrium_concentration,
                               self.quenched_equilibrium_concentration]):
            if c:
                if i == 0:
                    outs.append("++ Equilibrium concentration")
                    t = self.temperature
                    ef = self.equilibrium_ef - self.vbm
                else:
                    outs.append("++ Quenched equilibrium concentration")
                    t = self.quenched_temperature
                    ef = self.quenched_ef - self.vbm

                p = c["p"][1][None]
                n = c["n"][-1][None]
                out = [f"Temperature: {t} K.",
                       f"Fermi level from vbm: {round(ef, 2)} eV.",
                       f"{'p':>13}: {p:.1e} cm-3.",
                       f"{'n':>13}: {n:.1e} cm-3.",
                       f"{'p - n':>13}: {p - n:.1e} cm-3."]
                outs.extend(out)

                for name in c:
                    if name in ("p", "n"):
                        continue
                    outs.append("---")
                    for charge in c[name]:
                        for annotation, con in c[name][charge].items():
                            n = DefectName(name, charge, annotation)
                            outs.append(f"{str(n):>13}: {con:.1e} cm-3.")
                outs.append("")

        return "\n".join(outs)

    def calc_concentrations(self,
                            temperature: float,
                            fermi_range: list = None,
                            num_mesh: int = 100) -> None:
        """ Calculates defect formation energies from some files.

        Args:
            temperature (list):
                temperature in K.
            fermi_range (list):
                Range of Fermi level, [lower_limit, upper_limit]
            num_mesh (float):
                Number of mesh for Fermi level including boundary.
        """
        fermi_range = fermi_range or [self.vbm, self.cbm]
        self.fermi_mesh = np.linspace(fermi_range[0], fermi_range[1], num_mesh)

        if temperature:
            self.temperature = temperature
        else:
            raise ValueError("Temperature is not defined.")

        self.concentrations = []
        for e_f in self.fermi_mesh:
            self.concentrations.append(
                calc_concentration(defect_energies=self.defect_energies,
                                   temperature=self.temperature,
                                   e_f=e_f,
                                   vbm=self.vbm,
                                   cbm=self.cbm,
                                   total_dos=self.total_dos,
                                   volume=self.volume))

        if self.quenched_temperature:
            self.quenched_carrier_concentrations = []
            for e_f in self.fermi_mesh:
                self.quenched_carrier_concentrations.append(
                    calc_concentration(defect_energies=None,
                                       temperature=self.temperature,
                                       e_f=e_f,
                                       vbm=self.vbm,
                                       cbm=self.cbm,
                                       total_dos=self.total_dos,
                                       volume=self.volume))

    def calc_equilibrium_concentration(self,
                                       temperature: Optional[float] = None,
                                       verbose: bool = True) -> None:
        """
        Calculate equilibrium defect concentrations at given temperature.

        Args:
            temperature (list):
                temperature in K.
            verbose (book):
                Whether print the carrier and defect concentration during the
                seek of the selfconsistent results.
        """
        if temperature:
            self.temperature = temperature
        elif self.temperature is None:
            raise ValueError("Temperature is not defined.")

        self.equilibrium_ef, self.equilibrium_concentration = \
            calc_equilibrium_concentration(defect_energies=self.defect_energies,
                                           temperature=self.temperature,
                                           vbm=self.vbm,
                                           cbm=self.cbm,
                                           total_dos=self.total_dos,
                                           volume=self.volume,
                                           verbose=verbose)

    def calc_quenched_equilibrium_concentration(self,
                                                temperature: float = 298,
                                                verbose: bool = True) -> None:
        """Calculate defect concentrations quenched to low temperature.

        Args:
            temperature (float):
                Temperature in K.
            verbose (book):
                Whether print the carrier and defect concentration during the
                seek of the selfconsistent results.
        """
        if self.equilibrium_concentration is None:
            raise ValueError("To calculate the quenched equilibrium "
                             "concentrations, the equilibrium concentrations "
                             "are needed.")

        self.quenched_temperature = temperature

        self.quenched_ef, self.quenched_equilibrium_concentration = \
            calc_equilibrium_concentration(
                defect_energies=self.defect_energies,
                temperature=self.quenched_temperature,
                vbm=self.vbm,
                cbm=self.cbm,
                total_dos=self.total_dos,
                volume=self.volume,
                verbose=verbose,
                ref_concentration=self.equilibrium_concentration)

    def plot_carrier_concentrations(self,
                                    title: str = None,
                                    xlim: list = None,
                                    ylim: list = None,
                                    set_vbm_zero=True):
        """Get a matplotlib plot.

        Args:
            title (str):
                Title of the plot
            xlim (list):
                Specifies the x-axis limits. None for automatic determination.
            ylim (list):
                Specifies the y-axis limits. None for automatic determination.
            set_vbm_zero (bool):
                Set VBM to zero.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        if self.concentrations is None:
            raise ValueError("The defect and carrier concentrations need to "
                             "be calculated before their plot.")

        title = title or f"Temperature: {self.temperature}  K"
        plt.title(title)

        ax.set_xlabel("Fermi level (eV)")
        ax.set_ylabel("Concentration (cm-3)")
        ax.set_yscale("log", nonposy='clip')

        if xlim:
            plt.xlim(xlim[0], xlim[1])

        if ylim:
            ylim = [10 ** ylim[0], 10 ** ylim[1]]
        else:
            max_y = max([con for concentration in self.concentrations
                         for d in ("p", "n")
                         for c in concentration[d]
                         for con in concentration[d][c].values()])
            ylim = [10 ** 10, max_y * 2]

        plt.ylim(ylim[0], ylim[1])

        if set_vbm_zero is True:
            vbm = 0.0
            cbm = self.cbm - self.vbm
            fermi_mesh = [f - self.vbm for f in self.fermi_mesh]
            equilib_ef = self.equilibrium_ef - self.vbm
            quenched_ef = self.quenched_ef - self.vbm
        else:
            vbm = self.vbm
            cbm = self.cbm
            fermi_mesh = self.fermi_mesh
            equilib_ef = self.equilibrium_ef
            quenched_ef = self.quenched_ef

        plt.axvline(x=vbm, linewidth=1.0, linestyle='dashed')
        plt.axvline(x=cbm, linewidth=1.0, linestyle='dashed')

        hole1 = [v["p"][1][None] for v in self.concentrations]
        elec1 = [v["n"][-1][None] for v in self.concentrations]
        hole2 = [v["p"][1][None] for v in self.quenched_carrier_concentrations]
        elec2 = [v["n"][-1][None] for v in self.quenched_carrier_concentrations]

        ax.plot(fermi_mesh, hole1, '-', color="blue", label="p")
        ax.plot(fermi_mesh, elec1, '-', color="red", label="n")
        ax.plot(fermi_mesh, hole2, '--', color="blue", label="p")
        ax.plot(fermi_mesh, elec2, '--', color="red", label="n")

        plt.axvline(x=equilib_ef, linewidth=1.0, linestyle=':', color='g')
        plt.axvline(x=quenched_ef, linewidth=1.0, linestyle=':', color='g')

        ax.annotate(f"{round(self.temperature, 1)} K",
                    (equilib_ef, ylim[1]), fontsize=10, ha='center', color='g')
        ax.annotate(f"{round(self.quenched_temperature, 1)} K",
                    (quenched_ef, ylim[1]), fontsize=10, ha='center', color='g')

        # plt.arrow(x=equilibrium_ef, y=ylim[1] / 2,
        #           dx=quenched_ef - equilibrium_ef, dy=0,
        #           width=0.005,
        #           head_width=0.5,
        #           head_length=0.1,
        #           length_includes_head=True,
        #           color='green')

        return plt
