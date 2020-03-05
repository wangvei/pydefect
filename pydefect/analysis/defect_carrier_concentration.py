# -*- coding: utf-8 -*-

from collections import defaultdict
from copy import deepcopy
from typing import Optional, Tuple, Dict
import json

from matplotlib import pyplot as plt

from monty.serialization import loadfn
from monty.json import MontyEncoder, MSONable

import numpy as np

from pydefect.core.defect_name import DefectName
from pydefect.core.error_classes import NoConvergenceError
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.distribution_function import (
    maxwell_boltzmann_distribution, fermi_dirac_distribution)
from pydefect.util.logger import get_logger
from pydefect.util.tools import (
    flatten_dict, sanitize_keys_in_dict)

logger = get_logger(__name__)


def hole_concentration(temperature: float,
                       e_f: float,
                       total_dos: list,
                       vbm: float,
                       volume: float,
                       threshold: float = 0.05) -> float:
    """Hole carrier concentration at the given absolute fermi_level.

    Args:
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        total_dos (list):
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
                           total_dos: list,
                           cbm: float,
                           volume: float,
                           threshold: float = 0.05) -> float:
    """Electron carrier concentration at the given absolute fermi_level.

    Args:
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        total_dos (list):
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


def calc_concentration(defect_energies: Optional[Dict],
                       temperature: float,
                       e_f: float,
                       vbm: float,
                       cbm: float,
                       total_dos: list,
                       volume: float,
                       ref_concentration: Optional[dict] = None) -> dict:
    """Calculate defect/carrier concentrations at a temperature & a Fermi level.

    When the reference_concentration is provided, each defect concentration is
    fixed. For example, n(Va_Mg) = n(Va_Mg_0) + n(Va_Mg_-1) + n(Va_Mg_-2) is
    kept fixed even the temperature is quenched.

    Args:
        defect_energies (dict):
            Defect formation energies.
            defect_energies[name][charge] = DefectEnergy object.
        temperature (float):
            Temperature in K.
        e_f (float):
            Fermi level in the absolute scale.
        vbm (float):
            Valence band maximum in the unitcell in the absolute scale.
        cbm (float):
            Conduction band minimum in the unitcell in the absolute scale.
        total_dos (list):
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
    concentrations = defaultdict(dict)

    concentrations["p"][1] = \
        hole_concentration(temperature, e_f, total_dos, vbm, volume)
    concentrations["n"][-1] = \
        electron_concentration(temperature, e_f, total_dos, cbm, volume)

    if defect_energies is None:
        return dict(concentrations)

    for name in defect_energies:
        concentration_by_name = {}
        for charge, de in defect_energies[name].items():

            num_mag_conf = abs(de.magnetization) + 1
            degree_of_freedom = de.multiplicity * num_mag_conf

            energy = de.defect_energy + e_f * charge
            # volume unit conversion from [A^3] to [cm^3]
            concentration_by_name[charge] = \
                (maxwell_boltzmann_distribution(energy, temperature)
                 * degree_of_freedom / (volume / 10 ** 24))

        if ref_concentration:
            ref_total_concentration = sum(ref_concentration[name].values())
            total_concentration = sum(concentration_by_name.values())
            factor = (ref_total_concentration / total_concentration)

            concentration_by_name = \
                {k: v * factor for k, v in concentration_by_name.items()}

        concentrations[name] = concentration_by_name

    return dict(concentrations)


def calc_equilibrium_concentration(defect_energies: dict,
                                   temperature: float,
                                   vbm: float,
                                   cbm: float,
                                   total_dos: list,
                                   volume: float,
                                   ref_concentration: Optional[dict] = None,
                                   verbose: bool = False,
                                   max_iteration: int = 200,
                                   interval_decay_parameter: float = 0.7,
                                   threshold: float = 1e-5
                                   ) -> Tuple[float, dict]:
    """Calculates equilibrium carrier & defect concentration at a temperature.

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
        total_dos (list):
            Total density of states, [[dos1, dos2, ..], [energy1, energy2, ..]]
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
                 for charge, energy in defect_concentration[name].items()])

        if verbose:
            logger.info(f"- {iteration}th iteration ----")
            logger.info(f"Fermi level: {e_f:.2f} eV.")

            for name in defect_concentration:
                for charge in defect_concentration[name]:
                    concentration = defect_concentration[name][charge]
                    logger.info(f"{name:>8}  {charge:2d}:"
                                f"   {concentration:.1e} cm-3.")

            logger.info(f"Charge sum: {total_charge:.1e} cm-3.")

        interval *= interval_decay_parameter
        e_f = e_f + np.sign(total_charge) * interval

        max_concentration = \
            np.amax([v[-1] for v in flatten_dict(defect_concentration)])

        # This part controls the accuracy.
        if np.abs(total_charge / max_concentration) < threshold:
            return e_f, defect_concentration

    raise NoConvergenceError("Equilibrium condition has not been reached.")


class DefectConcentration(MSONable):
    """Class related to a set of carrier and defect concentration."""

    def __init__(self,
                 defect_energies: dict,
                 volume: float,
                 vbm: float,
                 cbm: float,
                 total_dos: list,
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
                energies[name][charge] = DefectEnergy object
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
                A set of carrier and defect concentrations in cm-3 as a
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
                          defect_energies: dict,
                          unitcell: UnitcellCalcResults,
                          round_magnetization: bool = False,
                          fractional_criterion: float = 0.1
                          ) -> "DefectConcentration":
        """Create instance object from DefectEnergies and UnitcellCalcResults.

        Args:
            defect_energies (dict):
                DefectEnergy as a function of name, charge, and annotation.
                energies[name][charge] = DefectEnergy object
            unitcell (UnitcellCalcResults):
                UnitcellDftResults object for volume, band edges and total_dos
            round_magnetization (bool)
                Whether to round the fractional magnetization.
            fractional_criterion (float):
                The criterion to determine if the magnetization is fractional.

        Returns:
            DefectConcentration object.
        """
        defect_energies = deepcopy(defect_energies)

        for name, charge, defect_energy in flatten_dict(defect_energies):
            defect_name = DefectName(name, charge, defect_energy.annotation)
            mag = defect_energy.magnetization
            rounded_mag = round(mag)
            if abs(mag - rounded_mag) > fractional_criterion:
                logger.warning(f"The total_magnetization of {defect_name} "
                               f"is {mag}, and not integer")

                if round_magnetization:
                    logger.warning(f"The magnetization of {defect_name} is "
                                   f"rounded to {rounded_mag}.")
                    defect_energies[name][charge].magnetization = rounded_mag

        for attr in ["volume", "band_edge", "total_dos"]:
            if unitcell.__getattribute__(attr) is None:
                logger.critical(f"{attr} needs to be set for the calculation of"
                                f" carrier/defect concentrations.")
                raise ValueError

        return cls(defect_energies=defect_energies,
                   volume=unitcell.volume,  # [A^3]
                   vbm=unitcell.band_edge[0],
                   cbm=unitcell.band_edge[1],
                   total_dos=unitcell.total_dos)

    @classmethod
    def from_dict(cls, d):
        """Construct a class object from a dictionary. """
        d["defect_energies"] = sanitize_keys_in_dict(d["defect_energies"])
        d["equilibrium_concentration"] = \
            sanitize_keys_in_dict(d["equilibrium_concentration"])
        d["quenched_equilibrium_concentration"] = \
            sanitize_keys_in_dict(d["quenched_equilibrium_concentration"])

        d.pop("@module", None)
        d.pop("@class", None)
        d.pop("@version", None)

        return cls(**d)

    def to_json_file(self, filename: str = "defect_concentrations.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename: str) -> "DefectConcentration":
        return loadfn(filename)

    def __repr__(self):

        outs = [f"volume (A^3): {round(self.volume, 2)}",
                f"vbm, cbm (eV): {round(self.vbm, 2)}, {round(self.cbm, 2)}",
                f"band gap (eV): {round(self.cbm - self.vbm, 2)}", ""]

        for equil_c, quenched_c in enumerate([self.equilibrium_concentration,
                                     self.quenched_equilibrium_concentration]):
            if quenched_c:
                if equil_c == 0:
                    outs.append("++ Equilibrium concentration")
                    t = self.temperature
                    ef = self.equilibrium_ef - self.vbm
                else:
                    outs.append("++ Quenched equilibrium concentration")
                    t = self.quenched_temperature
                    ef = self.quenched_ef - self.vbm

                p = quenched_c["p"][1]
                n = quenched_c["n"][-1]
                out = [f"Temperature: {t} K.",
                       f"Fermi level from vbm: {round(ef, 2)} eV.",
                       f"{'p':>13}: {p:.1e} cm-3.",
                       f"{'n':>13}: {n:.1e} cm-3.",
                       f"{'p - n':>13}: {p - n:.1e} cm-3."]
                outs.extend(out)

                for name in quenched_c:
                    if name in ("p", "n"):
                        continue
                    outs.append("---")
                    for charge, concentration in quenched_c[name].items():
                        annotation = \
                            self.defect_energies[name][charge].annotation
                        n = DefectName(name, charge, annotation)
                        outs.append(f"{str(n):>13}: {concentration:.1e} cm-3.")
                outs.append("")

        return "\n".join(outs)

    def calc_concentrations(self,
                            temperature: float,
                            fermi_range: list = None,
                            num_mesh: int = 100) -> None:
        """Calculate carrier & defect concentrations as function of Fermi level.

        Quenched concentrations are also calculated by default.

        Args:
            temperature (list):
                temperature in K.
            fermi_range (list):
                Range of Fermi level, [lower_limit, upper_limit]
            num_mesh (float):
                Number of mesh for Fermi level including boundary.

        Returns:
            None
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
        """Calculate equilibrium defect concentrations at a given temperature.

        Args:
            temperature (list):
                temperature in K.
            verbose (book):
                Whether print the carrier and defect concentration during the
                seek of the selfconsistent results.

        Returns:
            None
        """
        if temperature is not None:
            self.temperature = temperature
        else:
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

        Returns:
            None
        """
        if self.equilibrium_concentration is None:
            raise ValueError(
                "To calculate the quenched equilibrium concentrations, the "
                "equilibrium concentrations are needed.")

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
                                    set_vbm_zero=True) -> plt:
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

        Returns:
            Matplotlib pyplot.
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
                         for c, con in concentration[d].items()])
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
