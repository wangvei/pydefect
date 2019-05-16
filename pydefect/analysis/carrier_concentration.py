# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from pydefect.util.distribution_function import fermi_dirac_distribution

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class CarrierConcentration:
    """
    Carrier calc_defect_carrier_concentration
    """

    def __init__(self,
                 temperatures: list, vbm: float, cbm: float, fermi_levels: list,
                 ns: list, ps: list):
        """
        Args:
            temperatures (list): list of temperature considered.
            fermi_levels: Fermi levels considered in the absolute scale in eV.
            ns[temperature][Fermi level]: Carrier electron concentrations.
            ps[temperature][Fermi level]: Carrier hole concentrations.
        """
        self.temperatures = temperatures
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
    def from_unitcell(cls, temperatures, unitcell, e_range=None):
        """
        Calculates defect formation energies.
        Args:
            temperatures (list):
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

            ns.append(cls.n(temperatures, f, total_dos, cbm, volume))
            ps.append(cls.p(temperatures, f, total_dos, vbm, volume))

        return cls(temperatures, vbm, cbm, fermi_levels, ns, ps)

    @staticmethod
    def n(temperature, fermi_level, total_dos, cbm, volume, threshold=0.05):
        """ electron carrier concentration at the given absolute fermi_level.
        """
        mesh_distance = total_dos[1][1] - total_dos[1][0]
        n = sum(fermi_dirac_distribution(e, fermi_level, temperature) * td
                for td, e in zip(total_dos[0], total_dos[1])
                if e >= cbm - threshold)
        return n * mesh_distance / volume

    @staticmethod
    def p(temperature, fermi_level, total_dos, vbm, volume, threshold=0.05):
        """ hole carrier concentration at the given absolute fermi_level.
        """
        mesh_distance = total_dos[1][1] - total_dos[1][0]
        p = sum(fermi_dirac_distribution(fermi_level, e, temperature) * td
                for td, e in zip(total_dos[0], total_dos[1])
                if e <= vbm + threshold)
        return p * mesh_distance / volume

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