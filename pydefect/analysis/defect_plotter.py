# -*- coding: utf-8 -*-

from copy import deepcopy
from enum import Enum
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from monty.serialization import loadfn
import json
import warnings

from monty.json import MontyEncoder
from monty.serialization import loadfn

from collections import defaultdict, namedtuple
from itertools import combinations

from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.analysis.defect_concentration import DefectConcentration
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.input_maker.defect_set_maker import is_name_selected
from pydefect.util.carrier_concentration import maxwell_boltzmann_dist, \
    CarrierConcentration


class DefectPlotter:
    """
    """
    def __init__(self, defect_energies, defect_concentration=None):
        """
        Calculates defect formation energies.
        Args:
            defect_energies (DefectEnergies):
            defect_concentration (DefectConcentration):
        """
        # Objects of the Concentration named tuple for 1st and 2nd temperature.
        self._title = defect_energies.title
        self._energies = defect_energies.energies
        self._transition_levels = defect_energies.transition_levels
        self._vbm = defect_energies.vbm
        self._cbm = defect_energies.cbm
        self._supercell_vbm = defect_energies.supercell_vbm
        self._supercell_cbm = defect_energies.supercell_cbm
        self._band_gap = self._cbm - self._vbm

        self._e_f1 = None
        self._e_f2 = None
        self._t1 = None
        self._t2 = None

        if defect_concentration:
            self._e_f1 = defect_concentration.e_f
            self._t1 = defect_concentration.temperature

        if defect_concentration.previous_concentration:
            self._e_f2 = defect_concentration.previous_concentration.e_f
            self._t2 = defect_concentration.previous_concentration.temperature

    def plot_energy(self, filtering_words=None, x_range=None,
                    y_range=None, show_Fermi_level=True,
                    show_transition_levels=False):
        """
        Plots the defect formation energies as a function of the Fermi level.

        Args:
            filtering_words (str):
                Words used for filtering the defect types
            x_range (2x1 list):
                x range for the plot.
            y_range (2x1 list):
                y range for the plot.
            show_transition_levels (bool):
                Whether the transition levels are shown in the plot.
        """

        fig, ax = plt.subplots()

        plt.title(self._title, fontsize=15)

        ax.set_xlabel("Fermi level (eV)", fontsize=15)
        ax.set_ylabel("Formation energy (eV)", fontsize=15)

        # filtering specific words for plot.
        transition_levels = deepcopy(self._transition_levels)
        if filtering_words:
            for name in self._transition_levels.keys():
                if is_name_selected(name, filtering_words) is False:
                    transition_levels.pop(name)

        if transition_levels == {}:
            raise KeyError("No transition levels are detected.")

        if x_range:
            x_min = x_range[0]
            x_max = x_range[1]
        else:
            x_min = 0
            x_max = self._band_gap

        if y_range:
            y_min = y_range[0]
            y_max = y_range[1]
        else:
            y_min = min([min(np.array(tl.cross_points).transpose()[1])
                         for tl in transition_levels.values()])
            y_max = max([max(np.array(tl.cross_points).transpose()[1])
                         for tl in transition_levels.values()])
            if y_min > 0:
                y_min = 0
            # Make top and bottom space
            margin = (y_max - y_min) * 0.1
            y_min -= margin
            y_max += margin

        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

        # plot vbm, cbm, supercell_vbm, supercell_cbm lines.
        if x_range:
            plt.axvline(x=0, linewidth=1.0, linestyle='dashed')
            plt.axvline(x=self._band_gap, linewidth=1.0, linestyle='dashed')
            ax.annotate("cbm", (self._band_gap, y_min * 0.9 + y_max * 0.1),
                        fontsize=10)

            supercell_vbm = self._supercell_vbm - self._vbm
            supercell_cbm = self._supercell_cbm - self._vbm
            if supercell_vbm < - 0.05:
                plt.axvline(x=supercell_vbm, linewidth=1.0, linestyle='dotted')
                ax.annotate("supercell vbm",
                            (supercell_vbm, y_min * 0.9 + y_max * 0.1),
                            fontsize=10)
            if supercell_cbm > self._band_gap + 0.05:
                plt.axvline(x=supercell_cbm, linewidth=1.0, linestyle='dotted')
                ax.annotate("supercell cbm",
                            (supercell_cbm, y_min * 0.9 + y_max * 0.1),
                            fontsize=10)

        plt.axhline(y=0, linewidth=1.0, linestyle='dashed')

        # Lines for the equilibrium Fermi level
        if show_Fermi_level:
            if self._e_f1:
                e_f1 = self._e_f1 - self._vbm
                plt.axvline(x=e_f1, linewidth=2.0, linestyle='dashed',
                            color="red")
                ax.annotate("T$_1$=" + str(self._t1) + "K", (e_f1, y_min),
                            fontsize=7)
            if self._e_f2:
                e_f2 = self._e_f2 - self._vbm
                plt.axvline(x=e_f2, linewidth=2.0, linestyle='dashed',
                            color="purple")
                ax.annotate("T$_2$=" + str(self._t2) + "K", (e_f2, y_min),
                            fontsize=7)

        for i, (name, tl) in enumerate(transition_levels.items()):

            color = matplotlib.cm.hot(float(i) / len(self._transition_levels))

            # ---------- Calculate cross points --------------------------------
            energies = self._energies[name]
            cross_points = []
            shallow_points = []
            transition_points = []
            charge_set = set()

            # x_min
            c, y = DefectEnergies.min_e_at_ef(energies, x_min + self._vbm)
            cross_points.append([x_min, y])
            # Adding charge is necessary when transition level is absent.
            charge_set.add(c)
            if c < 0:
                shallow_points.append([x_min, y])

            for cp, c in zip(tl.cross_points, tl.charges):
                if x_min < cp[0] - self._vbm < x_max:
                    cross_points.append([cp[0] - self._vbm, cp[1]])
                    transition_points.append([cp[0] - self._vbm, cp[1]])
                    charge_set.add(c[0])
                    charge_set.add(c[1])
            # x_max
            c, y = DefectEnergies.min_e_at_ef(energies, x_max + self._vbm)
            cross_points.append([x_max, y])
            if c > 0:
                shallow_points.append([x_max, y])
            # ------------------------------------------------------------------

            # set the x and y arrays to be compatible with matplotlib style.
            # e.g., x = [0.0, 0.3, 0.5, ...], y = [2.1, 3.2, 1.2, ...]
            x, y = np.array(cross_points).transpose()
            line, = ax.plot(x, y, '-', color=color, label=name)
            line.set_label(name)

            # plot filled circles for transition levels.
            if transition_points:
                x, y = np.array(transition_points).transpose()
                ax.scatter(x, y, color=color)

            # plot unfilled circles for shallow levels.
            if shallow_points:
                x, y = np.array(shallow_points).transpose()
                ax.scatter(x, y, facecolor="white", edgecolor=color)

            # These determine the positions of the transition levels.
            margin_x = (x_max - x_min) * 0.025
            margin_y = (y_max - y_min) * 0.025
            if show_transition_levels:
                for cp in cross_points:
                    if cp[0] < x_min + 0.0001 or cp[0] > x_max - 0.0001:
                        continue
                    s = str(round(cp[0], 2)) + ", " + str(round(cp[1], 2))
                    ax.annotate(s, (cp[0] + margin_x, cp[1] - margin_y),
                                color=color, fontsize=8)

            # Arrange the charge states at the middle of the transition levels.
            charge_state_position = \
                [[(a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + margin_y]
                 for a, b in zip(cross_points, cross_points[1:])]
            for j, (x, y) in enumerate(charge_state_position):
                ax.annotate(str(sorted(charge_set)[j]), (x, y), color=color,
                            fontsize=13)

        ax.legend()
        fig.subplots_adjust(right=0.75)

        return plt

# TODO: draw thin lines for options

