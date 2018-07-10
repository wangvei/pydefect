# -*- coding: utf-8 -*-

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

    def plot_energy(self, file_name=None, filtering_words=None, x_range=None,
                    y_range=None, show_transition_levels=False):
        """
        Plots the defect formation energies as a function of the Fermi level.

        Args:
            file_name (str):
                File name for saving the plot.
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

        if self._transition_levels is None:
            raise NoTLError("Transition levels are not calculated yet.")

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
                         for tl in self._transition_levels.values()])
            y_max = max([max(np.array(tl.cross_points).transpose()[1])
                         for tl in self._transition_levels.values()])
            if y_min > 0:
                y_min = 0
            # Make top and bottom space
            margin = (y_max - y_min) * 0.08
            y_min -= margin
            y_max += margin

        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

        # support lines
        if x_range:
            plt.axvline(x=0, linewidth=1.0, linestyle='dashed')
            plt.axvline(x=self._band_gap, linewidth=1.0, linestyle='dashed')
            ax.annotate("cbm", (self._band_gap, y_min * 0.9 + y_max * 0.1),
                        fontsize=10)

            perfect_vbm = self._supercell_vbm - self._vbm
            perfect_cbm = self._supercell_cbm - self._vbm
            print(self._supercell_vbm)
            print(perfect_vbm)
            print(perfect_cbm)
            if perfect_vbm < - 0.05:
                plt.axvline(x=perfect_vbm, linewidth=1.0, linestyle='dotted')
                ax.annotate("supercell vbm",
                            (perfect_vbm, y_min * 0.9 + y_max * 0.1),
                            fontsize=10)
            if perfect_cbm > self._band_gap + 0.05:
                plt.axvline(x=perfect_cbm, linewidth=1.0, linestyle='dotted')
                ax.annotate("supercell cbm",
                            (perfect_cbm, y_min * 0.9 + y_max * 0.1),
                            fontsize=10)

        plt.axhline(y=0, linewidth=1.0, linestyle='dashed')

        # Lines for the equilibrium Fermi level(s)
        if self._e_f1:
            plt.axvline(x=self._e_f1, linewidth=2.0, linestyle='dashed', color="red")
            ax.annotate("T$_1$=" + str(self._t1) + "K", (self._e_f1, y_min),
                        fontsize=7)
            if self._e_f2:
                plt.axvline(x=self._e_f2, linewidth=2.0, linestyle='dashed',
                            color="purple")
                ax.annotate("T$_2$=" + str(self._t2) + "K",
                            (self._e_f2, self._t2), fontsize=7)

        for i, (name, tl) in enumerate(self._transition_levels.items()):

            color = matplotlib.cm.hot(float(i) / len(self._transition_levels))

            cross_points = tl.cross_points
            shallow = []
            transition_levels = []

            # Store the coords for the filled and open circles.
            for cp in cross_points:
                # tl_x_min = self._calculated_transition_level_range[0]
                # # tl_x_max = self._calculated_transition_level_range[1]
                # if (cp[0] < tl_x_min + 0.0001 and max(tl.charges) < 0) \
                #         or (cp[0] > tl_x_max - 0.0001 and min(tl.charges) > 0):
                #     shallow.append(cp)
                # elif not (cp[0] == tl_x_min or cp[0] == tl_x_max):
                #     transition_levels.append(cp)
                if (cp[0] < x_min + 0.0001 and max(tl.charges) < 0) \
                        or (cp[0] > x_max - 0.0001 and min(tl.charges) > 0):
                    shallow.append(cp)
                elif not (cp[0] == x_min or cp[0] == x_max):
                    transition_levels.append(cp)

            # set the x and y arrays to be compatible with matplotlib style.
            # e.g., x = [0.0, 0.3, 0.5, ...], y = [2.1, 3.2, 1.2, ...]
            # plot lines
            x, y = np.array(cross_points).transpose()
            line, = ax.plot(x, y, '-', color=color, label=name)
            line.set_label(name)

            # plot filled circles for transition levels.
            if transition_levels:
                x, y = np.array(transition_levels).transpose()
                ax.scatter(x, y, color=color)

            # plot unfilled circles for shallow levels.
            if shallow:
                x, y = np.array(shallow).transpose()
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
            middle_points = \
                reversed([[(a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + margin_y]
                          for a, b in zip(cross_points, cross_points[1:])])

            for j, (x, y) in enumerate(middle_points):
                ax.annotate(str(tl.charges[j]), (x, y), color=color,
                            fontsize=13)

        ax.legend()
        fig.subplots_adjust(right=0.75)

        if file_name:
            plt.savefig(file_name)
        else:
            plt.show()


class NoTLError(Exception):
    pass

