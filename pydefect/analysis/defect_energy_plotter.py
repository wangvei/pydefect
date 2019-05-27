# -*- coding: utf-8 -*-

from collections import namedtuple
from itertools import combinations
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.analysis.defect_carrier_concentration import DefectConcentration
from pydefect.input_maker.defect_initial_setting import SimpleDefectName

TransitionLevel = namedtuple("TransitionLevel", ("cross_points", "charges"))


class DefectEnergyPlotter:
    """ Plotter of defect formation energies as a function of Fermi level. """

    def __init__(self,
                 defect_energies: DefectEnergies):
#                 defect_concentration: DefectConcentration = None):
        """ Store data used for the plot.

        Args:
            defect_energies (DefectEnergies):
                DefectEnergies class object.
            # defect_concentration (DefectConcentration):
            #     DefectConcentration class object.
        """
        # Objects of the Concentration named tuple for 1st and 2nd temperature.
        self.title = defect_energies.title
        self.energies = defect_energies.energies
        self.vbm = defect_energies.vbm
        self.cbm = defect_energies.cbm
        self.supercell_vbm = defect_energies.supercell_vbm
        self.supercell_cbm = defect_energies.supercell_cbm
        self.band_gap = self.cbm - self.vbm

        self.transition_levels = {}
        self.e_f1 = None
        self.e_f2 = None
        self.t1 = None
        self.t2 = None

        # if defect_concentration:
        #     self.e_f1 = defect_concentration.e_f
        #     self.t1 = defect_concentration.temperature

            # if defect_concentration.previous_concentration:
            #     self.e_f2 = defect_concentration.previous_concentration.e_f
            #     self.t2 = \
            #         defect_concentration.previous_concentration.temperature

    def plot_energy(self,
                    filtering_words: list = None,
                    x_range: list = None,
                    y_range: list = None,
                    show_fermi_level: bool = True,
                    exclude_shallow_defects: bool = True,
                    show_transition_levels: bool = False,
                    show_all_lines: bool = False):
        """ Plots defect formation energies as a function of the Fermi level.

        Args:
            filtering_words (list):
                List of words used for filtering the defects
            x_range (list):
                1x2 list for determining x range.
            y_range (list):
                1x2 list for determining y range.
            show_fermi_level (bool):
                Whether to show the Fermi level(s) in the plot.
            show_transition_levels (bool):
                Whether to show the transition levels in the plot.
            show_all_lines (bool):
                Whether to show all lines in the plot.
        """
        fig, ax = plt.subplots()
        plt.title(self.title, fontsize=15)

        ax.set_xlabel("Fermi level (eV)", fontsize=15)
        ax.set_ylabel("Formation energy (eV)", fontsize=15)

        # e_of_c is energy as a function of charge: e_of_c[charge] = energy
        for name, e_of_c in self.energies.items():

            for c in e_of_c.keys():
                n = SimpleDefectName.from_str("_".join([name, str(c)]))
                if filtering_words \
                        and n.is_name_selected(filtering_words) is False:
                    e_of_c.pop(c)

            cross_points_charge = []

            for (c1, e1), (c2, e2) in combinations(e_of_c.items(), r=2):
                # The cross point between two charge states.
                x = - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                # The lowest energy among all the charge states to be compared.
                compared_energy = \
                    min([energy + c * x for c, energy in e_of_c.items()])

                if y < compared_energy + 1e-5:
                    cross_points_charge.append([[x, y], [c1, c2]])

            # need to sort the points along x-axis.
            cross_points = []
            charge = []
            for cp, c in sorted(cross_points_charge, key=lambda z: z[0][0]):
                cross_points.append(cp)
                charge.append(c)

            self.transition_levels[name] = \
                TransitionLevel(cross_points=cross_points, charges=charge)

        if x_range:
            x_min = x_range[0]
            x_max = x_range[1]
        else:
            x_min = 0
            x_max = self.band_gap

        # Lines for the equilibrium Fermi level
        # if show_fermi_level:
        #     if self.e_f1:
        #         e_f1 = self.e_f1 - self.vbm
        #         plt.axvline(x=e_f1, linewidth=2.0, linestyle='dashed',
        #                     color="red")
        #         ax.annotate("T$_1$=" + str(self.t1) + "K", (e_f1, y_min),
        #                     fontsize=7)
        #     if self.e_f2:
        #         e_f2 = self.e_f2 - self.vbm
        #         plt.axvline(x=e_f2, linewidth=2.0, linestyle='dashed',
        #                     color="purple")
        #         ax.annotate("T$_2$=" + str(self.t2) + "K", (e_f2, y_min),
        #                     fontsize=7)

        y_min = float("inf")
        y_max = -float("inf")

        for i, (name, tl) in enumerate(self.transition_levels.items()):

            color = matplotlib.cm.hot(float(i) / len(self.transition_levels))

            # ---------- Calculate cross points including edges ----------------
            energies = self.energies[name]
            cross_points = []
            shallow_points = []
            transition_points = []
            charge_set = set()

            # must keep sequence (x_min -> transition levels -> x_max) for plot
            # x_min
            c, y = min_e_at_ef(energies, x_min + self.vbm)
            cross_points.append([x_min, y])
            # Adding charge is necessary when transition level is absent.
            charge_set.add(c)
            if c < 0:
                shallow_points.append([x_min, y])
            if y < y_min:
                y_min = y
            if y > y_max:
                y_max = y

            for cp, c in zip(tl.cross_points, tl.charges):
                if x_min < cp[0] - self.vbm < x_max:
                    cross_points.append([cp[0] - self.vbm, cp[1]])
                    transition_points.append([cp[0] - self.vbm, cp[1]])
                    charge_set.add(c[0])
                    charge_set.add(c[1])
                    if cp[1] < y_min:
                        y_min = cp[1]
                    if cp[1] > y_max:
                        y_max = cp[1]

            # x_max
            c, y = min_e_at_ef(energies, x_max + self.vbm)
            cross_points.append([x_max, y])
            if c > 0:
                shallow_points.append([x_max, y])
            if y < y_min:
                y_min = y
            if y > y_max:
                y_max = y

            # ------------------------------------------------------------------

            # set x and y arrays to be compatible with matplotlib style.
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

            # margin_x and _y determine the positions of the transition levels.
            margin_x = (x_max - x_min) * 0.025
            margin_y = (y_max - y_min) * 0.025

            if show_transition_levels:
                for cp in cross_points:
                    if x_min + 1e-4 < cp[0] < x_max - 1e-4 is False:
                        continue
                    s = str(round(cp[0], 2)) + ", " + str(round(cp[1], 2))
                    pos_x = cp[0] + margin_x
                    pos_y = cp[1] - 1.8 * margin_y
                    ax.annotate(s, (pos_x, pos_y), color=color, fontsize=10)

            # Arrange the charge states at the middle of the transition levels.
            charge_pos = [[(a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + margin_y]
                          for a, b in zip(cross_points, cross_points[1:])]

            sorted_charge_set = sorted(charge_set, reverse=True)

            for j, (x, y) in enumerate(charge_pos):
                s = str(sorted_charge_set[j])
                ax.annotate(s, (x, y), color=color, fontsize=13)

            if show_all_lines:
                for c, e in self.energies[name].items():
                    y1 = e + c * (x_min + self.vbm)
                    y2 = e + c * (x_max + self.vbm)
                    ax.plot([x_min, x_max], [y1, y2], '-', linewidth=0.3,
                            color=color)

        if y_range:
            y_min = y_range[0]
            y_max = y_range[1]
        else:
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
            plt.axvline(x=self.band_gap, linewidth=1.0, linestyle='dashed')

            ax.annotate("vbm", (0, y_min * 0.9 + y_max * 0.1),
                        fontsize=10)
            ax.annotate("cbm", (self.band_gap, y_min * 0.9 + y_max * 0.1),
                        fontsize=10)

            supercell_vbm = self.supercell_vbm - self.vbm
            supercell_cbm = self.supercell_cbm - self.vbm

            if supercell_vbm < - 0.05:
                plt.axvline(x=supercell_vbm, linewidth=1.0, linestyle='dotted')
                ax.annotate("supercell vbm",
                            (supercell_vbm, y_min * 0.9 + y_max * 0.1),
                            fontsize=10)

            if supercell_cbm > self.band_gap + 0.05:
                plt.axvline(x=supercell_cbm, linewidth=1.0, linestyle='dotted')
                ax.annotate("supercell cbm",
                            (supercell_cbm, y_min * 0.9 + y_max * 0.1),
                            fontsize=10)

        plt.axhline(y=0, linewidth=1.0, linestyle='dashed')
        ax.legend()
#        fig.subplots_adjust(right=0.75)

        return plt


def min_e_at_ef(ec, ef):
    # calculate each energy at the given Fermi level ef.
    d = {c: e + c * ef for c, e in ec.items()}
    # return the charge with the lowest energy, and its energy value
    return min(d.items(), key=lambda x: x[1])
