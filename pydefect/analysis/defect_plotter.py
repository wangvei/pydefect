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
                tl_x_min = self._calculated_transition_level_range[0]
                tl_x_max = self._calculated_transition_level_range[1]
                if (cp[0] < tl_x_min + 0.0001 and max(tl.charges) < 0) \
                        or (cp[0] > tl_x_max - 0.0001 and min(tl.charges) > 0):
                    shallow.append(cp)
                elif not (cp[0] == tl_x_min or cp[0] == tl_x_max):
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

# class DefectAll:
#     """
#     A class related to a set of defect formation energies.
#     """
#     def __init__(self, unitcell, perfect, defects, chem_pot, chem_pot_label,
#                  filtering_words=None, system_name=""):
#         """
#         Calculates defect formation energies.
#         Args:
#             unitcell (UnitcellDftResults):
#             perfect (SupercellDftResults)
#             defects (list of namedtuple Defect):
#                 [[Defect, ...]
#                 Defect = namedtuple("Defect", "defect_entry", "dft_results",
#                                     "correction")
#             chem_pot (ChemPot): Chempot class object.
#             chem_pot_label (str):
#             filtering_words (str): Words used for filtering the defect types
#             system_name (str): A system name used for the title
#         """
#         # Objects of the Concentration named tuple for 1st and 2nd temperature.
#         self._original_concent = None
#         self._quenched_concent = None
#         # Object of the TransitionLevel named tuple.
#         self._transition_levels = None
#         self._calculated_transition_level_range = None

        # # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy.
        # self._vbm, self._cbm = unitcell.band_edge
        # # The band edges in the perfect supercell calculation
        # self._perfect_cbm, self._perfect_vbm = \
        #     perfect.eigenvalue_properties[1:3]
        # self._vbm2, self._cbm2 = (None, None)

        # self._volume = unitcell.volume * 10 ** -24  # [A^3] -> [cm^3]
        # self._total_dos = unitcell.total_dos
        # self._title = system_name + " condition " + chem_pot_label
        # self._defects = defects

        # # chemical potentials
        # _, standard_e = chem_pot
        # relative_chem_pot = _[chem_pot_label]

        # # defect formation energies at the vbm
        # self._defect_energies = defaultdict(dict)
        # self._num_defect_sites = {}

        # for d in defects:
        #     name = d.defect_entry.name
        #     if filtering_words \
        #             and is_name_selected(name, filtering_words) is False:
        #         continue
        #     charge = d.defect_entry.charge
        #     element_diff = d.defect_entry.element_diff
        #     # calculate four terms for a defect formation energy.
        #     relative_energy = d.dft_results.relative_total_energy(perfect)
        #     correction_energy = d.correction.total_correction_energy
        #     electron_interchange_energy = self._vbm * charge
        #     element_interchange_energy = \
        #         - sum([v * (relative_chem_pot.elem_coords[k] + standard_e[k])
        #                for k, v in element_diff.items()])

            # self._defect_energies[name][charge] = \
            #     relative_energy + correction_energy + \
            #     electron_interchange_energy + element_interchange_energy

    # def calc_concentration(self, temperature, e_f, num_sites_filename):
    #     """
    #     Calculates defect concentrations. When the self._original_concent is set,
    #     the each defect calc_concentration is kept fixed. For example, the sum of
    #     Va_O1_0, Va_O1_1 and Va_O1_2 are fixed.

        # Args:
        #     temperature (float): Temperature in K.
        #     e_f (float): Fermi level in the absolute scale.
        #     num_sites_filename (str): Yaml file name for the number of sites in
        #                               a given volume. Example is
        #                                 Va_Mg1 : 2
        #                                 Va_O1 : 2
        #                                 Mg_i1 : 4
        # Return:
        #     concentrations[defect name][charge]:
        # """

        # num_sites = loadfn(num_sites_filename)
        # concentrations = defaultdict(dict)

        # for name in self._defect_energies.keys():
        #     if name not in num_sites.keys():
        #         raise KeyError(
        #             "{} is not in {}.".format(name, num_sites_filename))

            # if self._original_concent is None:
            #     for charge in self._defect_energies[name].keys():

            #         energy = self._defect_energies[name][charge] + e_f * charge
            #         concentrations[name][charge] = \
            #             maxwell_boltzmann_dist(energy, temperature) \
            #             / self._volume * num_sites[name]
            # else:
            #     if self._quenched_concent:
            #         raise ValueError(
            #             "Concentrations at the quenched temperature is already "
            #             "set. Unset the concentrations first.")

                # dc = self._original_concent.calc_concentration[name]
                # total = sum(dc.values())
                # ratio = {}
                # for charge in self._defect_energies[name].keys():
                #     energy = self._defect_energies[name][charge] + e_f * charge
                #     ratio[charge] = \
                #         maxwell_boltzmann_dist(energy, temperature) \
                #         / self._volume * num_sites[name]
                # ratio_sum = sum(ratio.values())
                # for charge in self._defect_energies[name].keys():
                #     concentrations[name][charge] = \
                #         total * ratio[charge] / ratio_sum

        # return concentrations

    # def calc_equilibrium_concentration(self, temperature, num_sites_filename,
    #                               max_iteration=50, threshold=1e-5,
    #                               initialize=False):
    #     """
    #     Calculates the equilibrium defect concentrations

        # Args:
        #     temperature (float): Temperature in K.
        #     num_sites_filename (str): Yaml file name for the number of sites in
        #                               a given volume.
        #     max_iteration (int): Max iteration number for seeking the
        #                          equilibrium defect concentrations
        #     threshold (float): Threshold for the convergence of the SCF
        #                        calculation of the carrier/defect concentrations.
        #     initialize (bool): Whether to reset the defect calc_concentration
        # """
        # if initialize is True:
        #     self.unset_concentration()

        # # In case the Fermi level locates in between vbm and cbm, the common
        # distance = self._cbm - self._vbm
        # e_f_abs = (self._vbm + self._cbm) / 2
        # for iteration in range(max_iteration):

            # # Use absolute Fermi level for p and n.
            # p = CarrierConcentration.p(temperature, e_f_abs, self._total_dos,
            #                            self._vbm, self._volume)
            # n = CarrierConcentration.n(temperature, e_f_abs, self._total_dos,
            #                            self._cbm, self._volume)
            # e_f = e_f_abs - self._vbm
            # calc_concentration = \
            #     self.calc_concentration(temperature, e_f, num_sites_filename)
            # defect_charge = \
            #     sum([sum([c * e for c, e in calc_concentration[d].items()])
            #          for d in calc_concentration])
            # charge_sum = p - n + defect_charge

            # # In case the Fermi level locates in between vbm and cbm, the common
            # # ratio 0.5 is okay. Otherwise, higher common ratio is needed.
            # # Therefore, 0.7 is set here.
            # distance *= 0.7
            # e_f_abs = e_f_abs + np.sign(charge_sum) * distance
            # # This line controls the accuracy.
            # max_concentration = np.amax([n, p, charge_sum])
            # if np.abs(charge_sum / max_concentration) < threshold:
            #     c = Concentration(temperature, e_f, n, p, calc_concentration)
            #     if self._original_concent is None:
            #         self._original_concent = c
            #     else:
            #         self._quenched_concent = c
            #     return True

        # print("Equilibrium condition has not been reached. ")
        # return False

    # def unset_concentration(self):
    #     self._original_concent = None
    #     self._quenched_concent = None

    # def set_vbm2_cbm2(self, vbm2, cbm2):
    #     self._vbm2 = vbm2
    #     self._cbm2 = cbm2

    # @staticmethod
    # def min_e_at_ef(ec, ef):
    #     # d[c1] = energy for charge c1
    #     d = {c: e + c * ef for c, e in ec.items()}
    #     # return charge for the lowest energy, and its energy value
    #     return min(d.items(), key=lambda x: x[1])

    # def calc_transition_levels(self, x_range=None):
    #     """
    #     Calculates a set of transition levels

        # Parameter in use:
        #     TransitionLevel (namedtuple):
        # """

        # transition_levels = {}

        # if x_range:
        #     x_min = x_range[0]
        #     x_max = x_range[1]
        # else:
        #     x_min = 0.0
        #     x_max = self.band_gap

        # # e_of_c[charge] = energy
        # for name, e_of_c in self._defect_energies.items():
        #     points = []
        #     charge = set()

            # # Estimate the lowest energies at the lowest Fermi level.
            # c_vbm, min_e_vbm = self.min_e_at_ef(e_of_c, x_min)
            # points.append([x_min, min_e_vbm])
            # charge.add(c_vbm)

            # # and at the highest Fermi level
            # c_cbm, min_e_cbm = self.min_e_at_ef(e_of_c, x_max)
            # points.append([x_max, min_e_cbm])
            # charge.add(c_cbm)

            # for (c1, e1), (c2, e2) in combinations(e_of_c.items(), r=2):
            #     # Estimate the cross point between two charge states
            #     x = - (e1 - e2) / (c1 - c2)
            #     y = (c1 * e2 - c2 * e1) / (c1 - c2)

                # compared_energy = self.min_e_at_ef(e_of_c, x)[1] + 0.00001

                # if x_min < x < x_max and y < compared_energy:
                #     points.append([x, y])
                #     charge.add(c1)
                #     charge.add(c2)

            # transition_levels[name] = \
            #     TransitionLevel(cross_points=sorted(points, key=lambda x: x[0]),
            #                     charges=sorted(charge))

        # self._transition_levels = transition_levels
        # self._calculated_transition_level_range = [x_min, x_max]

    # def plot_energy(self, file_name=None, x_range=None, y_range=None,
    #                 show_tls=False):
    #     """
    #     Plots the defect formation energies as a function of the Fermi level.

        # Args:
        #     file_name (str): File name for saving the plot.
        #     x_range (2x1 list): x range for the plot.
        #     y_range (2x1 list): y range for the plot.
        #     show_tls (bool): If transition level values are shown or not.
        # """

        # fig, ax = plt.subplots()

        # plt.title(self._title, fontsize=15)

        # ax.set_xlabel("Fermi level (eV)", fontsize=15)
        # ax.set_ylabel("Formation energy (eV)", fontsize=15)

        # if self.transition_levels is None:
        #     raise NoTLError("Transition levels are not calculated yet.")

        # if x_range:
        #     x_min = x_range[0]
        #     x_max = x_range[1]
        # else:
        #     x_min = 0
        #     x_max = self.band_gap

        # if y_range:
        #     y_min = y_range[0]
        #     y_max = y_range[1]
        # else:
        #     y_min = min([min(np.array(tl.cross_points).transpose()[1])
        #                  for tl in self.transition_levels.values()])
        #     y_max = max([max(np.array(tl.cross_points).transpose()[1])
        #                  for tl in self.transition_levels.values()])
        #     if y_min > 0:
        #         y_min = 0
        #     # Make top and bottom space
        #     margin = (y_max - y_min) * 0.08
        #     y_min -= margin
        #     y_max += margin

        # plt.xlim(x_min, x_max)
        # plt.ylim(y_min, y_max)

        # # support lines
        # if x_range:
        #     plt.axvline(x=0, linewidth=1.0, linestyle='dashed')
        #     plt.axvline(x=self.band_gap, linewidth=1.0, linestyle='dashed')
        #     ax.annotate("cbm", (self.band_gap, y_min * 0.9 + y_max * 0.1),
        #                 fontsize=10)

            # perfect_vbm = self._perfect_vbm - self._vbm
            # perfect_cbm = self._perfect_cbm - self._vbm
            # print(self._perfect_vbm)
            # print(perfect_vbm)
            # print(perfect_cbm)
            # if perfect_vbm < - 0.05:
            #     plt.axvline(x=perfect_vbm, linewidth=1.0, linestyle='dotted')
            #     ax.annotate("supercell vbm",
            #                 (perfect_vbm, y_min * 0.9 + y_max * 0.1),
            #                 fontsize=10)
            # if perfect_cbm > self.band_gap + 0.05:
            #     plt.axvline(x=perfect_cbm, linewidth=1.0, linestyle='dotted')
            #     ax.annotate("supercell cbm",
            #                 (perfect_cbm, y_min * 0.9 + y_max * 0.1),
            #                 fontsize=10)

        # plt.axhline(y=0, linewidth=1.0, linestyle='dashed')

        # # Lines for the equilibrium Fermi level(s)
        # if self._original_concent:
        #     t1 = self._original_concent.temperature
        #     x1 = self._original_concent.e_f
        #     plt.axvline(x=x1, linewidth=2.0, linestyle='dashed', color="red")
        #     ax.annotate("T$_1$=" + str(t1) + "K", (x1, y_min), fontsize=7)
        #     if self._quenched_concent:
        #         t2 = self._quenched_concent.temperature
        #         x2 = self._quenched_concent.e_f
        #         plt.axvline(x=x2, linewidth=2.0, linestyle='dashed',
        #                     color="purple")
        #         ax.annotate("T$_2$=" + str(t2) + "K", (x2, y_min), fontsize=7)

        # for i, (name, tl) in enumerate(self.transition_levels.items()):

            # color = matplotlib.cm.hot(float(i) / len(self.transition_levels))

            # cross_points = tl.cross_points
            # shallow = []
            # transition_levels = []

            # # Store the coords for the filled and open circles.
            # for cp in cross_points:
            #     tl_x_min = self._calculated_transition_level_range[0]
            #     tl_x_max = self._calculated_transition_level_range[1]
            #     if (cp[0] < tl_x_min + 0.0001 and max(tl.charges) < 0) \
            #             or (cp[0] > tl_x_max - 0.0001 and min(tl.charges) > 0):
            #         shallow.append(cp)
            #     elif not (cp[0] == tl_x_min or cp[0] == tl_x_max):
            #         transition_levels.append(cp)

            # # set the x and y arrays to be compatible with matplotlib style.
            # # e.g., x = [0.0, 0.3, 0.5, ...], y = [2.1, 3.2, 1.2, ...]
            # # plot lines
            # x, y = np.array(cross_points).transpose()
            # line, = ax.plot(x, y, '-', color=color, label=name)
            # line.set_label(name)

            # # plot filled circles for transition levels.
            # if transition_levels:
            #     x, y = np.array(transition_levels).transpose()
            #     ax.scatter(x, y, color=color)

            # # plot unfilled circles for shallow levels.
            # if shallow:
            #     x, y = np.array(shallow).transpose()
            #     ax.scatter(x, y, facecolor="white", edgecolor=color)

            # # These determine the positions of the transition levels.
            # margin_x = (x_max - x_min) * 0.025
            # margin_y = (y_max - y_min) * 0.025
            # if show_tls:
            #     for cp in cross_points:
            #         if cp[0] < x_min + 0.0001 or cp[0] > x_max - 0.0001:
            #             continue
            #         s = str(round(cp[0], 2)) + ", " + str(round(cp[1], 2))
            #         ax.annotate(s, (cp[0] + margin_x, cp[1] - margin_y),
            #                     color=color, fontsize=8)

            # # Arrange the charge states at the middle of the transition levels.
            # middle_points = \
            #     reversed([[(a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + margin_y]
            #               for a, b in zip(cross_points, cross_points[1:])])

            # for j, (x, y) in enumerate(middle_points):
            #     ax.annotate(str(tl.charges[j]), (x, y), color=color,
            #                 fontsize=13)

        # ax.legend()
        # fig.subplots_adjust(right=0.75)

        # if file_name:
        #     plt.savefig(file_name)
        # else:
        #     plt.show()

    # def show_concentration(self):
    #     if self._quenched_concent:
    #         cs = [self._original_concent, self._quenched_concent]
    #     else:
    #         cs = [self._original_concent]

        # for c in cs:
        #     print("--------")
        #     print("Tempereture: {} K.".format(c.temperature))
        #     print("Fermi level: {} eV.".format(c.e_f))
        #     print("p: {:.1e} cm-3.".format(c.p))
        #     print("n: {:.1e} cm-3.".format(c.n))
        #     print("p - n: {:.1e} cm-3.".format(c.p - c.n))
        #     for name, c_of_charge in c.calc_concentration.items():
        #         print("---")
        #         for charge, calc_concentration in c_of_charge.items():
        #             print("{} {}: {:.1e} cm-3.".format(name, charge,
        #                                                calc_concentration))

    # @property
    # def defect_energies(self):
    #     return self._defect_energies

    # @property
    # def transition_levels(self):
    #     return self._transition_levels

    # @property
    # def vbm(self):
    #     return self._vbm

    # @property
    # def cbm(self):
    #     return self._cbm

    # @property
    # def band_gap(self):
    #     return self._cbm - self._vbm


