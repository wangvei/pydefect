# -*- coding: utf-8 -*-

from collections import defaultdict
from copy import copy
from itertools import combinations
import json
from monty.json import MontyEncoder, MSONable
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from typing import List

from obadb.analyzer.chempotdiag.chem_pot_diag import ChemPotDiag
from pydefect.analysis.defect import Defect

from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_name import DefectName
from pydefect.util.logger import get_logger


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class TransitionLevel:
    def __init__(self,
                 cross_points: list,
                 charges: list):
        self.cross_points = cross_points
        self.charges = charges


class DefectEnergy(MSONable):
    def __init__(self,
                 energy: float,
                 convergence: bool,
                 is_shallow: bool):
        """
            energy (dict):
            convergence (dict):
            is_shallow (dict):
                Whether defect is shallow or not.
        """
        self.energy = energy
        self.convergence = convergence
        self.is_shallow = is_shallow


class DefectEnergies(MSONable):
    def __init__(self,
                 defect_energies: dict,
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 title: str = None):
        """ A class related to a set of defect formation energies.

        Args:
            defect_energies (dict):
                DefectEnergy as a function of name, charge, and annotation.
                defect_energies[name][charge][annotation] = DefectEnergy object
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            supercell_vbm (float):
                Valence band maximum in the perfect supercell.
            supercell_cbm (float):
                Conduction band minimum in the perfect supercell.
            title (str):
                Title of the system.
        """
        self.defect_energies = defect_energies
        self.vbm = vbm
        self.cbm = cbm
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.title = title

    @classmethod
    def from_files(cls,
                   unitcell: UnitcellCalcResults,
                   perfect: SupercellCalcResults,
                   defects: List[Defect],
                   chem_pot: ChemPotDiag,
                   chem_pot_label: str,
                   system: str = ""):
        """ Calculates defect formation energies from several objects.

        Note that all the energies are calculated at 0 eV in the absolute scale.
        Args:
            unitcell (UnitcellCalcResults):
                UnitcellCalcResults object for band edge.
            perfect (SupercellCalcResults):
                SupercellDftResults object of perfect supercell for band edge in
                supercell.
            defects (list of namedtuple Defect):
                List of the Defect namedtuple object.
                Defect = namedtuple(
                    "Defect", "defect_entry", "dft_results", "correction")
            chem_pot (ChemPot):
                Chemical potentials of the competing phases.
            chem_pot_label (str):
                Equilibrium point specified in ChemPot.
            system (str):
                System name used for the title.
        """
        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy.
        vbm, cbm = unitcell.band_edge
        supercell_cbm, supercell_vbm = perfect.eigenvalue_properties[1:3]

        title = system + " condition " + chem_pot_label

        # Chemical potentials
        relative_chem_pots, standard_e = chem_pot
        relative_chem_pot = relative_chem_pots[chem_pot_label]

        defect_energies = defaultdict(dict)

        for d in defects:
            name = d.defect_entry.name

            charge = d.defect_entry.charge
            annotation = d.defect_entry.annotation

            element_diff = d.defect_entry.changes_of_num_elements

            # Calculate defect formation energies at the vbm
            relative_energy = d.dft_results.relative_total_energy
            correction_energy = d.correction.correction_energy
            element_interchange_energy = \
                - sum([v * (relative_chem_pot.elem_coords[k] + standard_e[k])
                       for k, v in element_diff.items()])
            energy = \
                relative_energy + correction_energy + element_interchange_energy

            e = DefectEnergy(energy=energy,
                             convergence=d.dft_results.is_converged,
                             is_shallow=d.is_shallow)

            defect_energies[name].setdefault(charge, {}).update({annotation: e})

        return cls(defect_energies=dict(defect_energies),
                   vbm=vbm,
                   cbm=cbm,
                   supercell_vbm=supercell_vbm,
                   supercell_cbm=supercell_cbm,
                   title=title)

    def to_json_file(self, filename="defect_energy.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def __str__(self):
        outs = []
        for n in self.defect_energies.keys():
            for c in self.defect_energies[n].keys():
                for a, de in self.defect_energies[n][c].items():
                    outs.extend(
                        ["name: {}:".format(n),
                         "charge: {}".format(c),
                         "Annotation: {}".format(a),
                         "Energy @ ef = 0 (eV): {}".format(round(de.energy, 4)),
                         "Convergence: {}".format(de.convergence),
                         "Is shallow: {}".format(de.is_shallow)])
                    outs.append("")

        return "\n".join(outs)

    def u(self,
          name: str,
          charges: list,
          annotations: list = None):
        """ Return the U value among three sequential charge states.

        Args:
            name (str):
                Name of the defect.
            charges (list):
                1x3 list comprising three charge states.
            annotations (list)
                Assign specific annotations
        Return:
            U value.
        """
        if len(charges) != 3:
            assert ValueError("Length of charge states must be 3.")
        elif charges[2] - charges[1] != 1 or charges[1] - charges[0] != 1:
            assert ValueError("The charge states {} {} {} are not sequential."
                              .format(*charges))
        elif not charges[0] < charges[1] < charges[2]:
            assert ValueError("The charge states {} {} {} are not incremetal."
                              .format(*charges))

        energies = []
        if annotations is None:
            # annotations with lowest energies at given charge states
            annotations = []
            for c in charges:
                annotation, min_e = min(self.defect_energies[name][c].items(),
                                        key=lambda x: x[1])
                energies.append(min_e)
                annotations.append(annotation)
        else:
            for c, a in zip(charges, annotations):
                energies.append(self.defect_energies[name][c][a].energy)

        return energies[0] + energies[2] - 2 * energies[1], annotations

    @property
    def band_gap(self):
        return self.cbm - self.vbm

    def plot_energy(self,
                    filtering_words: list = None,
                    x_range: list = None,
                    y_range: list = None,
                    exclude_unconverged_defects: bool = True,
                    exclude_shallow_defects: bool = True,
                    show_fermi_level: bool = True,
                    show_transition_levels: bool = False,
                    show_all_energies: bool = False):
        """ Plots defect formation energies as a function of the Fermi level.

        Args:
            filtering_words (list):
                List of words used for filtering the defects
            x_range (list):
                1x2 list for determining x range.
            y_range (list):
                1x2 list for determining y range.
            exclude_unconverged_defects (bool):
                Whether to exclude the unconverged defects from the plot.
            exclude_shallow_defects (bool):
                Whether to exclude the shallow defects from the plot.
            show_fermi_level (bool):
                Whether to show the Fermi level(s) in the plot.
            show_transition_levels (bool):
                Whether to show values of transition levels in the plot.
            show_all_energies (bool):
                Whether to show all energies in the plot.
        """
        fig, ax = plt.subplots()
        plt.title(self.title, fontsize=15)

        ax.set_xlabel("Fermi level (eV)", fontsize=15)
        ax.set_ylabel("Formation energy (eV)", fontsize=15)

        # e_of_c is energy as a function of charge: e_of_c[charge] = energy
        transition_levels = []
        for name, e_of_c in copy(self.defect_energies).items():
            lowest_energy_at_name = []
            lowest_energy_annotation = []
            for c, e_of_a in copy(e_of_c).items():

                for a, e in copy(e_of_a).items():
                    n = DefectName(name=name, charge=c)
                    if n.is_name_matched(filtering_words) is False:
                        e_of_a.pop(c)
                    elif exclude_shallow_defects and self.defect_energies[name][c]:
                        logger.warning("{} is shallow, so omitted.".format(n))
                        e_of_a.pop(c)
                    elif (exclude_unconverged_defects
                          and self.convergence[name][c] is False):
                        logger.warning("{} is unconverged, so omitted.".format(n))
                    e_of_c.pop(c)

                annotation, min_e = min(e_of_a.items(), key=lambda x: x[1])
                lowest_energy_at_name[c] = min_e
                lowest_energy_annotation[c] = annotation

            cross_points_charge = []

            for (c1, e1), (c2, e2) in combinations(lowest_energy_at_name.items(), r=2):
                # The cross point between two charge states.
                x = - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                # The lowest energy among all the charge states to be compared.
                compared_energy = \
                    min([energy + c * x for c, energy in lowest_energy_at_name.items()])

                if y < compared_energy + 1e-5:
                    cross_points_charge.append([[x, y], [c1, c2]])

            # need to sort the points along x-axis.
            cross_points = []
            charge = []
            for cp, c in sorted(cross_points_charge, key=lambda z: z[0][0]):
                cross_points.append(cp)
                charge.append(c)

            transition_levels[name] = \
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

        for i, (name, tl) in enumerate(transition_levels.items()):

            color = matplotlib.cm.hot(float(i) / len(transition_levels))

            # ---------- Calculate cross points including edges ----------------
            energies = self.energies[name]
            cross_points = []
            shallow_points = []
            transition_points = []
            charge_set = set()

            # must keep sequence (x_min -> transition levels -> x_max) for
            # connecting the points.

            # Add the plot point at x_min
            c, y = min_e_at_ef(energies, x_min + self.vbm)
            cross_points.append([x_min, y])
            # Adding charge is a must
            charge_set.add(c)
            if c < 0:
                shallow_points.append([x_min, y])
            if y < y_min:
                y_min = y
            if y > y_max:
                y_max = y

            # Add points between x_min and x_max
            for cp, c in zip(tl.cross_points, tl.charges):
                if x_min < cp[0] - self.vbm < x_max:
                    cross_points.append([cp[0] - self.vbm, cp[1]])
                    transition_points.append([cp[0] - self.vbm, cp[1]])
                    # need to sort the charge.
                    charge_set.add(sorted(c, key=lambda x: x)[1])
                    if cp[1] < y_min:
                        y_min = cp[1]
                    if cp[1] > y_max:
                        y_max = cp[1]


            # Add the plot point at x_max
            c, y = min_e_at_ef(energies, x_max + self.vbm)
            cross_points.append([x_max, y])
            charge_set.add(c)
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

            # if show_all_energies:
            #     for c, e in self.energies[name].items():
            #         y1 = e + c * (x_min + self.vbm)
            #         y2 = e + c * (x_max + self.vbm)
            #         ax.plot([x_min, x_max], [y1, y2], '-', linewidth=0.3,
            #                 color=color)

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
