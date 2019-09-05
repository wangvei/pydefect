# -*- coding: utf-8 -*-

import json
from itertools import combinations
from itertools import groupby
from operator import itemgetter, attrgetter
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn
from pydefect.analysis.defect import Defect
from pydefect.core.config import COLOR
from pydefect.core.defect_name import DefectName
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.logger import get_logger
from pydefect.util.tools import (
    construct_obj_in_dict, sanitize_keys_in_dict, defaultdict_to_dict,
    flatten_dict)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class DefectEnergy(MSONable):
    def __init__(self,
                 defect_energy: float,
                 annotation: str,
                 multiplicity: int,
                 magnetization: float,
                 convergence: bool,
                 shallow: bool):
        """

        Args:
            defect_energy (float):
                Defect formation energy at Ef = 0.
            annotation (str):
                Defect annotation
            multiplicity (dict):
                Spatial multiplicity.
            magnetization (dict):
                Total magnetization in mu..
            convergence (bool):
                Ionic convergence.
            shallow (bool):
                Whether it is shallow or not.
        """
        self.defect_energy = defect_energy
        self.annotation = annotation
        self.multiplicity = multiplicity
        self.magnetization = magnetization
        self.convergence = convergence
        self.shallow = shallow


class DefectEnergies(MSONable):
    def __init__(self,
                 defect_energies: Dict[str, Dict[int, DefectEnergy]],
                 transition_levels: Dict[str, Dict[str, list]],
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 title: str = None):
        """ A class related to a set of defect formation energies.

        Args:
            defect_energies (dict):
                DefectEnergy as a function of name, charge, and annotation.
                energies[name][charge] = DefectEnergy object
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
        self.transition_levels = transition_levels
        self.vbm = vbm
        self.cbm = cbm
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.title = title

    @classmethod
    def from_objects(cls,
                     unitcell: UnitcellCalcResults,
                     perfect: SupercellCalcResults,
                     defects: List[Defect],
                     chem_pot: tuple,
                     chem_pot_label: str,
                     filtering_words: list = None,
                     exclude_unconverged_defects: bool = True,
                     exclude_shallow_defects: bool = True,
                     system: str = None):
        """ Calculates defect formation energies from several objects.

        All the energies are calculated at 0 eV in the absolute scale.

        Args:
            unitcell (UnitcellCalcResults):
                UnitcellCalcResults object for band edge.
            perfect (SupercellCalcResults):
                SupercellDftResults object of perfect supercell for band edge
                in supercell.
            defects (list of Defect objects):
                List of Defect objects, each of which has "defect_entry",
                "dft_results", and "correction" attributes,
            chem_pot (tuple):
                Return of ChemPotDiag.load_vertices_yaml method
            chem_pot_label (str):
                Equilibrium point specified in ChemPot.
            filtering_words (list):
                List of words used for filtering the defects
            exclude_unconverged_defects (bool):
                Whether to exclude the unconverged defects from the plot.
            exclude_shallow_defects (bool):
                Whether to exclude the shallow defects from the plot.
            system (str):
                System name used for the title.
        """
        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute scale.
        vbm, cbm = unitcell.band_edge

        system = system or perfect.final_structure.composition.reduced_formula
        title = f"{system} at condition {chem_pot_label}"

        # Chemical potentials
        relative_chem_pots, standard_e = chem_pot
        relative_chem_pot = relative_chem_pots[chem_pot_label]

        # defect_energies[name][charge] = DefectEnergy
        defect_energies = {}
        transition_levels = {}

        defects.sort(key=attrgetter("name"))
        for name, g_name in groupby(defects, key=attrgetter("name")):
            energy_by_c = {}
            g_name = sorted(list(g_name), key=attrgetter("charge"))
            for charge, g_charge in groupby(g_name, key=attrgetter("charge")):
                energy_by_annotation = {}

                for d in g_charge:
                    n = DefectName(name, charge, d.annotation)

                    if n.is_name_matched(filtering_words) is False:
                        logger.info(f"{n} filtered out, so omitted.")
                        continue
                    elif exclude_shallow_defects and d.is_shallow:
                        logger.info(f"{n} is shallow, so omitted.")
                        continue
                    elif exclude_unconverged_defects and not d.is_converged:
                        logger.info(f"{n} unconverged, so omitted.")
                        continue

                    energy_by_annotation[d] = \
                        d.relative_total_energy + d.correction_energy

                if not energy_by_annotation:
                    continue

                # get the lowest energy and its annotation
                defect, defect_energy = \
                    min(energy_by_annotation.items(), key=itemgetter(1))

                for el, diff in defect.changes_of_num_elements.items():
                    relative_e = relative_chem_pot.elem_coords[el]
                    standard = standard_e[el]
                    defect_energy -= diff * (relative_e + standard)

                energy_by_c[charge] = \
                    DefectEnergy(defect_energy=defect_energy,
                                 annotation=defect.annotation,
                                 multiplicity=defect.final_multiplicity,
                                 magnetization=defect.magnetization,
                                 convergence=defect.is_converged,
                                 shallow=defect.is_shallow)

            cross_points_with_charges = []
            for (c1, e1), (c2, e2) in combinations(energy_by_c.items(), r=2):
                # The cross point between two charge states.
                x = - (e1.defect_energy - e2.defect_energy) / (c1 - c2)
                y = (c1 * e2.defect_energy - c2 * e1.defect_energy) / (c1 - c2)

                # The lowest energy among all the charge states to be compared.
                compared_energy = min([e.defect_energy + c * x
                                       for c, e in energy_by_c.items()])

                if y < compared_energy + 1e-5:
                    cross_points_with_charges.append([[x, y], [c1, c2]])

            # need to sort the points along x-axis.
            cpwc = sorted(cross_points_with_charges, key=lambda z: z[0][0])
            transition_levels[name] = {"cross_points": [i[0] for i in cpwc],
                                       "charges": [i[1] for i in cpwc]}
            defect_energies[name] = energy_by_c

        defect_energies = defaultdict_to_dict(defect_energies)

        return cls(defect_energies=defect_energies,
                   transition_levels=transition_levels,
                   vbm=vbm,
                   cbm=cbm,
                   supercell_vbm=perfect.vbm,
                   supercell_cbm=perfect.cbm,
                   title=title)

    @classmethod
    def from_dict(cls, d):
        """ Construct a class object from a dictionary. """

        defect_energies = sanitize_keys_in_dict(d["defect_energies"])
        defect_energies = construct_obj_in_dict(defect_energies, DefectEnergy)

        return cls(defect_energies=defect_energies,
                   transition_levels=d["transition_levels"],
                   vbm=d["vbm"],
                   cbm=d["cbm"],
                   supercell_vbm=d["supercell_vbm"],
                   supercell_cbm=d["supercell_cbm"],
                   title=d["title"])

    def to_json_file(self, filename="defect_energies.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    def __repr__(self):
        outs = []
        for name, charge, defect_energy in flatten_dict(self.defect_energies):

            energy_vbm = defect_energy.defect_energy + charge * self.vbm
            energy_cbm = defect_energy.defect_energy + charge * self.cbm
            outs.extend(
                [f"name: {name}:",
                 f"charge: {charge}",
                 f"Annotation: {defect_energy.annotation}",
                 f"Energy @ vbm (eV): {round(energy_vbm, 4)}",
                 f"Energy @ cbm (eV): {round(energy_cbm, 4)}",
                 f"Convergence: {defect_energy.convergence}",
                 f"Is shallow: {defect_energy.shallow}",
                 f"multiplicity: {defect_energy.multiplicity}",
                 f"magnetization: {defect_energy.magnetization}"])
            outs.append("")

        return "\n".join(outs)

    def u(self, name: str, charges: list) -> Tuple[float, list]:
        """ Return the U value among three sequential charge states.

        Args:
            name (str):
                Name of the defect.
            charges (list):
                1x3 list comprising three charge states.
        Return:
            U value.
        """
        if len(charges) != 3:
            assert ValueError("Length of charge states must be 3.")
        elif charges[2] - charges[1] != 1 or charges[1] - charges[0] != 1:
            assert ValueError("The charge states {} {} {} are not sequential."
                              .format(*charges))
        elif not charges[0] < charges[1] < charges[2]:
            assert ValueError("The charge states {} {} {} are not incremental."
                              .format(*charges))

        energies = []
        names = []
        for charge in charges:
            defect = self.defect_energies[name][charge]
            energies.append(defect.defect_energy)
            names.append(str(DefectName(name, charge, defect.annotation)))

        return energies[0] + energies[2] - 2 * energies[1], names

    @property
    def band_gap(self):
        return self.cbm - self.vbm

    def plot_energy(self,
                    x_range: list = None,
                    y_range: list = None,
                    fermi_levels: list = None,
                    show_transition_levels: bool = False,
                    show_all_energies: bool = False,
                    color: list = None,
                    compile: bool = False):
        """ Plots defect formation energies as a function of the Fermi level.

        Args:
            x_range (list):
                1x2 list for determining x range.
            y_range (list):
                1x2 list for determining y range.
            fermi_levels (list):
                Ghe Fermi level(s) (E_f) to be shown in the plot.
                The structure must be as follows.
                [[temperature, E_f]] or
                [[temperature, E_f], [quenched temperature, E_f]]
            show_transition_levels (bool):
                Whether to show values of transition levels in the plot.
            show_all_energies (bool):
                Whether to show all energies in the plot.
            color (list)
                User favorite color scheme.
        """
        if compile:
            fig, (ax, ax2) = plt.subplots(2, sharex=True)
        else:
            fig, ax = plt.subplots()
        plt.title(self.title, fontsize=15)

        color = color or COLOR

        ax.set_xlabel("Fermi level (eV)", fontsize=15)
        ax.set_ylabel("Formation energy (eV)", fontsize=15)

        x_min, x_max = x_range if x_range else (0, self.band_gap)
        y_min, y_max = float("inf"), -float("inf")

        def min_e_at_ef(ec: Dict[int, DefectEnergy], ef):
            # calculate each energy at the given Fermi level ef.
            d = {cc: ee.defect_energy + cc * ef for cc, ee in ec.items()}
            # return the charge with the lowest energy, and its energy value
            return min(d.items(), key=itemgetter(1))

        # Note: len(self.energies) <= len(transition_levels)
        for i, (name, tl) in enumerate(self.transition_levels.items()):
            # lowest_energies[name] is {} when all results are omitted above.
            line_type = '-' if i < 5 else '--' if i < 10 else '-.'
            # ---------- Calculate cross points including edges ---------------
            cross_points = []
            charge_set = set()

            # keep x_min -> transition levels -> x_max for connecting points.
            # Add the plot point at x_min
            charge, y = min_e_at_ef(self.defect_energies[name], x_min + self.vbm)
            cross_points.append([x_min, y])
            # Adding charge is a must for plot
            charge_set.add(charge)
            y_min, y_max = min([y, y_min]), max([y, y_max])
            # Unfilled and filled circles for shallow and transition levels
            if max(charge_set) < 0:
                ax.plot(x_min, y, marker="o", mec=color[i], mfc="white")

            # Add points between x_min and x_max
            for cp, charges in zip(tl["cross_points"], tl["charges"]):
                if x_min < cp[0] - self.vbm < x_max:
                    cross_points.append([cp[0] - self.vbm, cp[1]])
                    # need to sort the charge.
                    charge_set.add(sorted(charges)[1])
                    y_min, y_max = min([cp[1], y_min]), max([cp[1], y_max])
                    ax.scatter(cp[0] - self.vbm, cp[1], marker='o',
                               color=color[i])

            # Add the plot point at x_max
            charge, y = min_e_at_ef(self.defect_energies[name], x_max + self.vbm)
            cross_points.append([x_max, y])
            charge_set.add(charge)
            y_min, y_max = min([y, y_min]), max([y, y_max])
            if min(charge_set) > 0:
                ax.plot(x_max, y, marker="o", mec=color[i], mfc="white")

            # -----------------------------------------------------------------
            # set x and y arrays to be compatible with matplotlib style.
            # e.g., x = [0.0, 0.3, 0.5, ...], y = [2.1, 3.2, 1.2, ...]
            x, y = np.array(cross_points).transpose()
            line, = ax.plot(x, y, line_type, color=color[i], label=name)
            line.set_label(name)

            # margin_x and _y determine the positions of the transition levels.
            margin_y = (y_max - y_min) * 0.1

            if show_transition_levels:
                for cp in cross_points:
                    if x_min + 1e-4 < cp[0] < x_max - 1e-4 is False:
                        continue
                    s = f"{round(cp[0], 2)}, {round(cp[1], 2)}"
                    pos_x = cp[0]
                    pos_y = cp[1] - margin_y
                    ax.annotate(s, (pos_x, pos_y), color=color[i], fontsize=10)

            # Arrange the charge states at the middle of the transition levels.
            margin_y = 0
            charge_pos = [[(a[0] + b[0]) / 2, (a[1] + b[1] + margin_y) / 2]
                          for a, b in zip(cross_points, cross_points[1:])]

            sorted_charge_set = sorted(charge_set, reverse=True)

            for j, (x, y) in enumerate(charge_pos):
                charge = sorted_charge_set[j]
                annotation = self.defect_energies[name][charge].annotation
                chg = f"{charge}: {annotation}" if annotation else str(charge)

                ax.annotate(
                    chg, (x, y), color=color[i], fontsize=9,
                    bbox=dict(facecolor='white', edgecolor=color[i], pad=1))

            if show_all_energies:
                for c, e in self.defect_energies[name].items():
                    y1 = e.defect_energy + c * (x_min + self.vbm)
                    y2 = e.defect_energy + c * (x_max + self.vbm)
                    ax.plot([x_min, x_max], [y1, y2], line_type, linewidth=0.3,
                            color=color[i])

        if y_range:
            y_min, y_max = y_range
        else:
            y_min = 0 if y_min > 0 else y_min
            margin = (y_max - y_min) * 0.1
            y_min, y_max = y_min - margin, y_max + margin

        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

        # plot vbm, cbm, supercell_vbm, supercell_cbm lines.
        plt.axvline(x=0, linewidth=1.0, linestyle="-", color='b')
        plt.axvline(x=self.band_gap, linewidth=1.0, linestyle="-", color='b')
        ax.annotate("vbm", (0, y_min * 0.9 + y_max * 0.1),
                    fontsize=10, color='b')
        ax.annotate("cbm", (self.band_gap, y_min * 0.9 + y_max * 0.1),
                    fontsize=10, color='b')

        supercell_vbm = self.supercell_vbm - self.vbm
        supercell_cbm = self.supercell_cbm - self.vbm

#        if supercell_vbm < - 0.05:
        plt.axvline(x=supercell_vbm, linewidth=1.0, linestyle='-', color='r')
        ax.annotate("supercell vbm",
                    (supercell_vbm, y_min * 0.8 + y_max * 0.2),
                    fontsize=10, color='r')

#        if supercell_cbm > self.band_gap + 0.05:
        plt.axvline(x=supercell_cbm, linewidth=1.0, linestyle='-', color='r')
        ax.annotate("supercell cbm",
                    (supercell_cbm, y_min * 0.8 + y_max * 0.2),
                    fontsize=10, color='r')

        plt.axhline(y=0, linewidth=1.0, linestyle='dashed')

        if fermi_levels:
            y = y_min * 0.85 + y_max * 0.15
            for i, (t, f) in enumerate(fermi_levels):
                plt.axvline(x=f - self.vbm, linewidth=1.0, linestyle=':',
                            color='g')
                ax.annotate(t, ((f - self.vbm), y), fontsize=10, color='green')
                if i == 1:
                    y = y_min * 0.88 + y_max * 0.12
                    plt.arrow(x=fermi_levels[0][1] - self.vbm, y=y,
                              dx=f - fermi_levels[0][1], dy=0,
                              width=0.005,
                              head_width=0.5,
                              head_length=0.1,
                              length_includes_head=True,
                              color='green')

        # Shrink current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height * 0.9])
        ax.legend(bbox_to_anchor=(1, 0.5), loc='center left')
        #        fig.subplots_adjust(right=0.75)

        if compile:
            return plt, ax2
        else:
            return plt
