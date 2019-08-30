# -*- coding: utf-8 -*-

import json
from collections import defaultdict
from copy import copy
from itertools import combinations
from itertools import groupby
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn
from pydefect.analysis.defect import Defect
from pydefect.core.defect_name import DefectName
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.logger import get_logger
from pydefect.util.tools import (
    construct_obj_in_dict, sanitize_keys_in_dict, defaultdict_to_dict, flatten_dict)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class DefectEnergy(MSONable):
    def __init__(self,
                 defect_energy: float,
                 multiplicity: int,
                 magnetization: float,
                 convergence: bool,
                 shallow: bool):
        """

        Args:
            defect_energy (float):
                Defect formation energy at Ef = 0.
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
        self.multiplicity = multiplicity
        self.magnetization = magnetization
        self.convergence = convergence
        self.shallow = shallow


class DefectEnergies(MSONable):
    def __init__(self,
                 defect_energies: Dict[str, Dict[int, Dict[str, DefectEnergy]]],
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 title: str = None):
        """ A class related to a set of defect formation energies.

        Args:
            defect_energies (dict):
                DefectEnergy as a function of name, charge, and annotation.
                energies[name][charge][annotation] = DefectEnergy object
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
    def from_objects(cls,
                     unitcell: UnitcellCalcResults,
                     perfect: SupercellCalcResults,
                     defects: List[Defect],
                     chem_pot: tuple,
                     chem_pot_label: str,
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
            system (str):
                System name used for the title.
        """
        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute scale.
        vbm, cbm = unitcell.band_edge

        system = system if system else str(perfect.final_structure.composition)
        title = system + " condition " + chem_pot_label

        # Chemical potentials
        relative_chem_pots, standard_e = chem_pot
        relative_chem_pot = relative_chem_pots[chem_pot_label]

        # defect_energies[name][charge][annotation] = DefectEnergy
        defect_energies = \
            defaultdict(lambda: defaultdict(lambda: defaultdict(DefectEnergy)))

        for name, g in groupby(defects, key=lambda n: n.name):
            for d in g:
                # Calculate defect formation energies at the vbm
                atom_exchange_energy = 0
                for el, diff in d.changes_of_num_elements.items():
                    relative_e = relative_chem_pot.elem_coords[el]
                    standard = standard_e[el]
                    atom_exchange_energy -= diff * (relative_e + standard)

                energy = (d.relative_total_energy + d.correction_energy +
                          atom_exchange_energy)

                defect_name = DefectName(d.name, d.charge, d.annotation)
                if d.final_multiplicity.is_integer():
                    multiplicity = int(d.final_multiplicity)
                else:
                    logger.warning(
                        f"Multiplicity of {defect_name} is invalid."
                        f"initial sym: {d.initial_symmetry}, "
                        f"final sym: {d.final_symmetry}.")
                    multiplicity = d.final_multiplicity

                if not d.magnetization.is_integer() and not d.is_shallow:
                    logger.warning(
                        f"{defect_name} is not shallow but with "
                        f"fractional magnetization: {d.magnetization}")

                defect_energy = DefectEnergy(defect_energy=energy,
                                             multiplicity=multiplicity,
                                             magnetization=d.magnetization,
                                             convergence=d.is_converged,
                                             shallow=d.is_shallow)

                defect_energies[name][d.charge][d.annotation] = defect_energy

        defect_energies = defaultdict_to_dict(defect_energies)

        return cls(defect_energies=defect_energies, vbm=vbm, cbm=cbm,
                   supercell_vbm=perfect.vbm, supercell_cbm=perfect.cbm,
                   title=title)

    @classmethod
    def from_dict(cls, d):
        """ Construct a class object from a dictionary. """

        defect_energies = sanitize_keys_in_dict(d["defect_energies"])
        defect_energies = construct_obj_in_dict(defect_energies, DefectEnergy)

        return cls(defect_energies=defect_energies,
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
        for name, charge, annotation, defect_energy \
                in flatten_dict(self.defect_energies):
            outs.extend(
                [f"name: {name}:",
                 f"charge: {charge}",
                 f"Annotation: {annotation}",
                 f"Energy @ef=0 (eV): {round(defect_energy.defect_energy, 4)}",
                 f"Convergence: {defect_energy.convergence}",
                 f"Is shallow: {defect_energy.shallow}",
                 f"multiplicity: {defect_energy.multiplicity}",
                 f"magnetization: {defect_energy.magnetization}"])
            outs.append("")

        return "\n".join(outs)

    def u(self,
          name: str,
          charges: list,
          annotations: list = None) -> Tuple[float, list]:
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
            assert ValueError("The charge states {} {} {} are not incremental."
                              .format(*charges))

        energies = []
        if annotations is None:
            # annotations with lowest energies at given charge states
            annotations = []
            for c in charges:
                try:
                    annotation, min_e = \
                        min(self.defect_energies[name][c].items(),
                            key=lambda z: z[1].defect_energy)
                except KeyError:
                    print(f"Charge {c} does not exist for {name}.")
                    raise
                energies.append(min_e.defect_energy)
                annotations.append(annotation)
        else:
            for c, a in zip(charges, annotations):
                energies.append(self.defect_energies[name][c][a].defect_energy)

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
                    fermi_levels: list = None,
                    show_transition_levels: bool = False,
                    show_all_energies: bool = False,
                    color: list = None):
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
        fig, ax = plt.subplots()
        plt.title(self.title, fontsize=15)

        if color is None:
            from pydefect.core.config import COLOR
            color = COLOR

        ax.set_xlabel("Fermi level (eV)", fontsize=15)
        ax.set_ylabel("Formation energy (eV)", fontsize=15)

        # e_of_c is energy as a function of charge and annotation.
        # e_of_c[charge][annotation] = DefectEnergy object
        transition_levels = {}
        lowest_energies = defaultdict(dict)
        lowest_energy_annotations = defaultdict(dict)

        for name, energy_by_charge in self.defect_energies.items():
            for charge, energy_by_annotation in copy(energy_by_charge).items():
                for annotation, defect_energy \
                        in copy(energy_by_annotation).items():

                    n = DefectName(name, charge, annotation)
                    is_shallow = defect_energy.shallow
                    convergence = defect_energy.convergence

                    if n.is_name_matched(filtering_words) is False:
                        logger.info(f"{n} filtered out, so omitted.")
                    elif exclude_shallow_defects and is_shallow:
                        logger.info(f"{n} is shallow, so omitted.")
                    elif exclude_unconverged_defects and convergence is False:
                        logger.info(f"{n} unconverged, so omitted.")
                    else:
                        continue
                    energy_by_annotation.pop(annotation)

                # Store lowest energy and its annotation wrt name and charge.
                if energy_by_annotation:
                    annotation, min_defect_energy = \
                        min(energy_by_annotation.items(),
                            key=lambda z: z[1].defect_energy)
                    lowest_energies[name][charge] = \
                        min_defect_energy.defect_energy
                    lowest_energy_annotations[name][charge] = annotation

            # Store cross point coord and its relevant two charge states.
            cross_points_with_charges = []

            for (c1, e1), (c2, e2) \
                    in combinations(lowest_energies[name].items(), r=2):
                # The cross point between two charge states.
                x = - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                # The lowest energy among all the charge states to be compared.
                compared_energy = min([energy + c * x for c, energy
                                       in lowest_energies[name].items()])

                if y < compared_energy + 1e-5:
                    cross_points_with_charges.append([[x, y], [c1, c2]])

            # need to sort the points along x-axis.
            cpwc = sorted(cross_points_with_charges, key=lambda z: z[0][0])
            transition_levels[name] = {"cross_points": [i[0] for i in cpwc],
                                       "charges": [i[1] for i in cpwc]}

        x_min, x_max = x_range if x_range else (0, self.band_gap)
        y_min, y_max = float("inf"), -float("inf")

        def min_e_at_ef(ec, ef):
            # calculate each energy at the given Fermi level ef.
            d = {c: e + c * ef for c, e in ec.items()}
            # return the charge with the lowest energy, and its energy value
            return min(d.items(), key=lambda x: x[1])

        # Note: len(self.energies) <= len(transition_levels)
        for i, (name, tl) in enumerate(transition_levels.items()):
            # lowest_energies[name] is {} when all results are omitted above.
            if not lowest_energies[name]:
                continue
            line_type = '-' if i < 5 else '--' if i < 10 else '-.'
            # ---------- Calculate cross points including edges ---------------
            cross_points = []
            charge_set = set()

            # keep x_min -> transition levels -> x_max for connecting points.
            # Add the plot point at x_min
            charge, y = min_e_at_ef(lowest_energies[name], x_min + self.vbm)
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
                    charge_set.add(sorted(charges, key=lambda z: z)[1])
                    y_min, y_max = min([cp[1], y_min]), max([cp[1], y_max])
                    ax.scatter(cp[0] - self.vbm, cp[1], marker='o',
                               color=color[i])

            # Add the plot point at x_max
            charge, y = min_e_at_ef(lowest_energies[name], x_max + self.vbm)
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
                    s = str(round(cp[0], 2)) + ", " + str(round(cp[1], 2))
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
                annotation = lowest_energy_annotations[name][charge]
                chg = f"{charge}: {annotation}" if annotation else str(charge)

                ax.annotate(
                    chg, (x, y), color=color[i], fontsize=9,
                    bbox=dict(facecolor='white', edgecolor=color[i], pad=1))

            if show_all_energies:
                for c, e in lowest_energies[name].items():
                    y1 = e + c * (x_min + self.vbm)
                    y2 = e + c * (x_max + self.vbm)
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

        return plt

