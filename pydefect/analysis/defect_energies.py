# -*- coding: utf-8 -*-

import numpy as np
import os

from collections import defaultdict, namedtuple
from scipy.spatial import HalfspaceIntersection

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults


class DefectEnergies:

    def __init__(self, unitcell, perfect, defects, chem_pot, chem_pot_label):
#        is_lower_energy=False):
        """
        Args:
        Defect = namedtuple("Defect", "defect_entry", "dft_results", "correction")
            unitcell (UnitcellDftResults):
            chem_pot (ChemPot):
            perfect (SupercellDftResults)
            defects (list of namedtuple Defect):
            [[DefectEntry, SupercellDftResults, Correction], ...]
        """

        self._vbm, self._cbm = unitcell.band_edge
        self._band_gap = self._cbm - self._vbm

        # TODO: check if exists
#        self._vbm2, self._cbm2 = unitcell.band_edge2

        chem_pot = chem_pot[chem_pot_label]

        self._defect_energies = defaultdict(lambda: defaultdict(float))
        for d in defects:
            name = d.defect_entry.name

            charge = d.defect_entry.charge
            element_diff = d.defect_entry.element_diff

            relative_energy = d.dft_results.relative_total_energy(perfect)
            correction_energy = d.correction.total_correction_energy
            electron_interchange_energy = self._vbm * charge
            element_interchange_energy = \
                - sum([v * chem_pot[k] for k, v in element_diff.items()])

            energy = relative_energy + correction_energy + \
                     electron_interchange_energy + element_interchange_energy

            self._defect_energies[name][charge] = energy

#        self._defect_energies = dict(defect_energies)

#            if self._defect_energies[name][charge]:
#                if is_lower_energy is True:
#                    if self._defect_energies[name][charge] > energy:
#                        self._defect_energies[name][charge] = energy
#                else:
#                    raise ArithmeticError


    # @classmethod
    # def from_directories(cls, unitcell_dir, chem_pot_dir, perfect_dir, defect_dirs):
    #     """
    #     Args:
    #     """

        # unitcell = UnitcellDftResults.\
        #     json_load(os.path.join(unitcell_dir, "unitcell.json"))
        # chem_pot = "xx"
        # perfect = SupercellDftResults.json_load(os.path.join(perfect_dir, "perfect.json"))

        # defects = []
        # for d_dir in defect_dirs:
        #     defect_entry = DefectEntry.json_load(os.path.join(d_dir, "defect_entry.json"))
        #     dft_results = SupercellDftResults.json_load(os.path.join(d_dir, "defect_entry.json"))
        #     correction = "aaa"

            # defects.append([defect_entry, dft_results, correction])

        # return cls(unitcell, chem_pot, perfect, defects)


    def calc_transition_levels(self):
        points = {}
        # value = {charge: energy at the vbm}
        for name, value in self._defect_energies.items():
            half_spaces = [[- charge, 1, - energy_at_vbm]
                           for charge, energy_at_vbm in value.items()]
            # Estimate the energies at the vbm and cbm
            energies = np.array(
                [[energy_at_vbm, energy_at_vbm + charge * self._band_gap]
                 for charge, energy_at_vbm in value.items()])

            min_energy = np.amin(energies)
            small_number = 0.1
            lowest_energy = min_energy - small_number

            half_spaces.append([-1, 0, 0]) # x >= 0
            half_spaces.append([1, 0, - self._band_gap]) # x <= band_gap
            half_spaces.append([0, -1, lowest_energy]) # y >= lowest_energy
            half_spaces = np.array(half_spaces)

            feasible_point = np.array([small_number, min_energy])
            hs = HalfspaceIntersection(half_spaces, feasible_point)

            points_by_name = []
            # hs.dual_facets[0]: two half_spaces comprising 1st crossing point.
            #                    e.g., [0, 1] --> 0th and 1st lines
            # hs.dual_facets[1]: two half_spaces comprising 2nd crossing point.
            for i, comprising_lines in enumerate(hs.dual_facets):

                if comprising_lines[0] < len(half_spaces) - 3 and \
                        comprising_lines[1] < len(half_spaces) - 3:
                    charges = [int(-half_spaces[comprising_lines[0]][0]),
                               int(-half_spaces[comprising_lines[1]][0])]

                elif comprising_lines[0] == len(half_spaces) - 1:
                    charges = [int(-half_spaces[comprising_lines[0]][0])]

                elif comprising_lines[1] == len(half_spaces) - 1:
                    charges = [int(-half_spaces[comprising_lines[1]][0])]

                else:
                    break
                points_by_name.append([hs.intersections[i], charges])

            points[name] = np.array(points_by_name)

        return points


    def plot_energy(self, does_plot_all=False):
        self._lowest_defect_energies
