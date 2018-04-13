# -*- coding: utf-8 -*-

import numpy as np
import os

from collections import defaultdict, namedtuple
from itertools import combinations
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
        # value = {charge: energy at the vbm}

        transition_levels = {}

        for name, e_vs_c in self._defect_energies.items():
            points = {}
            charge = set()
            # Estimate the lowest energies for each defect at the vbm and cbm

            min_e_vbm = float("inf")
            for c, e in e_vs_c.items():
                if e < min_e_vbm:
                    c_vbm = c
                    min_e_vbm = e
            points[0.0] = min_e_vbm
            charge.add(c_vbm)

            min_e_cbm = float("inf")
            for c, e in e_vs_c.items():
                e_cbm = e + self._band_gap * c
                if e_cbm < min_e_cbm:
                    c_cbm = c
                    min_e_cbm = e_cbm
            points[self._band_gap] = min_e_cbm
            charge.add(c_cbm)

            points[self._band_gap] = \
                min([e + self._band_gap * c for c, e in e_vs_c.items()])

            for (c1, e1), (c2, e2) in combinations(e_vs_c.items(), r=2):
                x = - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                compared_energy = float("inf")
                for c3, e3 in e_vs_c.items():
                    if c3 == c1 or c3 == c2:
                        continue
                    else:
                        e = e3 + c3 * x
                        if e < compared_energy:
                            compared_energy = e

                if 0 < x < self._band_gap and y < compared_energy:
                    points[x] = y
                    charge.add(c1)
                    charge.add(c2)

            transition_levels[name] = [points, charge]

        return transition_levels


    def plot_energy(self, does_plot_all=False):
        self._lowest_defect_energies
