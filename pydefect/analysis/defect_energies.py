# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.path as path

from collections import defaultdict, namedtuple
from itertools import combinations
from scipy.spatial import HalfspaceIntersection

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults


class DefectEnergies:
    """
    A class related to a set of defect formation energies.
    """
    def __init__(self, unitcell, perfect, defects, chem_pot, chem_pot_label):
        """
        Calculates defect formation energies.
        Args:

            unitcell (UnitcellDftResults):
            perfect (SupercellDftResults)
            defects (list of namedtuple Defect):
                [[Defect, ...]
                Defect = namedtuple("Defect", "defect_entry", "dft_results",
                                "correction")
            chem_pot (ChemPot):
            chem_pot_label (str):
        """

        self._vbm, self._cbm = unitcell.band_edge
        # TODO: check if exists
        # self._vbm2, self._cbm2 = unitcell.band_edge2
        self._transition_levels = {}

        chem_pot = chem_pot[chem_pot_label]

        self._defect_energies = defaultdict(dict)

        for d in defects:
            name = d.defect_entry.name
            charge = d.defect_entry.charge
            element_diff = d.defect_entry.element_diff

            # calculate four terms for a defect formation energy.
            relative_energy = d.dft_results.relative_total_energy(perfect)
            correction_energy = d.correction.total_correction_energy
            electron_interchange_energy = self._vbm * charge
            element_interchange_energy = \
                - sum([v * chem_pot[k] for k, v in element_diff.items()])

            self._defect_energies[name][charge] = \
                relative_energy + correction_energy + \
                electron_interchange_energy + element_interchange_energy

    @staticmethod
    def min_e_at_ef(ec, ef):
        d = {c: e + c * ef for c, e in ec.items()}
        return min(d.items(), key=lambda x: x[1])

    def calc_transition_levels(self):
        """"""
        transition_levels = {}
        TransitionLevel = namedtuple("TransitionLevel",
                                     ("cross_points", "charges"))

        for name, e_of_c in self._defect_energies.items():
            points = []
            charge = set()

            # Estimate the lowest energies for each defect at the vbm
            c_vbm, min_e_vbm = self.min_e_at_ef(e_of_c, 0)
            points.append([0.0, min_e_vbm])
            charge.add(c_vbm)

            # Estimate the lowest energies for each defect at the cbm
            c_cbm, min_e_cbm = self.min_e_at_ef(e_of_c, self.band_gap)
            points.append([self.band_gap, min_e_cbm])
            charge.add(c_cbm)

            for (c1, e1), (c2, e2) in combinations(e_of_c.items(), r=2):
                # Estimate the cross point between two charge states
                x= - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                compared_energy = self.min_e_at_ef(e_of_c, x)[1] + 0.00001

                if 0 < x < self.band_gap and y < compared_energy:
                    points.append([x, y])
                    charge.add(c1)
                    charge.add(c2)

            transition_levels[name] = \
                TransitionLevel(cross_points=sorted(points, key=lambda x: x[0]),
                                charges=sorted(charge))

        self._transition_levels = transition_levels

    def plot_energy(self):
        fig = plt.figure()
#        ax = fig.add_subplot(111)
#        y_max =
#        y_min = np.array(tl.cross_points).transpose()

        for name, tl in self.transition_levels.items():
            x, y = np.array(tl.cross_points).transpose()
            plt.plot(x, y, 'ro-')

        plt.show()

#        ax.plot()

    @property
    def transition_levels(self):
        """
        transition_levels[name] = [[E_F, energy], [charge1, charge2]]

        """
        return self._transition_levels

    @property
    def vbm(self):
        return self._vbm

    @property
    def cbm(self):
        return self._cbm

    @property
    def band_gap(self):
        return self._cbm - self._vbm



