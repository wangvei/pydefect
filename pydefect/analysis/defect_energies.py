# -*- coding: utf-8 -*-

import os

from collections import defaultdict, namedtuple

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
        self._lowest_defect_energies = {}
        for d in self._defect_energies:
            points = {}
            # Estimate the lowest energies for each defect at the vbm and cbm
            points[0.0] = min([d[i] for i in d])
            points[band_gap] = min([d[i] + band_gap * j for j in defect[i]])
            for j in defect[i]:
                for k in defect[i]:
                    if j > k:
                        x = - (defect[i][j] - defect[i][k]) / (j - k)
                        y = (j * defect[i][k] - k * defect[i][j]) / (j - k)
                        if 0 < x < band_gap and \
                                y - min([defect[i][l] + l * x for l in
                                         defect[i]]) < 0.0001:
                            points[x] = y

            outputfile.write(str(i) + "\n")
            for key in sorted(points):
                outputfile.write(str(key) + " " + str(points[key]) + "\n")
            outputfile.write("\n\n")

    def plot_energy(self, does_plot_all=False):
        self._lowest_defect_energies
