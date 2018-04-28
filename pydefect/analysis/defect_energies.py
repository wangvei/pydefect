# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings

from collections import defaultdict, namedtuple
from itertools import combinations

from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.input_maker.input_maker import is_name_selected


Defect = namedtuple("Defect", ("defect_entry", "dft_results", "correction"))
TransitionLevel = namedtuple("TransitionLevel", ("cross_points", "charges"))


class DefectEnergies:
    """
    A class related to a set of defect formation energies.
    """
    def __init__(self, unitcell, perfect, defects, chem_pot, chem_pot_label,
                 filtering_words=None, system_name=""):
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
        self.title = system_name + " condition " + chem_pot_label
        # TODO: check if exists
        # self._vbm2, self._cbm2 = unitcell.band_edge2
        self._transition_levels = {}

        chem_pot = chem_pot[chem_pot_label]

        self._defect_energies = defaultdict(dict)

        for d in defects:
            name = d.defect_entry.name
            if filtering_words \
                    and is_name_selected(name, filtering_words) is False:
                continue
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

    def print_energies(self):
        pass

    @staticmethod
    def min_e_at_ef(ec, ef):
        d = {c: e + c * ef for c, e in ec.items()}
        return min(d.items(), key=lambda x: x[1])

    def calc_transition_levels(self):
        """
        Calculates a set of transition levels

        Parameter in use:
            TransitionLevel (namedtuple):
        """

        transition_levels = {}

        for name, e_of_c in self._defect_energies.items():
            points = []
            charge = set()

            # Estimate the lowest energies for each defect at the vbm
            c_vbm, min_e_vbm = self.min_e_at_ef(e_of_c, 0)
            points.append([0.0, min_e_vbm])
            charge.add(c_vbm)

            # Same but at the cbm
            c_cbm, min_e_cbm = self.min_e_at_ef(e_of_c, self.band_gap)
            points.append([self.band_gap, min_e_cbm])
            charge.add(c_cbm)

            for (c1, e1), (c2, e2) in combinations(e_of_c.items(), r=2):
                # Estimate the cross point between two charge states
                x = - (e1 - e2) / (c1 - c2)
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

    def plot_energy(self, file_name=None, x_range=None, y_range=None,
                    show_tls=False):
        """
        Plots the defect formation energies as a function of the Fermi level.

        Args:
            file_name (str): File name for saving the plot.
            x_range (2x1 list): x range for the plot.
            y_range (2x1 list): y range for the plot.
            show_tls (bool): If transition level values are shown or not.
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)

        plt.title(self.title)

        ax.set_xlabel("Fermi level (eV)")
        ax.set_ylabel("Formation energy (eV)")

        x_max = max([max(np.array(tl.cross_points).transpose()[0])
                     for tl in self.transition_levels.values()])
        x_min = 0
        y_max = max([max(np.array(tl.cross_points).transpose()[1])
                     for tl in self.transition_levels.values()])
        y_min = min([min(np.array(tl.cross_points).transpose()[1])
                     for tl in self.transition_levels.values()])
        if y_min > 0:
            y_min = 0

        if x_range:
            plt.xlim(x_range[0], x_range[1])
        if y_range:
            plt.ylim(y_range[0], y_range[1])
        else:
            margin = (y_max - y_min) * 0.08
            plt.ylim(y_min - margin, y_max + margin)

        margin_plot_x = (x_max - x_min) * 0.025
        margin_plot_y = (y_max - y_min) * 0.025

        # support lines
        plt.axvline(x=0, linewidth=1.0, linestyle='dashed')
        plt.axvline(x=x_max, linewidth=1.0, linestyle='dashed')
        plt.axhline(y=0, linewidth=1.0, linestyle='dashed')

        for i, (name, tl) in enumerate(self.transition_levels.items()):

            color = matplotlib.cm.hot(float(i) / len(self.transition_levels))

            cross_points = tl.cross_points
            shallow = []
            transition_levels = []

            for cp in cross_points:
                if cp[0] == 0.0 or cp[0] == self.band_gap:
                    if cp[0] == 0.0 and max(tl.charges) < 0:
                        shallow.append(cp)
                    elif cp[0] == self.band_gap and min(tl.charges) > 0:
                        shallow.append(cp)
                else:
                    transition_levels.append(cp)

            # set the x and y arrays to be compatible with matplotlib style.
            # e.g., x = [0.0, 0.3, 0.5, ...], y = [2.1, 3.2, 1.2, ...]
            # plot lines
            x, y = np.array(cross_points).transpose()
            ax.plot(x, y, '-', color=color, label=name)

            # plot filled circles for transition levels.
            if transition_levels:
                x, y = np.array(transition_levels).transpose()
                ax.scatter(x, y, color=color)

            # plot unfilled circles for shallow levels.
            if shallow:
                x, y = np.array(shallow).transpose()
                ax.scatter(x, y, facecolor="white", edgecolor=color)

            if show_tls:
                for cp in cross_points:
                    s = str(round(cp[0], 2)) + ", " + str(round(cp[1], 2))
                    ax.annotate(s,
                                (cp[0] + margin_plot_x, cp[1] - margin_plot_y),
                                color=color, fontsize=8)

            ax.legend(bbox_to_anchor=(1.05, 0.0), loc="lower left")

            # Arrange the charge states at the middle of the transition levels.
            middle_points = \
                reversed([[(a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + margin_plot_y]
                          for a, b in zip(cross_points, cross_points[1:])])

            for j, (x, y) in enumerate(middle_points):
                ax.annotate(str(tl.charges[j]), (x, y), color=color)

        fig.subplots_adjust(right=0.75)

        if file_name:
            plt.savefig(file_name)
        else:
            plt.show()

    @property
    def defect_energies(self):
        return self._defect_energies

    @property
    def transition_levels(self):
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


def main():
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--perfect", dest="perfect",
                        default="perfect", type=str)
    parser.add_argument("--defects", dest="defects",
                        default="", type=str)

    opts = parser.parse_args()

    perfect = SupercellDftResults.load_json(opts.perfect + "dft_results.json")

    if not opts.defects:
        from glob import glob
        defects_dirs = glob('*[0-9]/')
    else:
        defects_dirs = opts.defects

    defects = []
    for d in defects_dirs:
        try:
            defect_dft_results = \
                SupercellDftResults.load_json(d + "dft_results.json")
            defects.append(defect_dft_results)
        except:
            warnings.warn(message="Parsing data in " + d + " is failed.")



if __name__ == "__main__":
    main()
