# -*- coding: utf-8 -*-
import numpy as np

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.io.vasp import Kpoints, BSVasprun

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class ModBSPlotter(BSPlotter):
    def add_effective_masses(self):
        pass

    def add_band_gap(self):
        pass

    def show(self, zero_to_efermi=True, ylim=None, smooth=False,
             vbm_cbm_marker=True, smooth_tol=None):
        """
        Overwrite to support vbm_cbm_marker flag.
        Show the plot using matplotlib.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
            smooth_tol (float) : tolerance for fitting spline to band data.
                Default is None such that no tolerance will be used.
        """
        plt = self.get_plot(zero_to_efermi, ylim, smooth, vbm_cbm_marker)
        plt.show()

    def save_plot(self, filename, img_format="eps", ylim=None,
                  zero_to_efermi=True, smooth=False, vbm_cbm_marker=True):
        """
        Overwrite to support vbm_cbm_marker flag.
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi,
                            smooth=smooth, vbm_cbm_marker=vbm_cbm_marker)
        plt.savefig(filename, format=img_format)
        plt.close()


def plot_band_structure(kpoints_name, vasprun_name):

    kpoints = Kpoints.from_file(kpoints_name)
    vasprun = BSVasprun(filename=vasprun_name)

    eigenvalues = vasprun.eigenvalues
    lattice_rec = vasprun.lattice_rec
    efermi = vasprun.efermi

    first_index_wo_weight = 0
    for w in kpoints.kpts_weights:
        if w > 0:
            first_index_wo_weight += 1
            continue

    kpts_wo_weight = kpoints.kpts[first_index_wo_weight:]

    eigenvalues_wo_weight = {}
    for s in eigenvalues:
        # transpose is essential
        # When parsing vasprun.xml using BSVasprun,
        # eigenvalues[kpt-index][band-index] = [energy, occupation]
        # For BSPlotter
        # eigenvalues[band-index][kpt-index] = energy
        eigenvalues_wo_weight[Spin(s)] = \
            eigenvalues[Spin(s)][first_index_wo_weight:, :, 0].transpose()

    # Store label except for "None".
    labels_dict = {}
    for i, label in enumerate(kpoints.labels):
        if label != "None" and label is not None:
            # TODO: Add more greek letters
            if label == "GAMMA":
                labels_dict["\u0393"] = kpoints.kpts[i]
            else:
                labels_dict[label] = kpoints.kpts[i]

    band_structure = BandStructureSymmLine(kpts_wo_weight,
                                           eigenvalues_wo_weight,
                                           lattice_rec, efermi, labels_dict)

    bs_plotter = ModBSPlotter(band_structure)
    bs_plotter.show()


