# -*- coding: utf-8 -*-
import numpy as np

from pymatgen.electronic_structure.plotter import BSPlotter, BSPlotterProjected
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.io.vasp import Kpoints, BSVasprun
from pymatgen.core.structure import Structure

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class ModBSPlotter(BSPlotter):

    # def get_plot(self, zero_to_efermi=True, ylim=None, smooth=False,
    #              vbm_cbm_marker=True, smooth_tol=None):

        # plt = super().get_plot(zero_to_efermi, ylim, smooth, vbm_cbm_marker,
        #                        smooth_tol)

        # plt.ylabel(r'$\mathrm{Energy\ (eV)}$', fontsize=30)

        # if self._bs.is_metal() is False:
        #     band_gap = self._bs.get_band_gap()["energy"]
        #     x = (plt.xlim()[0] + plt.xlim()[1]) / 2
        #     vbm = self.bs_plot_data(zero_to_efermi=zero_to_efermi)["vbm"][0][1]
        #     plt = self.add_band_gap(plt, vbm, band_gap, x)

        # return plt

    def get_plot(self, zero_to_efermi=True, ylim=None, smooth=False,
                 vbm_cbm_marker=True, smooth_tol=None):
        """
        COPIED FROM PYMATGEN.2018.5.22
        Get a matplotlib object for the bandstructure plot.
        Blue lines are up spin, red lines are down spin.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
            vbm_cbm_marker:
            smooth_tol (float) : tolerance for fitting spline to band data.
                Default is None such that no tolerance will be used.
        """
        from pymatgen.util.plotting import pretty_plot
        plt = pretty_plot(12, 8)
        import scipy.interpolate as scint

        # main internal config options
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        band_linewidth = 1

        data = self.bs_plot_data(zero_to_efermi)
        if not smooth:
            for d in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    plt.plot(data['distances'][d],
                             [data['energy'][d][str(Spin.up)][i][j]
                              for j in range(len(data['distances'][d]))], 'b-',
                             linewidth=band_linewidth)
                    if self._bs.is_spin_polarized:
                        plt.plot(data['distances'][d],
                                 [data['energy'][d][str(Spin.down)][i][j]
                                  for j in range(len(data['distances'][d]))],
                                 'r--', linewidth=band_linewidth)
        else:
            # Interpolation failure can be caused by trying to fit an entire
            # band with one spline rather than fitting with piecewise splines
            # (splines are ill-suited to fit discontinuities).
            #
            # The number of splines used to fit a band is determined by the
            # number of branches (high symmetry lines) defined in the
            # BandStructureSymmLine object (see BandStructureSymmLine._branches)

            warning = "WARNING! Distance / branch {d}, band {i} cannot be " + \
                      "interpolated.\n" + \
                      "See full warning in source.\n" + \
                      "If this is not a mistake, try increasing " + \
                      "smooth_tol.\nCurrent smooth_tol is {s}."

            for d in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    tck = scint.splrep(
                        data['distances'][d],
                        [data['energy'][d][str(Spin.up)][i][j]
                         for j in range(len(data['distances'][d]))],
                        s=smooth_tol)
                    step = (data['distances'][d][-1]
                            - data['distances'][d][0]) / 1000

                    xs = [x * step + data['distances'][d][0]
                          for x in range(1000)]

                    ys = [scint.splev(x * step + data['distances'][d][0],
                                      tck, der=0)
                          for x in range(1000)]

                    for y in ys:
                        if np.isnan(y):
                            print(warning.format(d=str(d), i=str(i),
                                                 s=str(smooth_tol)))
                            break

                    plt.plot(xs, ys, 'b-', linewidth=band_linewidth)

                    if self._bs.is_spin_polarized:
                        tck = scint.splrep(
                            data['distances'][d],
                            [data['energy'][d][str(Spin.down)][i][j]
                             for j in range(len(data['distances'][d]))],
                            s=smooth_tol)
                        step = (data['distances'][d][-1]
                                - data['distances'][d][0]) / 1000

                        xs = [x * step + data['distances'][d][0]
                              for x in range(1000)]

                        ys = [scint.splev(
                            x * step + data['distances'][d][0],
                            tck, der=0)
                            for x in range(1000)]

                        for y in ys:
                            if np.isnan(y):
                                print(warning.format(d=str(d), i=str(i),
                                                     s=str(smooth_tol)))
                                break

                        plt.plot(xs, ys, 'r--', linewidth=band_linewidth)

        self._maketicks(plt)

        # Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
        plt.ylabel(r'$\mathrm{Energy\ (eV)}$', fontsize=30)

        # Draw Fermi energy, only if not zero and metal
        if not zero_to_efermi and self._bs.is_metal():
            ef = self._bs.efermi
            plt.axhline(ef, linewidth=0.5, color='k')

        # X range (K)
        # last distance point
        x_max = data['distances'][-1][-1]
        plt.xlim(0, x_max)

        if ylim is None:
            if self._bs.is_metal():
                # Plot A Metal
                if zero_to_efermi:
                    plt.ylim(e_min, e_max)
                else:
                    plt.ylim(self._bs.efermi + e_min, self._bs.efermi + e_max)
            else:
                if vbm_cbm_marker:
                    for cbm in data['cbm']:
                        plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                    s=100)
                    for vbm in data['vbm']:
                        plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                    s=100)
                plt.ylim(data['vbm'][0][1] + e_min,
                         data['cbm'][0][1] + e_max)
        else:
            plt.ylim(ylim)
            if not self._bs.is_metal() and vbm_cbm_marker:
                for cbm in data['cbm']:
                        plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                    s=100)
                for vbm in data['vbm']:
                        plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                    s=100)

        plt.tight_layout()

        if self._bs.is_metal() is False:
            band_gap = self._bs.get_band_gap()["energy"]
            x = (plt.xlim()[0] + plt.xlim()[1]) / 2
            vbm = self.bs_plot_data(zero_to_efermi=zero_to_efermi)["vbm"][0][1]
            plt = self.add_band_gap(plt, vbm, band_gap, x)

        return plt

    def add_band_gap(self, plt, vbm, band_gap, annotate_x, line_color="purple"):
        plt.hlines([vbm, vbm + band_gap], plt.xlim()[0], plt.xlim()[1],
                   line_color, linestyles='dashed')

        plt.annotate('', xy=(annotate_x, vbm),
                     xytext=(annotate_x, vbm + band_gap),
                     fontsize=7, color='black',
                     arrowprops=dict(edgecolor='black', arrowstyle='<|-|>',
                                     shrinkA=0, shrinkB=0))
        # TODO: When adding the following annotation, the Gamma notation
        # TODO: disappears at x=0.
        plt.annotate(str(round(vbm + band_gap, 2)) + " eV",
                     (annotate_x * 1.03, (vbm + band_gap) / 2), fontsize=15)

        return plt

    def plot_compare(self, other_plotter, ylim=None, legend=True,
                     zero_to_efermi=True):
        """
        Args:
            other_plotter: another BSPlotter class object.
            legend: If show figure legend or not
            zero_to_efermi: If the energy is aligned or not.
        """

        plt = self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi)
        data_orig = self.bs_plot_data(zero_to_efermi=zero_to_efermi)
        data = other_plotter.bs_plot_data(zero_to_efermi=zero_to_efermi)

        second_band_linewidth = 2
        for i in range(other_plotter._nb_bands):
            for d in range(len(data_orig['distances'])):
                plt.plot(data_orig['distances'][d],
                         [e[str(Spin.up)][i] for e in data['energy']][d],
                         'c-', linewidth=second_band_linewidth)
                if other_plotter._bs.is_spin_polarized:
                    plt.plot(data_orig['distances'][d],
                             [e[str(Spin.down)][i] for e in data['energy']][d],
                             'm--', linewidth=second_band_linewidth)

        if other_plotter._bs.is_metal() is False:
            band_gap = other_plotter._bs.get_band_gap()["energy"]
            x = (plt.xlim()[0] * 0.4  + plt.xlim()[1] * 0.6)
            vbm = other_plotter.bs_plot_data(
                zero_to_efermi=zero_to_efermi)["vbm"][0][1]
            plt = self.add_band_gap(plt, vbm, band_gap, x)

        import matplotlib.lines as mlines
        if legend:
            if self._bs.is_spin_polarized:
                handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                         label='bs 1 up'),
                           mlines.Line2D([], [], linewidth=2, color='r',
                                         label='bs 1 down', linestyle="--")]
            else:
                handles = [mlines.Line2D([], [], linewidth=2, color='b',
                                         label='bs 1')]
            if other_plotter._bs.is_spin_polarized:
                handles.extend([mlines.Line2D([], [], linewidth=2, color='c',
                                              label='bs 2 up'),
                                mlines.Line2D([], [], linewidth=2, color='m',
                                              linestyle="--",
                                              label='bs 2 down')])
            else:
                handles.extend([mlines.Line2D([], [], linewidth=2, color='c',
                                              label='bs 2')])
            plt.legend(handles=handles)

        return plt


class VaspBandStructureSymmLine(BandStructureSymmLine):

    def __init__(self, kpoints_name, vasprun_name, is_projection=False,
                 poscar_name=None):

        kpoints = Kpoints.from_file(kpoints_name)

        if is_projection:
            vasprun = BSVasprun(filename=vasprun_name,
                                parse_projected_eigen=True)
        else:
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

        if is_projection:
            v = vasprun.projected_eigenvalues
            projections = {}
            for s in eigenvalues:
                # swapaxes is essential
                # When parsing vasprun.xml using BSVasprun,
                #  Vasprun.projected_eigenvalues[spin][kpoint index][band index]
                #                               [atom index][orbital_index]
                # For BSPlotter
                # projections: dict of orbital projections as {spin: ndarray}.
                # The indices of the ndarrayare [band_index, kpoint_index,
                #                                orbital_index, ion_index].

                projections[Spin(s)] = v[Spin(s)].swapaxes(0, 1).swapaxes(2, 3)
                print(projections)
            structure = Structure.from_file(poscar_name)
            super().__init__(kpts_wo_weight, eigenvalues_wo_weight, lattice_rec,
                             efermi, labels_dict, structure, projections)

        else:
            super().__init__(kpts_wo_weight, eigenvalues_wo_weight, lattice_rec,
                             efermi, labels_dict)


