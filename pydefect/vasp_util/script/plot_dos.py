# -*- coding: utf-8 -*-

import numpy as np

from pymatgen.electronic_structure.plotter import DosPlotter
from collections import OrderedDict
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import Vasprun

from pydefect.input_maker.defect_initial_setting import SYMPREC

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class ModDosPlotter(DosPlotter):

    def get_plot(self, xlim=None, ylims=None, cbm_vbm=None, legend=True):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylims: Specifies the y-axes limits. Two types of input.
                [[y1min, y1max], [y2min, y2max], ..]
            legend:
        """

        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)
        import palettable
        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        y = None
        all_densities = []
        all_energies = []

        for key, dos in self._doses.items():
            energies = dos['energies']
            densities = dos['densities']
            if not y:
                y = {Spin.up: np.zeros(energies.shape),
                     Spin.down: np.zeros(energies.shape)}
            all_energies.append(energies)
            all_densities.append(densities)

        # Make groups to be shown in the same figure.
        keys = list(self._doses.keys())
        grouped_keys = OrderedDict()
        for k in keys:
            first_word = k.split()[0]
            if first_word in grouped_keys:
                grouped_keys[first_word].append(k)
            else:
                grouped_keys[first_word] = [k]

        import matplotlib.pyplot as plt
        num_figs = len(grouped_keys)
        fig, axs = plt.subplots(num_figs, 1, sharex=True)

        if xlim:
            axs[0].set_xlim(xlim)

        n = 0
        for i, gk in enumerate(grouped_keys):
            all_pts = []
            for j, key in enumerate(grouped_keys[gk]):
                x = []
                y = []
                for spin in [Spin.up, Spin.down]:
                    if spin in all_densities[n]:
                        densities = list(int(spin) * all_densities[n][spin])
                        energies = list(all_energies[n])
                        x.extend(energies)
                        y.extend(densities)
                all_pts.extend(list(zip(x, y)))
                axs[i].plot(x, y, color=colors[j % ncolors], label=str(key),
                            linewidth=2)
                n += 1

            # plot vertical lines for band edges or Fermi level
            if self.zero_at_efermi:
                # plot a line
                axs[i].axvline(0, color="black", linestyle="--", linewidth=0.5)
                if cbm_vbm:
                    axs[i].axvline(cbm_vbm[0] - cbm_vbm[1], color="black",
                                   linestyle="--", linewidth=0.5)
            else:
                axs[i].axvline(self._doses[key]['efermi'],
                               color="black", linestyle="--", linewidth=0.5)
                if cbm_vbm:
                    axs[i].axvline(cbm_vbm[0], color="black", linestyle="--",
                                   linewidth=0.5)

            if legend:
                axs[i].legend(loc="best", markerscale=0.1)
#                axs[i].legend(bbox_to_anchor=(1.1, 0.8), loc="best")
                leg = axs[i].get_legend()
                for legobj in leg.legendHandles:
                    legobj.set_linewidth(1.2)
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize=7)
            else:
                axs[i].set_title(key, fontsize=7)

            axs[i].axhline(0, color="black", linewidth=0.5)

        if ylims and len(ylims) not in (num_figs, 2):
            raise ValueError("The number of y-ranges is not proper.")

        if ylims and len(ylims) == 2:
            axs[0].set_ylim(ylims[0])
            for i in range(1, len(axs)):
                axs[i].set_ylim(ylims[1])
        elif ylims:
            for i in range(len(axs)):
                axs[i].set_ylim(ylims[i])
        # else:
        #     for i in range(len(axs)):
        #         ylim = axs[i].get_ylim()
        #         print(ylim)
        #         relevanty = [p[1] for p in all_pts
        #                      if ylim[0] < p[0] < ylim[1]]
        #         axs[i].set_ylim((min(relevanty), max(relevanty)))

        axs[-1].set_xlabel('Energy (eV)')
        plt.tight_layout()

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                            wspace=0.2, hspace=0.2)

        return plt


def get_dos_plot(vasprun_file, sites=None, orbital=True, xlim=None, ymaxs=None,
                 zero_at_efermi=True, legend=True, symprec=SYMPREC):
    v = Vasprun(vasprun_file, ionic_step_skip=True, parse_eigen=False)
    complete_dos = v.complete_dos

    # check cbm

    if complete_dos.get_gap() > 0.1:
        cbm_vbm = complete_dos.get_cbm_vbm()
    else:
        cbm_vbm = None

    structure = v.final_structure

    dos = OrderedDict()
    dos["Total"] = complete_dos

    if sites is None:
        s = v.final_structure
        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        symmetrized_structure = \
            SpacegroupAnalyzer(s, symprec=symprec).get_symmetrized_structure()

        equiv_indices = symmetrized_structure.equivalent_indices
        sites = [indices[0] for indices in equiv_indices]
    else:
        # increment site index
        sites = [i - 1 for i in sites]

    if True:
        orbital_set = ["s", "p", "d", "f"]
        for s in sites:
            site = structure[s]
            if orbital:
                for i, o in enumerate(orbital_set):
                    name = "Site:" + str(s + 1) + " " + site.specie.symbol \
                           + "-" + o
                    try:
                        dos[name] = \
                            complete_dos.get_site_spd_dos(site)[OrbitalType[o]]
                    except:
                        # print("{} orbital does not exist.".format(o))
                        del orbital_set[i]
            else:
                dos["Site:" + str(s + 1) + " " + site.specie.symbol] = \
                    complete_dos.get_site_dos(site)

    plotter = ModDosPlotter(zero_at_efermi=zero_at_efermi)
    plotter.add_dos_dict(dos)

    if "ISPIN" in v.incar and v.incar["ISPIN"] == 2:
        spin = True
    else:
        spin = False

    if ymaxs:
        if spin:
            ylims = [[-y, y] for y in ymaxs]
        else:
            ylims = [[0, y] for y in ymaxs]
    else:
        ylims = None

    return plotter.get_plot(xlim=xlim, ylims=ylims, cbm_vbm=cbm_vbm,
                            legend=legend)
