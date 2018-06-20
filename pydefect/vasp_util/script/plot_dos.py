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

    def get_plot(self, xlim=None, ylims=None):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylims: Specifies the y-axes limits.
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

        keys = list(self._doses.keys())
        grouped_keys = {}
        for k in keys:
            first_word = k.split()[0]
            if first_word in grouped_keys:
                grouped_keys[first_word].append(k)
            else:
                grouped_keys[first_word] = [k]

        import matplotlib.pyplot as plt
        num_figs = len(grouped_keys)
        fig, axs = plt.subplots(num_figs, 1, sharex=True)
        fig.subplots_adjust(hspace=0)

        if ylims and num_figs != len(ylims):
            assert ValueError("The number of y-ranges is not proper.")

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

            # plot zeros for all the DOS
            if self.zero_at_efermi:
                # plot a line
                axs[i].axvline(0, color="black", linestyle="--", linewidth=0.5)
            else:
                axs[i].axvline(self._doses[key]['efermi'],
                               color="black", linestyle="--", linewidth=0.5)
                # ylim = axs[i].get_ylim()
                # axs[i].plot([self._doses[key]['efermi'],
                #              self._doses[key]['efermi']], ylim,
                #             color=colors[j % ncolors], linestyle='--',
                #             linewidth=2)

            # if xlim:
            #     axs[i].set_xlim(xlim)
            # if ylim:
            #     axs[i].set_ylim(ylim)
            # else:
            #     xlim = axs[i].get_ylim()
            #     relevanty = [p[1] for p in all_pts
            #                  if xlim[0] < p[0] < xlim[1]]
            #     axs[i].set_ylim((min(relevanty), max(relevanty)))

            axs[i].axhline(0, color="black", linewidth=0.5)

    #        axs[i].set_title(key)

            axs[i].legend()
#            axs[i].legend(bbox_to_anchor=(1.1, 0.8), loc="left")
    #            leg = axs[i].get_legend()
    #            ltext = leg.get_texts()  # all the text. Text instance in the legend
    #            axs[i].setp(ltext, fontsize=20)

        axs[-1].set_xlabel('Energy (eV)')
        # plt.tight_layout()
#        plt.subplots_adjust(right=0.8)

        return plt


def get_dos_plot(vasprun_file, sites=False, orbital=True, symprec=SYMPREC):
    v = Vasprun(vasprun_file)
    complete_dos = v.complete_dos

    structure = v.final_structure

    dos = OrderedDict()
    dos["Total"] = complete_dos

    if sites:
       pass
    else:
        s = v.final_structure
        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        symmetrized_structure = \
            SpacegroupAnalyzer(s, symprec=symprec).get_symmetrized_structure()

        print(symmetrized_structure.site_labels)
        print(symmetrized_structure.equivalent_sites)
        print(symmetrized_structure.equivalent_indices)
        print(symmetrized_structure.wyckoff_symbols)

        equiv_indices = symmetrized_structure.equivalent_sites
        equiv_sites = symmetrized_structure.equivalent_sites

        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        from collections import defaultdict
        num_irreducible_sites = defaultdict(int)

        # irreducible_sites (list): a set of IrreducibleSite class objects
        last_index = 0

        for i, (index, site) in enumerate(zip(equiv_indices, equiv_sites)):
            # set element name of equivalent site
            element = site.species_string

            # increment number of inequivalent sites for element
            num_irreducible_sites[element] += 1

            # the following np.array type must be converted to list
            # to keep the consistency of the IrreducibleSite object.
            irreducible_name = element + str(num_irreducible_sites[element])

            dos =


        all_dos["Site " + str(i) + " " + site.specie.symbol] = \
            complete_dos.get_site_dos(site)


        # for i in range(len(s)):
        #     site = s[i]

        # if element:
        #     syms = [tok.strip() for tok in element[0].split(",")]
        #     all_dos = {}
        #     for el, dos in complete_dos.get_element_dos().items():
        #         if el.symbol in syms:
        #             all_dos[el] = dos
        # if orbital:
        #     all_dos = complete_dos.get_spd_dos()







    # if True:
    #     for s in sites:
    #         site = structure[s]
    #         print(complete_dos.get_site_spd_dos(site)[OrbitalType.s])

#             if orbital:
#                 orbital_set = ["s", "p", "d", "f"]
# #                if forbital:
# #                    orbital_set = ["s", "p", "d", "f"]
# #                else:
# #                    orbital_set = ["s", "p", "d"]

            #     for o in orbital_set:
            #         name = "Site:" + str(s) + " " + site.specie.symbol + "-" + o
            #         try:
            #             dos[name] = \
            #                 complete_dos.get_site_spd_dos(site)[OrbitalType[o]]
            #         except:
            #             print("{} does not exist.".format(name))
            # else:
            #     dos["Site:" + str(s) + " " + site.specie.symbol] = \
            #         complete_dos.get_site_dos(site)

    # plotter = ModDosPlotter(zero_at_efermi=False)
    # plotter.add_dos_dict(dos)

    # return plotter.get_plot()



