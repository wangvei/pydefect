# -*- coding: utf-8 -*-

import numpy as np

from pymatgen.electronic_structure.plotter import DosPlotter
from collections import OrderedDict
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure
from pymatgen.electronic_structure.core import Spin
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

    def get_plot(self, xlim=None, ylim=None):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
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


        from pymatgen.util.plotting import pretty_plot
        plt = pretty_plot(12, 8)

        keys = list(self._doses.keys())
        all_pts = []
        for i, key in enumerate(keys):
            x = []
            y = []
            for spin in [Spin.up, Spin.down]:
                if spin in all_densities[i]:
                    densities = list(int(spin) * all_densities[i][spin])
                    energies = list(all_energies[i])
                    x.extend(energies)
                    y.extend(densities)
            all_pts.extend(list(zip(x, y)))
            plt.plot(x, y, color=colors[i % ncolors],
                     label=str(key), linewidth=3)
            # plot zeros for all the DOS
            if not self.zero_at_efermi:
                ylim = plt.ylim()
                plt.plot([self._doses[key]['efermi'],
                          self._doses[key]['efermi']], ylim,
                         color=colors[i % ncolors],
                         linestyle='--', linewidth=2)

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in all_pts
                         if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        plt.axhline(0, color="black", linewidth=0.5)
        if self.zero_at_efermi:
            # plot a line
            plt.axvline(0, color="black", linewidth=0.5)

        plt.xlabel('Energy (eV)')
        plt.ylabel('Density of states')

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text. Text instance in the legend
        plt.setp(ltext, fontsize=20)
        plt.tight_layout()

        return plt


def get_dos_plot(vasprun_file, sites=False, element=False, orbital=True,
                 symprec=SYMPREC):
    v = Vasprun(vasprun_file)
    complete_dos = v.complete_dos

    structure = v.final_structure

    dos = OrderedDict()
    dos["Total"] = complete_dos

    if sites:
        for s in sites:
            site = structure[s]
            dos["Site " + str(s) +" " + site.specie.symbol] = \
                complete_dos.get_site_dos(site)

    plotter = ModDosPlotter(zero_at_efermi=False)
    plotter.add_dos_dict(dos)

    return plotter.get_plot()



    # if site:
    #     s = v.final_structure
    #     # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
    #     symmetrized_structure = \
    #         SpacegroupAnalyzer(s, symprec=symprec).get_symmetrized_structure()
    #     equiv_sites = symmetrized_structure.equivalent_sites
    #     print(equiv_sites)
    #     # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
    #     from collections import defaultdict
    #     num_irreducible_sites = defaultdict(int)

        # # irreducible_sites (list): a set of IrreducibleSite class objects
        # irreducible_sites = []

        # last_index = 0

        # for i, equiv_site in enumerate(equiv_sites):
        #     # set element name of equivalent site
        #     element = equiv_site[0].species_string

            # # increment number of inequivalent sites for element
            # num_irreducible_sites[element] += 1

            # # the following np.array type must be converted to list
            # # to keep the consistency of the IrreducibleSite object.
            # irreducible_name = element + str(num_irreducible_sites[element])

            # all_dos["Site " + str(i) + " " + site.specie.symbol] = \
            #     complete_dos.get_site_dos(site)


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
