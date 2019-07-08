# -*- coding: utf-8 -*-

import json
import numpy as np

from monty.json import MSONable, MontyEncoder

import matplotlib.pyplot as plt

from pymatgen.electronic_structure.core import Spin

from pydefect.analysis.defect import Defect
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.error_classes import UnitcellCalcResultsError
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class DefectEigenvalue(MSONable):
    """ A class related to eigenvalues in a defect calculation. """

    def __init__(self,
                 name: str,
                 charge: int,
                 kpoint_coords: list,
                 perfect_kpoint_coords: list,
                 kpoint_weights: list,
                 eigenvalues: np.array,
                 perfect_eigenvalues: np.array,
                 perfect_symmops: list,
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 fermi_level: float,
                 total_magnetization: float,
                 orbital_character: dict = None,
                 eigenvalue_correction: dict = None,
                 band_edges: dict = None):
        """
        Args:
            name (str):
                Name of a defect
            charge (int):
                Defect charge.
            kpoint_coords (list):
                List of k-point coordinates
            kpoint_weights (list):
                List of k-point weights.
            eigenvalues (N_spin x N_kpoint x N_band np.array):
                Numpy array for the electron eigenvalues in absolute scale.
                e.g., eigenvalues[Spin.up][0][0] = array([-8.3171,  1.    ])
                                                           energy  occupation
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            supercell_vbm (float):
                Valence band maximum in the perfect supercell.
            supercell_cbm (float):
                Conduction band minimum in the perfect supercell.
            fermi_level (float):
                Fermi level.
            total_magnetization (float):
                Total total_magnetization.
            eigenvalue_correction (dict):
                Dict with key of correction name and value of correction value.
            band_edges (dict):
                Band edge states at each spin channel.
                None: no in gap state.
                "Donor PHS": Donor-type perturbed host state (PHS).
                "Acceptor PHS": Acceptor-type PHS.
                "Localized state": With in-gap localized state.
                    ex. {Spin.up: None, Spin.down:"Localized state}
        """
        self.name = name
        self.charge = charge
        self.kpoint_coords = kpoint_coords
        self.perfect_kpoint_coords = perfect_kpoint_coords
        self.kpoint_weights = kpoint_weights
        self.eigenvalues = eigenvalues
        self.perfect_eigenvalues = perfect_eigenvalues
        self.perfect_symmops = perfect_symmops
        self.vbm = vbm
        self.cbm = cbm
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.fermi_level = fermi_level
        self.total_magnetization = total_magnetization
        self.orbital_character = orbital_character
        self.eigenvalue_correction = \
            dict(eigenvalue_correction) if eigenvalue_correction else None
        self.band_edges = band_edges

    @classmethod
    def from_files(cls,
                   unitcell: UnitcellCalcResults,
                   perfect: SupercellCalcResults,
                   defect: Defect):
        """ Parse defect eigenvalues.

        Args:
            unitcell (UnitcellCalcResults):
                UnitcellCalcResults object for band edge.
            perfect (SupercellCalcResults):
                SupercellDftResults object of perfect supercell for band edge in
                supercell.
            defect (Defect):
                Defect namedtuple object of a defect supercell DFT calculation
        """
        if unitcell.is_set_all is False:
            raise UnitcellCalcResultsError(
                "All the unitcell-related property is not set yet. ")

        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy.
        vbm, cbm = unitcell.band_edge
        supercell_cbm, supercell_vbm = perfect.eigenvalue_properties[1:3]

        return cls(name=defect.defect_entry.name,
                   charge=defect.defect_entry.charge,
                   kpoint_coords=defect.dft_results.kpoint_coords,
                   perfect_kpoint_coords=perfect.kpoint_coords,
                   kpoint_weights=defect.dft_results.kpoint_weights,
                   eigenvalues=defect.dft_results.eigenvalues,
                   perfect_eigenvalues=perfect.eigenvalues,
                   perfect_symmops=perfect.symmops,
                   vbm=vbm,
                   cbm=cbm,
                   supercell_vbm=supercell_vbm,
                   supercell_cbm=supercell_cbm,
                   fermi_level=defect.dft_results.fermi_level,
                   total_magnetization=defect.dft_results.total_magnetization,
                   orbital_character=defect.dft_results.orbital_character,
                   band_edges=defect.dft_results.band_edges)

    def plot(self, yrange=None, title=None):
        """ Plots the defect eigenvalues.
        Args:
            yrange (list):
                1x2 list for determining y energy range.
        """
        num_figure = len(self.eigenvalues.keys())
        fig, axs = plt.subplots(nrows=1, ncols=num_figure, sharey=True)
        fig.subplots_adjust(wspace=0)
        title = \
            "_".join([self.name, str(self.charge)]) if title is None else title
        fig.suptitle(title, fontsize=12)

        plt.title(title, fontsize=15)

        axs[0].set_ylabel("Eigenvalues (eV)", fontsize=12)

        if yrange is None:
            yrange = [self.supercell_vbm - 3, self.supercell_cbm + 3]

        mapping = self.kpt_mapping_to_perfect_kpt
        k_index = self.add_eigenvalues_to_plot(
            axs, self.eigenvalues, self.perfect_eigenvalues, mapping)

        x_labels = ["\n".join([str(i) for i in k]) for k in self.kpoint_coords]

        for i, s in enumerate(self.band_edges):
            axs[i].set_title(s.name.upper() + ": " + str(self.band_edges[s]))
            axs[i].set_ylim(yrange[0], yrange[1])
            axs[i].set_xlim(-1, k_index + 1)

            axs[i].get_xaxis().set_tick_params(direction='out')
            axs[i].xaxis.set_ticks_position('bottom')
            axs[i].set_xticks(np.arange(0, k_index + 1))

            axs[i].set_xticklabels(x_labels)

            axs[i].axhline(y=self.vbm, linewidth=0.3)
            axs[i].axhline(y=self.cbm, linewidth=0.3)
            axs[i].axhline(y=self.supercell_vbm, linewidth=0.3,
                           linestyle='dashed')
            axs[i].axhline(y=self.supercell_cbm, linewidth=0.3,
                           linestyle='dashed')
            axs[i].axhline(y=self.fermi_level, linewidth=1, linestyle="--")

        axs[num_figure - 1].annotate(
            "supercell\nvbm", xy=(k_index + 1, self.supercell_vbm - 0.2),
            fontsize=10)
        axs[num_figure - 1].annotate(
            "supercell\ncbm", xy=(k_index + 1, self.supercell_cbm - 0.2),
            fontsize=10)
        axs[num_figure - 1].annotate(
            "Fermi\nlevel", xy=(k_index + 1, self.fermi_level - 0.2),
            fontsize=10)
        axs[0].annotate("vbm", xy=(-1, self.vbm), fontsize=10)
        axs[0].annotate("cbm", xy=(-1, self.cbm), fontsize=10)


#         # def set_axis_style(ax, labels):
#         #     ax.get_xaxis().set_tick_params(direction='out')
#         #     ax.xaxis.set_ticks_position('bottom')
#         #     ax.set_xticks(np.arange(1, len(labels) + 1))
#         #     ax.set_xticklabels(labels)
#         #     ax.set_xlim(0.25, len(labels) + 0.75)

        plt.show()

    @property
    def kpt_mapping_to_perfect_kpt(self):
        mapping = []
        for kc in self.kpoint_coords:
            x = False
            for so in self.perfect_symmops:
                for i, pkc in enumerate(self.perfect_kpoint_coords):
                    if so.are_symmetrically_related(kc, pkc):
                        mapping.append(i)
                        x = True
                        break

                if x is True:
                    break

            else:
                raise ValueError("kpt in defects cannot be mapped.")

        return mapping

    @staticmethod
    def add_eigenvalues_to_plot(axs, eigenvalues, perfect_eigenvalues, mapping,
                                x_offset=0.3):
        occupied_eigenvalues = []
        occupied_x = []
        unoccupied_eigenvalues = []
        unoccupied_x = []
        partial_occupied_eigenvalues = []
        partial_occupied_x = []

        perfect_occupied_eigenvalues = []
        perfect_occupied_x = []
        perfect_unoccupied_eigenvalues = []
        perfect_unoccupied_x = []
        perfect_partial_occupied_eigenvalues = []
        perfect_partial_occupied_x = []

        for i, s in enumerate(eigenvalues.keys()):
            ax = axs[i]
            for k_index, eigen in enumerate(eigenvalues[s]):
                for band_index, e in enumerate(eigen):
                    occupancy = e[1]
                    if occupancy < 0.1:
                        unoccupied_eigenvalues.append(e[0])
                        unoccupied_x.append(k_index)

                    elif occupancy > 0.9:
                        occupied_eigenvalues.append(e[0])
                        occupied_x.append(k_index)
                    else:
                        partial_occupied_eigenvalues.append(e[0])
                        partial_occupied_x.append(k_index)

#                    if k_index == 1 and e[0] - eigen[band_index-1][0] > 0.1:
                    if band_index < len(eigen) - 1:
                        if eigen[band_index + 1][0] - e[0] > 0.1:
                            ax.annotate(str(band_index + 1),
                                        xy=(k_index + 0.05, e[0]),
                                        fontsize=10)

                for e in perfect_eigenvalues[Spin.up][mapping[k_index]]:
                    occupancy = e[1]
                    if occupancy < 0.1:
                        perfect_unoccupied_eigenvalues.append(e[0])
                        perfect_unoccupied_x.append(k_index + x_offset)

                    elif occupancy > 0.9:
                        perfect_occupied_eigenvalues.append(e[0])
                        perfect_occupied_x.append(k_index + x_offset)
                    else:
                        perfect_partial_occupied_eigenvalues.append(e[0])
                        perfect_partial_occupied_x.append(k_index + x_offset)

            ax.plot(occupied_x, occupied_eigenvalues, 'o')
            ax.plot(unoccupied_x, unoccupied_eigenvalues, 'o')
            ax.plot(partial_occupied_x, partial_occupied_eigenvalues, 'o')
            ax.plot(perfect_occupied_x, perfect_occupied_eigenvalues, 'o')
            ax.plot(perfect_unoccupied_x, perfect_unoccupied_eigenvalues, 'o')
            ax.plot(perfect_partial_occupied_x,
                    perfect_partial_occupied_eigenvalues, 'o')

        return k_index

    def to_json_file(self, filename):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def __repr__(self):
        pass
