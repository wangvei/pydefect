# -*- coding: utf-8 -*-

import json
from collections import defaultdict
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MSONable, MontyEncoder
from pydefect.analysis.defect import Defect
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.logger import get_logger
from pydefect.core.defect_name import DefectName

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class DefectEigenvalue(MSONable):
    """ A class related to eigenvalues in a defect calculation. """

    def __init__(self,
                 name: str,
                 charge: int,
                 annotation: str,
                 kpoint_coords: list,
                 eigenvalues: np.array,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 fermi_level: float,
                 magnetization: float,
                 vbm: float = None,
                 cbm: float = None,
                 orbital_character: dict = None,
                 eigenvalue_correction: dict = None,
                 band_edge_states: dict = None):
        """
        Args:
            name (str):
                Name of a defect
            charge (int):
                Defect charge.
            annotation (str):
                Annotation
            kpoint_coords (list):
                List of k-point coordinates
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
            magnetization (float):
                Total magnetization.
            eigenvalue_correction (dict):
                Dict with key of correction name and value of correction value.
            band_edge_states (dict):
                Band edge states at each spin channel.
                None: no in gap state.
                "Donor PHS": Donor-type perturbed host state (PHS).
                "Acceptor PHS": Acceptor-type PHS.
                "Localized state": With in-gap localized state.
                    ex. {Spin.up: None, Spin.down:"Localized state}
        """
        self.name = name
        self.charge = charge
        self.annotation = annotation
        self.kpoint_coords = kpoint_coords
        self.eigenvalues = eigenvalues
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.fermi_level = fermi_level
        self.magnetization = magnetization
        self.vbm = vbm
        self.cbm = cbm
        self.orbital_character = orbital_character
        self.eigenvalue_correction = eigenvalue_correction
        self.band_edge_states = band_edge_states

    @classmethod
    def from_files(cls,
                   unitcell: UnitcellCalcResults,
                   defect: Defect):
        """ Parse defect eigenvalues.

        Args:
            unitcell (UnitcellCalcResults):
                UnitcellCalcResults object for band edge.
            defect (Defect):
                Defect namedtuple object of a defect supercell DFT calculation
        """
        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy scale.
        return cls(name=defect.name,
                   charge=defect.charge,
                   annotation=defect.annotation,
                   kpoint_coords=defect.kpoint_coords,
                   eigenvalues=defect.eigenvalues,
                   supercell_vbm=defect.supercell_vbm,
                   supercell_cbm=defect.supercell_cbm,
                   fermi_level=defect.fermi_level,
                   magnetization=defect.magnetization,
                   vbm=unitcell.band_edge[0],
                   cbm=unitcell.band_edge[1],
                   orbital_character=defect.orbital_character,
                   band_edge_states=defect.band_edge_states)

    def plot(self,
             y_range: Optional[str] = None,
             title: Optional[str] = None,
             filename: Optional[str] = None) -> None:
        """Plots the defect eigenvalues.

        Args:
            y_range (list):
                1x2 list for determining y energy range.
            title (str):
                Title of the plot
            filename (str):
                Filename when the plot is saved; otherwise show plot.

        Returns:
            Matplotlib pyplot object.
        """
        num_figure = len(self.eigenvalues.keys())
        fig, axs = plt.subplots(nrows=1, ncols=num_figure, sharey='all')
        fig.subplots_adjust(wspace=0)

        title = title or DefectName(self.name, self.charge, self.annotation)
        fig.suptitle(title, fontsize=12)
        plt.title(title, fontsize=15)

        axs[0].set_ylabel("Eigenvalues (eV)", fontsize=12)

        if y_range is None:
            y_range = [self.supercell_vbm - 3, self.supercell_cbm + 3]

        last_k_index = self._add_eigenvalues_to_plot(axs)

        x_labels = ["\n".join([str(i) for i in k]) for k in self.kpoint_coords]

        for i, s in enumerate(self.band_edge_states):
            # show band-edge states
            axs[i].set_title(f"{s.name.upper()}: {self.band_edge_states[s]}")
            axs[i].set_ylim(y_range[0], y_range[1])
            axs[i].set_xlim(-1, last_k_index + 1)

            axs[i].get_xaxis().set_tick_params(direction='out')
            axs[i].xaxis.set_ticks_position('bottom')
            axs[i].set_xticks(np.arange(0, last_k_index + 1))

            axs[i].set_xticklabels(x_labels)

            axs[i].axhline(y=self.supercell_vbm, linewidth=0.7, linestyle="-",
                           color='r')
            axs[i].axhline(y=self.supercell_cbm, linewidth=0.7, linestyle="-",
                           color='r')
            axs[i].axhline(y=self.fermi_level, linewidth=1, linestyle="--",
                           color='g')
            axs[i].axhline(y=self.vbm, linewidth=0.7, linestyle="-", color='b')
            axs[i].axhline(y=self.cbm, linewidth=0.7, linestyle="-", color='b')

        xy = [last_k_index + 1, self.fermi_level - 0.2]
        axs[num_figure - 1].annotate("Fermi\nlevel", fontsize=8, color='g',
                                     xy=xy)
        axs[0].annotate("vbm", xy=(-1, self.vbm), fontsize=8, color='b')
        axs[0].annotate("cbm", xy=(-1, self.cbm), fontsize=8, color='b')

        plt.savefig(filename) if filename else plt.show()

    def _add_eigenvalues_to_plot(self, axs):

        k_index = 0
        for i, s in enumerate(self.eigenvalues.keys()):
            eigenvalues = defaultdict(list)
            k_indices = defaultdict(list)

            ax = axs[i]
            for k_index, eigen in enumerate(self.eigenvalues[s]):
                for band_index, (energy, occupation) in enumerate(eigen):
                    if occupation < 0.1:
                        eigenvalues["unoccupied"].append(energy)
                        k_indices["unoccupied"].append(k_index)

                    elif occupation > 0.9:
                        eigenvalues["occupied"].append(energy)
                        k_indices["occupied"].append(k_index)
                    else:
                        eigenvalues["partial"].append(energy)
                        k_indices["partial"].append(k_index)

                    if band_index < len(eigen) - 1:
                        if (energy < self.fermi_level and
                            eigen[band_index + 1][0] - energy > 0.2) or \
                                (energy > self.fermi_level and
                                 energy - eigen[band_index - 1][0] > 0.2):
                            ax.annotate(str(band_index + 1),
                                        xy=(k_index + 0.05, energy),
                                        va='center',
                                        fontsize=10)

            ax.plot(k_indices["occupied"], eigenvalues["occupied"], 'o')
            ax.plot(k_indices["unoccupied"], eigenvalues["unoccupied"], 'o')
            ax.plot(k_indices["partial"], eigenvalues["partial"], 'o')

        return k_index

    def to_json_file(self, filename):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

