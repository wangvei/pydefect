# -*- coding: utf-8 -*-

import json
import numpy as np

from collections import defaultdict
from monty.json import MSONable, MontyEncoder

import matplotlib.pyplot as plt

from pymatgen.electronic_structure.core import Spin

from pydefect.analysis.defect_energies import Defect
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.error_classes import UnitcellCalcResultsError
from pydefect.vasp_util.util import calc_orbital_similarity
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
                 participation_ratio: dict = None,
                 orbital_character: dict = None,
                 perfect_orbital_character: dict = None,
                 eigenvalue_correction: dict = None,
                 deep_states: dict = None,
                 shallow: dict = None):
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
            deep_states (list):
                Band indices corresponding to the deep defect states.
                ex. {Spin.up: ["vbm", "cbm"]}
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
        self.participation_ratio = participation_ratio
        self.orbital_character = orbital_character
        self.perfect_orbital_character = perfect_orbital_character
        self.eigenvalue_correction = \
            dict(eigenvalue_correction) if eigenvalue_correction else None
        self.deep_states = deep_states if deep_states else defaultdict(list)
        self.shallow = shallow if shallow else defaultdict(dict)

    @classmethod
    def from_files(cls,
                   unitcell: UnitcellCalcResults,
                   perfect: SupercellCalcResults,
                   defect: Defect,
                   diagnose_band_edges: bool = True):
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
                   participation_ratio=defect.dft_results.participation_ratio,
                   orbital_character=defect.dft_results.orbital_character,
                   perfect_orbital_character=perfect.orbital_character)

    def diagnose_shallow_states(self,
                                localized_criterion: float = 0.25,
                                similarity_criterion: float = 0.3):
        """

        Args:
             localized_criterion (float):
                 Determines whether the eigenstate is localized.
             similarity_criterion:
                 Determines whether the eigenstate is a host state.
        """
        if all([self.participation_ratio, self.orbital_character,
                self.perfect_orbital_character]) is False:
            logger.warning("Diagnosing shallow states is impossible.")
            return False

        for spin in self.participation_ratio:
            for band_edge in "vbm", "cbm":
                participation_ratio = self.participation_ratio[spin][band_edge]

                # Consider the situation where perfect has only Spin.up.
                try:
                    perfect = self.perfect_orbital_character[spin][band_edge]
                except KeyError:
                    perfect = self.perfect_orbital_character[Spin.up][band_edge]

                dissimilarity = \
                    calc_orbital_similarity(
                        self.orbital_character[spin][band_edge], perfect)

                print("participation_ratio, dissimilarity")
                print(participation_ratio, dissimilarity)
                print("spin, band_edge")
                print(spin, band_edge)

                if (participation_ratio < localized_criterion) \
                        and (dissimilarity < similarity_criterion):
                    # band edge is a host state.
                    pass
                elif (participation_ratio > localized_criterion) \
                        and (dissimilarity > similarity_criterion):
                    # band edge is a deep state
                    self.deep_states[spin].append(band_edge)
                elif (participation_ratio < localized_criterion) \
                        and (dissimilarity > similarity_criterion):
                    # band edge is a perturbed host sate
                    self.shallow[spin][band_edge] = True
                else:
                    logger.warning("Orbitals are localized but similar to the "
                                   "band edge of the perfect supercell.")

    def plot(self, yrange=None, add_perfect=True, title=None):
        """
        Plots the defect eigenvalues.
        Args:
            yrange (list):
                1x2 list for determining y energy range.
        """
        fig, ax = plt.subplots()
        title = \
            "_".join([self.name, str(self.charge)]) if title is None else title

        plt.title(title, fontsize=15)

        ax.set_xlabel("K points", fontsize=15)
        ax.set_ylabel("Eigenvalues (eV)", fontsize=15)

        if yrange is None:
            yrange = [self.supercell_vbm - 2, self.supercell_cbm + 2]

        plt.axhline(y=self.vbm, linewidth=0.3)
        plt.axhline(y=self.cbm, linewidth=0.3)
        plt.axhline(y=self.fermi_level, linewidth=1, linestyle="--")

        mapping = self.kpt_mapping_to_perfect_kpt
        k_index = self.add_eigenvalues_to_plot(ax, self.eigenvalues,
                                               self.perfect_eigenvalues, mapping)
        ax.set_ylim(yrange[0], yrange[1])

        ax.get_xaxis().set_tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(np.arange(1, k_index + 1))
        ax.set_xticklabels(self.kpoint_coords * 2)
        #        ax.set_xticklabels([]

        ax.annotate("vbm", xy=(k_index + 1, self.vbm), fontsize=10)
        ax.annotate("cbm", xy=(k_index + 1, self.cbm), fontsize=10)

        if self.supercell_vbm < self.vbm - 0.05:
            plt.axhline(y=self.supercell_vbm, linewidth=0.3,
                        linestyle='dashed')
            ax.annotate("supercell vbm",
                        xy=(k_index + 1, self.supercell_vbm),
                        fontsize=10)
        if self.supercell_cbm > self.cbm + 0.05:
            plt.axhline(y=self.supercell_cbm, linewidth=0.3,
                        linestyle='dashed')
            ax.annotate("supercell cbm",
                        xy=(k_index + 1, self.supercell_cbm),
                        fontsize=10)

        # def set_axis_style(ax, labels):
        #     ax.get_xaxis().set_tick_params(direction='out')
        #     ax.xaxis.set_ticks_position('bottom')
        #     ax.set_xticks(np.arange(1, len(labels) + 1))
        #     ax.set_xticklabels(labels)
        #     ax.set_xlim(0.25, len(labels) + 0.75)

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
    def add_eigenvalues_to_plot(ax, eigenvalues, perfect_eigenvalues, mapping,
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

        k_index = 0

        for s in Spin.up, Spin.down:
            for k, eigen in enumerate(eigenvalues[s]):
                k_index += 1
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
                            ax.annotate(str(band_index),
                                        xy=(k_index + 0.05, e[0]),
                                        fontsize=10)

                for e in perfect_eigenvalues[Spin.up][mapping[k]]:
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

            plt.axvline(x=k_index + 0.5, linewidth=1, linestyle="--")

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

    def is_shallow(self, magnetization_threshold=0.1, parse_wavecar=False):
        if abs(self.total_magnetization) > magnetization_threshold or \
                self.fermi_level > self.supercell_cbm or \
                self.fermi_level < self.supercell_vbm:
            return False

    def __str__(self):
        pass
