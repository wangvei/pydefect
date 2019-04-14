# -*- coding: utf-8 -*-

import json
import numpy as np

from monty.json import MSONable, MontyEncoder
from monty.serialization import loadfn

import matplotlib.pyplot as plt

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Kpoints

from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults


class DefectEigenvalue(MSONable):
    """
    A class related to eigenvalues in a defect calculation.
    """
    def __init__(self,
                 name: str,
                 kpoints: Kpoints,
                 eigenvalues: np.array,
                 vbm: float,
                 cbm: float,
                 band_gap: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 fermi_level: float,
                 total_magnetization: float,
                 eigenvalue_correction: dict = None,
                 deep_states: list = None):
        """
        Args:
            name (str):
                Name of a defect
            kpoints (Kpoints):
                Parsed IBZKPT.
            eigenvalues (N_spin x N_kpoint x N_band np.array):
                Numpy array for the electron eigenvalues in absolute scale.
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            band_gap (float):
                Band gap.
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
        """
        self.name = name
        self.kpoints = kpoints
        self.eigenvalues = eigenvalues
        self.vbm = vbm
        self.cbm = cbm
        self.band_gap = band_gap
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.fermi_level = fermi_level
        self.total_magnetization = total_magnetization
        if eigenvalue_correction:
            self.eigenvalue_correction = dict(eigenvalue_correction)
        else:
            self.eigenvalue_correction = None
        self.deep_states = list(deep_states) if deep_states else deep_states

    # def plot(self, energy_range=None):
    #     """
    #     Plots the defect eigenvalues.
    #     Args:
    #         energy_range (list):
    #             1x2 list for determining y energy range.
    #     """
    #     fig, ax = plt.subplots()
    #     plt.title(self.title, fontsize=15)

        # ax.set_xlabel("K points", fontsize=15)
        # ax.set_ylabel("Eigenvalues (eV)", fontsize=15)

        # plt.axhline(y=self.vbm, linewidth=0.3)
        # plt.axhline(y=self.cbm, linewidth=0.3)
        # plt.axhline(y=self.fermi_level, linewidth=2)

        # print(self.vbm, self.cbm, self.fermi_level)

        # if self.supercell_vbm < self.vbm - 0.05:
        #     plt.axhline(y=self.supercell_vbm, linewidth=0.3,
        #                 linestyle='dashed')
        #     # ax.annotate("supercell vbm",
        #     #             (supercell_vbm, y_min * 0.9 + y_max * 0.1),
        #     #             fontsize=10)
        # if self.supercell_cbm > self.cbm + 0.05:
        #     plt.axhline(y=self.supercell_cbm, linewidth=0.3,
        #                 linestyle='dashed')
        #     # ax.annotate("supercell cbm",
        #     #             (supercell_cbm, y_min * 0.9 + y_max * 0.1),
        #     #             fontsize=10)

        # data = [i[0] for i in self.eigenvalues[Spin.up][0]]
        # x_axis = [1 for i in range(len(data))]
        # ax.plot(x_axis, data, 'o')

        # def set_axis_style(ax, labels):
        #     ax.get_xaxis().set_tick_params(direction='out')
        #     ax.xaxis.set_ticks_position('bottom')
        #     ax.set_xticks(np.arange(1, len(labels) + 1))
        #     ax.set_xticklabels(labels)
        #     ax.set_xlim(0.25, len(labels) + 0.75)

        # set_axis_style(ax, "A")

        # plt.show()

    # @classmethod
    # def from_files(cls, unitcell, perfect, defect, system=""):
    #     """
    #     Args:
    #         unitcell (UnitcellCalcResults):
    #             UnitcellCalcResults object for band edge.
    #         perfect (SupercellCalcResults)
    #             SupercellDftResults object of perfect supercell for band edge in
    #             supercell.
    #         defect (Defect):
    #             Defect namedtuple object of a defect supercell DFT calculation
    #         system (str):
    #             System name used for the title.
    #     """
    #     # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy.
    #     vbm, cbm = unitcell.band_edge
    #     supercell_cbm, supercell_vbm = perfect.eigenvalue_properties[1:3]
    #     eigenvalues = defect.dft_results.eigenvalues
    #     fermi_level = defect.dft_results.fermi_level
    #     return cls(eigenvalues, vbm, cbm, supercell_vbm, supercell_cbm,
    #                fermi_level)

    # @classmethod
    # def from_dict(cls, d):
    #     """
    #     Constructs a class object from a dictionary.
    #     """
    #     # Programmatic access to enumeration members in Enum class.
    #     return cls(d["eigenvalues"], d["vbm"], d["cbm"], d["supercell_vbm"],
    #                d["supercell_cbm"], d["fermi_level"], d["title"])

    # @classmethod
    # def load_json(cls, filename):
    #     """
    #     Constructs a class object from a json file.
    #     """
    #     d = loadfn(filename)
    #     return cls.from_dict(d)

    # def as_dict(self):
    #     """
    #     Dict representation of the class object.
    #     """
    #     d = {"eigenvalues":   self.eigenvalues,
    #          "vbm":           self.vbm,
    #          "cbm":           self.cbm,
    #          "supercell_vbm": self.supercell_vbm,
    #          "supercell_cbm": self.supercell_cbm,
    #          "fermi_level":   self.fermi_level,
    #          "title":         self.title}
    #     return d

    def to_json_file(self, filename):
        """
        Returns a json file, named dft_results.json.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def is_shallow(self, magnetization_threshold=0.1, parse_wavecar=False):
        if abs(self.total_magnetization) > magnetization_threshold or \
                self.fermi_level > self.supercell_cbm or \
                self.fermi_level < self.supercell_vbm:
            return False

    def __str__(self):
        pass


