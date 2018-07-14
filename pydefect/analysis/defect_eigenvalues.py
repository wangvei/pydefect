# -*- coding: utf-8 -*-

import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import Spin

from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults


class DefectEigenvalue:
    """
    A class related to eigenvalues in a defect calculation.
    """
    def __init__(self, eigenvalues, vbm, cbm, supercell_vbm, supercell_cbm,
                 fermi_level, title=None):
        """
        Args:
            eigenvalues (N_spin x N_kpoint x N_band np.array):
                Numpy array for the electron eigenvalues in absolute scale.
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            supercell_vbm (float):
                Valence band maximum in the perfect supercell.
            supercell_cbm (float):
                Conduction band minimum in the perfect supercell.
            title (str):
                Title of the system.
        """
        #TODO: Add correction
        self._eigenvalues = eigenvalues
        self._vbm = vbm
        self._cbm = cbm
        self._supercell_vbm = supercell_vbm
        self._supercell_cbm = supercell_cbm
        self._fermi_level = fermi_level
        self._title = title

    def plot(self, energy_range=None):
        """
        Plots the defect eigenvalues.
        Args:
            energy_range (list):
                1x2 list for determining y energy range.
        """
        fig, ax = plt.subplots()
        plt.title(self._title, fontsize=15)

        ax.set_xlabel("K points", fontsize=15)
        ax.set_ylabel("Eigenvalues (eV)", fontsize=15)

        plt.axhline(y=self._vbm, linewidth=0.3)
        plt.axhline(y=self._cbm, linewidth=0.3)
        plt.axhline(y=self._fermi_level, linewidth=2)

        print(self._vbm, self._cbm, self._fermi_level)

        if self._supercell_vbm < self._vbm - 0.05:
            plt.axhline(y=self._supercell_vbm, linewidth=0.3,
                        linestyle='dashed')
            # ax.annotate("supercell vbm",
            #             (supercell_vbm, y_min * 0.9 + y_max * 0.1),
            #             fontsize=10)
        if self._supercell_cbm > self._cbm + 0.05:
            plt.axhline(y=self._supercell_cbm, linewidth=0.3,
                        linestyle='dashed')
            # ax.annotate("supercell cbm",
            #             (supercell_cbm, y_min * 0.9 + y_max * 0.1),
            #             fontsize=10)

        data = [i[0] for i in self._eigenvalues[Spin.up][0]]
        x_axis = [1 for i in range(len(data))]
        ax.plot(x_axis, data, 'o')

        def set_axis_style(ax, labels):
            ax.get_xaxis().set_tick_params(direction='out')
            ax.xaxis.set_ticks_position('bottom')
            ax.set_xticks(np.arange(1, len(labels) + 1))
            ax.set_xticklabels(labels)
            ax.set_xlim(0.25, len(labels) + 0.75)

        set_axis_style(ax, "A")

        plt.show()

    @classmethod
    def from_files(cls, unitcell, perfect, defect, system=""):
        """
        Args:
            unitcell (UnitcellDftResults):
                UnitcellDftResults object for band edge.
            perfect (SupercellDftResults)
                SupercellDftResults object of perfect supercell for band edge in
                supercell.
            defect (Defect):
                Defect namedtuple object of a defect supercell DFT calculation
            system (str):
                System name used for the title.
        """
        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy.
        vbm, cbm = unitcell.band_edge
        supercell_cbm, supercell_vbm = perfect.eigenvalue_properties[1:3]
        eigenvalues = defect.dft_results.eigenvalues
        fermi_level = defect.dft_results.fermi_level
        return cls(eigenvalues, vbm, cbm, supercell_vbm, supercell_cbm,
                   fermi_level, title=system)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        # Programmatic access to enumeration members in Enum class.
        return cls(d["eigenvalues"], d["vbm"], d["cbm"], d["supercell_vbm"],
                   d["supercell_cbm"], d["fermi_level"], d["title"])

    @classmethod
    def load_json(cls, filename):
        """
        Constructs a class object from a json file.
        """
        d = loadfn(filename)
        return cls.from_dict(d)

    def as_dict(self):
        """
        Dict representation of the class object.
        """
        d = {"eigenvalues":   self._eigenvalues,
             "vbm":           self._vbm,
             "cbm":           self._cbm,
             "supercell_vbm": self._supercell_vbm,
             "supercell_cbm": self._supercell_cbm,
             "fermi_level":   self._fermi_level,
             "title":         self._title}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file, named dft_results.json.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def __str__(self):
        pass

    @property
    def eigenvalues(self):
        return self._eigenvalues

    @property
    def vbm(self):
        return self._vbm

    @property
    def cbm(self):
        return self._cbm

    @property
    def supercell_vbm(self):
        return self._supercell_vbm

    @property
    def supercell_cbm(self):
        return self._supercell_cbm

    @property
    def title(self):
        return self._title

