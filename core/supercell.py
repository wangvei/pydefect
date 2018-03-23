#!/usr/bin/env python

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.DFT_results import SupercellDftResults
from pydefect.input_maker.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def get_nions(defect_structure):
    """
    Return numbers of ions for elements in defect_structure.
    Example: Al1Mg63O6
        return: [1, 63, 64]

    """
    nions = [int(i)
             for i in defect_structure.to(fmt="poscar").split("\n")[6].split()]
    return nions


class Supercell(metaclass=ABCMeta):
    """
    """
    def __init__(self, dft_results=None):
        self._dft_results = dft_results

    @classmethod
    def from_vasp_results(cls, directory_path, contcar_name="/CONTCAR",
                          outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        dft_results = \
            SupercellDftResults.from_vasp_files(
                directory_path, contcar_name, outcar_name, vasprun_name)
        return cls(dft_results)

    @property
    def final_structure(self):
        return self._dft_results.final_structure

    @property
    def total_energy(self):
        return self._dft_results.total_energy

    @property
    def eigenvalues(self):
        return self._dft_results.eigenvalues

    @property
    def electrostatic_potential(self):
        return self._dft_results.electrostatic_potential

    @property
    def ewald_param(self):
        return self._ewald_param

    @ewald_param.setter
    def ewald_param(self, ewald_param):
        self._ewald_param = ewald_param

    def set_vasp_results(self, directory_path, contcar_name="/CONTCAR",
                         outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        """
        Used to set some vasp results a posteriori.
        """
        self._dft_results = \
            SupercellDftResults.from_vasp_files(
                directory_path, contcar_name, outcar_name, vasprun_name)


class Perfect(Supercell):
    pass


class Defect(Supercell):
    """
    This class object holds some properties related to a defect.
    Args:

    """
    def __init__(self, dft_results=None, defect_entry=None):
        super().__init__(dft_results)
        self._defect_entry = defect_entry

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        # Expansion of _dft_results is necessary,
        _dft_results = []
        for i in d["irreducible_sites"]:
            irreducible_sites.append(IrreducibleSite.from_dict(i))

        return cls(d["initial_structure"], d["removed_atom_index"],
                   d["inserted_atom_index"], d["defect_coords"], d["in_name"],
                   d["out_name"], d["charge"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a Defect class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def initial_structure(self):
        return self._defect_entry.initial_structure

    @property
    def removed_atom_index(self):
        return self._defect_entry.removed_atom_index

    @property
    def inserted_atom_index(self):
        return self._defect_entry.inserted_atom_index

    @property
    def defect_coords(self):
        return self._defect_entry.defect_coords

    @property
    def in_name(self):
        return self._defect_entry.in_name

    @property
    def out_name(self):
        return self._defect_entry.out_name

    @property
    def charge(self):
        return self._defect_entry.charge
