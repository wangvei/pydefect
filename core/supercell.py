#!/usr/bin/env python

from abc import ABCMeta, abstractmethod
import json

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.dft_results import SupercellDftResults
from pydefect.input_maker.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


class Supercell(metaclass=ABCMeta):
    """
    Abstract class that is subclassed by a calculation done with a supercell.

    Args:
        dft_results (SupercellDftResults): SupercellDftResults class object
            |- final_structure (Structure)
            |- total_energy (float)
            |- eigenvalues (N_k x N_band numpy array)
             - electrostatic_potential (list)
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

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a Defect class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

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

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def set_vasp_results(self, directory_path, contcar_name="/CONTCAR",
                         outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        """
        Used to set some vasp results a posteriori.
        """
        self._dft_results = \
            SupercellDftResults.from_vasp_files(
                directory_path, contcar_name, outcar_name, vasprun_name)

    @abstractmethod
    def as_dict(self):
        pass

    @classmethod
    @abstractmethod
    def from_dict(cls, d):
        pass


class Perfect(Supercell):
    def as_dict(self):
        """
        Dictionary representation of Perfect.
        """
        d = {"dft_results": self._dft_results.as_dict()}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a  class object from a dictionary.
        """
        dft_results = SupercellDftResults.from_dict(d["dft_results"])

        return cls(dft_results)


class Defect(Supercell):
    """
    This class object holds some properties related to a defect.
    Args:
        dft_results (SupercellDftResults): SupercellDftResults class object
        defect_entry (DefectEntry): DefectEntry class object
            |- initial_structure
            |- removed_atom_index
            |- inserted_atom_index
            |- defect_coords
            |- in_name
            |- out_name
             - charge
    """
    def __init__(self, defect_entry, dft_results=None):
        self._defect_entry = defect_entry
        super().__init__(dft_results)

    def as_dict(self):
        """
        Dictionary representation of Defect.
        """
        d = {"defect_entry": self._defect_entry.as_dict(),
             "dft_results": self._dft_results.as_dict()}

        return d

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        defect_entry = DefectEntry.from_dict(d["defect_entry"])
        dft_results = SupercellDftResults.from_dict(d["dft_results"])

        return cls(defect_entry, dft_results)

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
