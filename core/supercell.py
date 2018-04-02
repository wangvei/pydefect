#!/usr/bin/env python

from abc import ABCMeta, abstractmethod
import json
import operator
import warnings

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
            |- eigenvalues (N_k-point x N_band numpy array)
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
        if self._dft_results:
            return self._dft_results.final_structure
        else:
            return None

    @property
    def total_energy(self):
        if self._dft_results:
            return self._dft_results.total_energy
        else:
            return None

    @property
    def eigenvalues(self):
        if self._dft_results:
            return self._dft_results.eigenvalues
        else:
            return None

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
            |- removed_atoms
            |- inserted_atoms
            |- changes_of_num_elements
            |- in_name
            |- out_name
    """

    def __init__(self, defect_entry, dft_results=None):
        self._defect_entry = defect_entry
        super().__init__(dft_results)
        self._relative_total_energy = None
        self._relative_potential = None
        self._relative_atomic_displacements = None

    def as_dict(self):
        """
        Dictionary representation of Defect class object.
        """
        d = {"defect_entry": self._defect_entry.as_dict(),
             "dft_results":  self._dft_results.as_dict()}

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
    def removed_atoms(self):
        return self._defect_entry.removed_atoms

    @property
    def inserted_atoms(self):
        return self._defect_entry.inserted_atoms

    @property
    def changes_of_num_elements(self):
        return self._defect_entry.changes_of_num_elements

    @property
    def charge(self):
        return self._defect_entry.charge

    @property
    def in_name(self):
        return self._defect_entry.in_name

    @property
    def out_name(self):
        return self._defect_entry.out_name

    @property
    def relative_total_energy(self):
        if self._relative_total_energy:
            return self._relative_total_energy
        else:
            warnings.warn("relative_total_energy is not set.")
            return None

    @property
    def relative_potential(self):
        if self._relative_potential:
            return self._relative_potential
        else:
            warnings.warn("relative_potential is not set.")
            return None

#    @property
#    def relative_atomic_displacements(self):
#        if self._relative_atomic_displacements:
#            return self._relative_atomic_displacements
#        else:
#            warnings.warn("relative_atomic_displacements is not set.")
#            return None

    def set_relative_values(self, perfect):
        """
        Set relative values with respect to those of perfect calc.
        Args:
            perfect (Perfect):
        """
        self._relative_total_energy = self.total_energy - perfect.total_energy
        self._relative_potential = list(map(operator.sub,
                                            self.electrostatic_potential,
                                            perfect.electrostatic_potential))
        # TODO: write below
#        self._relative_atomic_displacements = \
