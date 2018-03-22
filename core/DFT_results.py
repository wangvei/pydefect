#!/usr/bin/env python

from abc import ABCMeta
import json
from monty.json import MSONable

import numpy as np
from monty.json import MontyEncoder
from monty.serialization import loadfn
from copy import deepcopy

import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.electronic_structure.core import Spin

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


class SupercellDftResults:
    """
    DFT results for supercell systems both w/ and w/o a defect.
    """

    def __init__(self, final_structure, total_energy, eigenvalues,
                 electrostatic_potential, ewald_param=None):
        self._final_structure = final_structure
        self._total_energy = total_energy
        self._eigenvalues = eigenvalues
        self._electrostatic_potential = electrostatic_potential
        self._ewald_param = ewald_param

    @classmethod
    def from_vasp_files(cls, directory_path, contcar_name="/CONTCAR",
                        outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        """
        Args:
            directory_path (str): path of directory.
            contcar_name (str): Name of converged CONTCAR file.
                                Defaults to CONTCAR.
            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
            vasprun_name (str): Name of vasprun.xml file.
                                Defaults to vasprun.xml.
        """
        contcar = Poscar.from_file(directory_path + contcar_name)
        outcar = Outcar(directory_path + outcar_name)
        vasprun = Vasprun(directory_path + vasprun_name)

        final_structure = contcar.structure
        total_energy = outcar.final_energy
        eigenvalues = vasprun.eigenvalues
        electrostatic_potential = outcar.electrostatic_potential

        return cls(final_structure, total_energy, eigenvalues,
                   electrostatic_potential, ewald_param=None)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        return cls(d["final_structure"], d["total_energy"], eigenvalues,
                   d["electrostatic_potential"], d["ewald_param"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a Defect class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def final_structure(self):
        return self._final_structure

    @property
    def total_energy(self):
        return self._total_energy

    @property
    def eigenvalues(self):
        return self._eigenvalues

    @property
    def electrostatic_potential(self):
        return self._electrostatic_potential

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        Json-serializable dict representation.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self._eigenvalues.items()}

        d = {"final_structure":         self._final_structure,
             "total_energy":            self._total_energy,
             "eigenvalues":             eigenvalues,
             "electrostatic_potential": self._electrostatic_potential,
             "ewald_param":             self._ewald_param}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)


class UnitcellDftResults:
    """
    DFT result of unitcell systems w/o any defect.
    Args:
        final_structure (Structure): pmg Structure/IStructure class object.
        total_energy (float):
        static_dielectric_tensor (3x3 numpy array):
        ionic_dielectric_tensor (3x3 numpy array):
        eigenvalues (Nx3 numpy array):
    """

    def __init__(self, final_structure, total_energy, static_dielectric_tensor,
                 ionic_dielectric_tensor, eigenvalues):
        """

        """

        self._final_structure = final_structure
        self._total_energy = total_energy
        self._static_dielectric_tensor = static_dielectric_tensor
        self._ionic_dielectric_tensor = ionic_dielectric_tensor
        self._eigenvalues = eigenvalues

    @classmethod
    def from_vasp_files(cls, directory_path, contcar_name="/CONTCAR",
                        outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        """
        Args:
            directory_path (str): path of directory.
            contcar_name (str): Name of CONTCAR file. Defaults to CONTCAR.
            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
            vasprun_name (str): Name of vasprun.xml file.
                                Defaults to vasprun.xml.
        """
        contcar = Poscar.from_file(directory_path + contcar_name)
        outcar = Outcar(directory_path + outcar_name)
        vasprun = Vasprun(directory_path + vasprun_name)

        final_structure = contcar.structure
        total_energy = outcar.final_energy
#        static_dielectric_tensor = np.array(outcar.dielectric_tensor)
#        ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)
        static_dielectric_tensor = None
        ionic_dielectric_tensor = None
        eigenvalues = vasprun.eigenvalues

        return cls(final_structure, total_energy, static_dielectric_tensor,
                   ionic_dielectric_tensor, eigenvalues)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a  class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        return cls(d["final_structure"], d["total_energy"],
                   d["static_dielectric_tensor"], d["ionic_dielectric_tensor"],
                   d["eigenvalues"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a Defect class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    def set_dielectric_constants_from_outcar(
            self, directory_path, outcar_name="/OUTCAR"):
        outcar = Outcar(directory_path + outcar_name)
        self._static_dielectric_tensor = np.array(outcar.dielectric_tensor)
        self._ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)

    @property
    def structure(self):
        return self._structure

    @property
    def total_energy(self):
        return self._total_energy

    @property
    def static_dielectric_tensor(self):
        return self._static_dielectric_tensor

    @static_dielectric_tensor.setter
    def static_dielectric_tensor(self, static_dielectric_tensor):
        """ The matrix format is a 3 x 3 list. """
        self._static_dielectric_tensor = static_dielectric_tensor

    @property
    def ionic_dielectric_tensor(self):
        return self._ionic_dielectric_tensor

    @ionic_dielectric_tensor.setter
    def ionic_dielectric_tensor(self, ionic_dielectric_tensor):
        self._ionic_dielectric_tensor = ionic_dielectric_tensor

    @property
    def total_dielectric_tensor(self):
        return self._static_dielectric_tensor + self._ionic_dielectric_tensor

    @property
    def eigenvalues(self):
        return self._eigenvalues

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self._eigenvalues.items()}

        d = {"final_structure":           self._final_structure,
             "total_energy":              self._total_energy,
             "static_dielectric_tensor": self._static_dielectric_tensor,
             "ionic_dielectric_tensor":   self._ionic_dielectric_tensor,
             "eigenvalues":               eigenvalues}

        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)
