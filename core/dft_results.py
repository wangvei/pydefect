#!/usr/bin/env python

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn
from pymatgen import Spin
from pymatgen.io.vasp import Poscar, Outcar, Vasprun
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


class DftResults(metaclass=ABCMeta):
    def __init__(self, final_structure, total_energy, eigenvalues,
                 electrostatic_potential):
        self._final_structure = final_structure
        self._total_energy = total_energy
        self._eigenvalues = eigenvalues
        self._electrostatic_potential = electrostatic_potential

    @classmethod
    def from_vasp_files(cls, directory_path, contcar_name="/CONTCAR",
                        outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        # TODO: change "/POSCAR" to "POSCAR"
        # Then, bot format w/ and w/o "/" for directory_path must be allowed.
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
                   electrostatic_potential)

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
                   d["electrostatic_potential"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a Defect class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    def as_dict(self):
        """
        Dict representation of DefectInitialSetting class object.
        Json-serializable dict representation.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self._eigenvalues.items()}

        d = {"final_structure":         self._final_structure,
             "total_energy":            self._total_energy,
             "eigenvalues":             eigenvalues,
             "electrostatic_potential": self._electrostatic_potential}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @property
    def final_structure(self):
        return self._final_structure

    @property
    def eigenvalues(self):
        return self._eigenvalues

    @property
    def total_energy(self):
        return self._total_energy

    @property
    def electrostatic_potential(self):
        return self._electrostatic_potential


class SupercellDftResults(DftResults):
    """
    DFT results for supercell systems both w/ and w/o a defect.
    """
    pass


class UnitcellDftResults(DftResults):
    """
    DFT result of unitcell systems w/o any defect.
    Args:
        final_structure (Structure): pmg Structure/IStructure class object.
        total_energy (float):
        eigenvalues (Nx3 numpy array):
        static_dielectric_tensor (3x3 numpy array):
        ionic_dielectric_tensor (3x3 numpy array):
    """

    def __init__(self, final_structure, total_energy, eigenvalues,
                 electrostatic_potential, static_dielectric_tensor=None,
                 ionic_dielectric_tensor=None):
        """ """
        super().__init__(final_structure, total_energy, eigenvalues,
                         electrostatic_potential)

        self._static_dielectric_tensor = static_dielectric_tensor
        self._ionic_dielectric_tensor = ionic_dielectric_tensor

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a  class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        return cls(d["final_structure"], d["total_energy"], eigenvalues,
                   d["electrostatic_potential"], d["static_dielectric_tensor"],
                   d["ionic_dielectric_tensor"])

    def set_dielectric_constants_from_outcar(self, directory_path,
                                             outcar_name="/OUTCAR"):
        outcar = Outcar(directory_path + outcar_name)
        self._static_dielectric_tensor = np.array(outcar.dielectric_tensor)
        self._ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)

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

    def as_dict(self):
        """
        Dict representation of DefectInitialSetting class object.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self._eigenvalues.items()}

        d = {"final_structure":          self._final_structure,
             "total_energy":             self._total_energy,
             "eigenvalues":              eigenvalues,
             "electrostatic_potential": self._electrostatic_potential,
             "static_dielectric_tensor": self._static_dielectric_tensor,
             "ionic_dielectric_tensor":  self._ionic_dielectric_tensor}
        return d
