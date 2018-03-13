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


# def eigenvalues_to_list(eigenvalues):


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
        path_contcar = directory_path + contcar_name
        path_outcar = directory_path + outcar_name
        path_vasprun = directory_path + vasprun_name

        contcar = Poscar.from_file(path_contcar)
        outcar = Outcar(path_outcar)
        vasprun = Vasprun(path_vasprun)

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
        print(d["eigenvalues"])

        for spin, v in d["eigenvalues"].items():
            if spin == "1":
                eigenvalues[Spin.up] = np.array(v)
            elif spin == "2":
                eigenvalues[Spin.down] = np.array(v)

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

    @property
    def ewald_param(self):
        return self._ewald_param

    @ewald_param.setter
    def ewald_param(self, ewald_param):
        self._ewald_param = ewald_param

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        Json-serializable dict representation.
        To make a dictionary for eigenvalus, we use the same way with pymatgen.
        See. e.g.,
        pymatgen/io/vasp/outputs.html#BSVasprun
        if self.projected_eigenvalues:
            vout['projected_eigenvalues'] = {
                str(spin): v.tolist()
                for spin, v in self.projected_eigenvalues.items()}
        """
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
            for k, v in self.as_dict().items():
                json.dump({k: v}, fw, indent=2, cls=MontyEncoder)


# class UnitcellDFTResults:
#    """
#    DFT result of unitcell without any defect.
#    """
#    def __init__(self, structure, total_energy, static_dielectric_tensor,
#                 ionic_dielectric_tensor, eigenvalues):
#
#        self._structure = structure
#        self._total_energy = total_energy
#        self._static_dielectric_tensor = static_dielectric_tensor
#        self._ionic_dielectric_tensor = ionic_dielectric_tensor
#        self._eigenvalues = eigenvalues
#
#    @classmethod
#    def from_vasp_files(cls, directory_path, contcar_name="/CONTCAR",
#                        outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
#        """
#        Args:
#            directory_path (str): path of directory.
#            contcar_name (str): Name of CONTCAR file. Defaults to CONTCAR.
#            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
#            vasprun_name (str): Name of vasprun.xml file.
#                                Defaults to vasprun.xml.
#        """
#        path_contcar = directory_path + contcar_name
#        path_outcar = directory_path + outcar_name
#        path_vasprun = directory_path + vasprun_name
#        contcar = Poscar.from_file(path_contcar)
#        outcar = Outcar(path_outcar)
#        vasprun = Vasprun(path_vasprun)
#        structure = contcar.structure
#        total_energy = outcar.final_energy
#        static_dielectric_tensor = np.array(outcar.dielectric_tensor)
#        ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)
#        eigenvalues = vasprun.eigenvalues
#
#        return cls(structure, total_energy, static_dielectric_tensor,
#                   ionic_dielectric_tensor, eigenvalues)
#
#    @classmethod
#    def from_dict(cls, d):
#        """
#        Constructs a  class object from a dictionary.
#        """
#        return cls(d)
#
#    @property
#    def structure(self):
#        return self._structure
#
#    @property
#    def total_energy(self):
#        return self._total_energy
#
#    @property
#    def static_dielectric_tensor(self):
#        return self._static_dielectric_tensor
#
#    @property
#    def ionic_dielectric_tensor(self):
#        return self._ionic_dielectric_tensor
#
#    @property
#    def total_dielectric_tensor(self):
#        return self._static_dielectric_tensor + self._ionic_dielectric_tensor
#
#    @property
#    def eigenvalue(self):
#        return self._eigenvalues

#    def as_dict(self):
#        """
#        Dict representation of DefectSetting class object.
#        """
#        d = {"structure":                 self._structure,
#             "total_energy":              self._total_energy,
#             "static_dielectric_tensor ": self._static_dielectric_tensor,
#             "ionic_dielectric_tensor":   self._ionic_dielectric_tensor,
#             "eigenvalues":                self._eigenvalues}
#
#        return d
