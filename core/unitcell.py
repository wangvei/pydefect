#!/usr/bin/env python
import json

from pydefect.core.dft_results import UnitcellDftResults

from monty.json import MontyEncoder
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "February 19, 2018"


class Unitcell:
    """
    This class object holds some properties related to a unit cell calculation.
    Args:
        dft_results (UnitcellDftResults): SupercellDftResults class object
    """

    def __init__(self, dft_results=None):
        self._dft_results = dft_results

    @classmethod
    def from_vasp_results(cls, directory_path, contcar_name="CONTCAR",
                          outcar_name="OUTCAR", vasprun_name="vasprun.xml"):
        dft_results = \
            UnitcellDftResults.from_vasp_files(
                directory_path, contcar_name, outcar_name, vasprun_name)

        return cls(dft_results)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        dft_results = UnitcellDftResults.from_dict(d["dft_results"])

        return cls(dft_results)

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a DefectSupercell class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def final_structure(self):
        if self._dft_results.final_structure:
            return self._dft_results.final_structure
        else:
            return None

    @property
    def total_energy(self):
        if self._dft_results.total_energy:
            return self._dft_results.total_energy
        else:
            return None

    @property
    def eigenvalues(self):
        if self._dft_results.eigenvalues:
            return self._dft_results.eigenvalues
        else:
            return None

    @property
    def static_dielectric_tensor(self):
        return self._dft_results.static_dielectric_tensor

    @property
    def ionic_dielectric_tensor(self):
        return self._dft_results.ionic_dielectric_tensor

    @property
    def total_dielectric_tensor(self):
        return self.static_dielectric_tensor + self.ionic_dielectric_tensor

    def set_vasp_results(self, directory_path, contcar_name="CONTCAR",
                         outcar_name="OUTCAR", vasprun_name="vasprun.xml"):
        self._dft_results = UnitcellDftResults.from_vasp_files(
            directory_path, contcar_name, outcar_name, vasprun_name)

    def set_dielectric_constants_from_outcar(self, directory_path,
                                             outcar_name="OUTCAR"):
        self._dft_results.set_dielectric_constants_from_outcar(directory_path,
                                                               outcar_name)

    def as_dict(self):
        """
        Dictionary representation of DefectSupercell class object.
        """
        d = {"dft_results":  self._dft_results.as_dict()}

        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)
