#!/usr/bin/env python

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.cell import PerfectSupercell, DefectSupercell
from pydefect.core.unitcell import Unitcell

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


class DefectSet:
    """
    """
    def __init__(self, unitcell=None, perfect=None, defects=[],
                 ewald_param=None):
        self._unitcell = unitcell
        self._perfect = perfect
        self._defects = defects
        self._ewald_param = ewald_param

    @classmethod
    def from_json_files(cls, unitcell_directory_path, perfect_directory_path,

#    @classmethod
#    def from_vasp_results(cls, unitcell_directory_path, perfect_directory_path,
#                          defect_directory_paths, contcar_name="CONTCAR",
#                          outcar_name="OUTCAR", vasprun_name="vasprun.xml"):
#        """
#        Constructs a class object from a set of directories.
#        """
#        unitcell = Unitcell.from_vasp_results(unitcell_directory_path,
#                                              contcar_name,
#                                              outcar_name,
#                                              vasprun_name)
#        perfect = PerfectSupercell.from_vasp_results(perfect_directory_path,
#                                            contcar_name,
#                                            outcar_name,
#                                            vasprun_name)
#        defects = []
#        for defect_directory_path in defect_directory_paths:
#            defects.append(DefectSupercell.from_vasp_results(defect_directory_path,
#                                                    contcar_name,
#                                                    outcar_name,
#                                                    vasprun_name))
#
#        cls(unitcell, perfect, defects)

#    def set_unitcell_from_vasp(self, unitcell_directory_path,
#                               contcar_name="CONTCAR",
#                               outcar_name="OUTCAR",
#                               vasprun_name="vasprun.xml"):
#        self._unitcell = Unitcell.from_vasp_results(unitcell_directory_path,
#                                                    contcar_name,
#                                                    outcar_name,
#                                                    vasprun_name)
#
#    def set_perfect_from_vasp(self):
#
    def set_defect_set_from_vasp(self):
        pass

    @classmethod
    def from_json(cls, unitcell_json, perfect_json, defect_json):
        return cls(Unitcell.json_load(unitcell_json),
                   PerfectSupercell.json_load(perfect_json),
                   DefectSupercell.json_load(defect_json))

    @property
    def ewald_param(self):
        return self._ewald_param

    @ewald_param.setter
    def ewald_param(self, ewald_param):
        self._ewald_param = ewald_param

    def data_check(self):
        pass

