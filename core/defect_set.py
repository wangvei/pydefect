#!/usr/bin/env python

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.cell import Unitcell, PerfectSupercell, DefectSupercell

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
    def __init__(self,
                 unitcell_directory_path,
                 defect_directory_paths,
                 perfect_directory_path="perfect",
                 contcar_name="CONTCAR",
                 outcar_name="OUTCAR",
                 vasprun_name="vasprun.xml",
                 does_update_json=True):

        self.unitcell = unitcell
        self.perfect = perfect
        self.defects = defects
        self.ewald_param = ewald_param

#    @classmethod
#    def from_json_files(cls, unitcell_directory_path, perfect_directory_path):


    @classmethod
    def from_vasp_results(cls,
                          unitcell_directory_path,
                          defect_directory_paths,
                          perfect_directory_path="perfect",
                          contcar_name="CONTCAR",
                          outcar_name="OUTCAR",
                          vasprun_name="vasprun.xml",
                          does_update_json=True):
        """
        Constructs a class object from a set of directories.
        """
        unitcell = Unitcell.from_vasp_results(unitcell_directory_path,
                                              contcar_name,
                                              outcar_name,
                                              vasprun_name)
        perfect = PerfectSupercell.from_vasp_results(perfect_directory_path,
                                                     contcar_name,
                                                     outcar_name,
                                                     vasprun_name)

        defects = []
        for path in defect_directory_paths:
            d = DefectSupercell.from_vasp_results(path,
                                                  contcar_name,
                                                  outcar_name,
                                                  vasprun_name)

        defects = [DefectSupercell.from_vasp_results(path,
                                                     contcar_name,
                                                     outcar_name,
                                                     vasprun_name)
                   for path in defect_directory_paths]

        if does_update_json:
            unitcell.to_json_file(unitcell_directory_path + "/unitcell.json")
            perfect.to_json_file(unitcell_directory_path + "/perfect.json")
            for defect, path in zip(defects, defect_directory_paths):
                defect.to_json_file(path + "/defect.json")

        cls(unitcell, perfect, defects)

    def set_directory_paths(self):

    # def set_unitcell_from_vasp(self,
    #                            unitcell_directory_path,
    #                            contcar_name="CONTCAR",
    #                            outcar_name="OUTCAR",
    #                            vasprun_name="vasprun.xml"):
    #     self.unitcell = Unitcell.from_vasp_results(unitcell_directory_path,
    #                                                contcar_name,
    #                                                outcar_name,
    #                                                vasprun_name)

    # def set_perfect_from_vasp(self,
    #                           perfect_directory_path,
    #                           contcar_name="CONTCAR",
    #                           outcar_name="OUTCAR",
    #                           vasprun_name="vasprun.xml"):
    #     self.perfect = PerfectSupercell.\
    #         from_vasp_results(perfect_directory_path,
    #                           contcar_name,
    #                           outcar_name,
    #                           vasprun_name)

    # @classmethod
    # def from_json(cls, unitcell_json, perfect_json, defect_json):
    #     return cls(Unitcell.json_load(unitcell_json),
    #                PerfectSupercell.json_load(perfect_json),
    #                DefectSupercell.json_load(defect_json))

    # def to_json(self):


    @property
    def ewald_param(self):
        return self.ewald_param

    @ewald_param.setter
    def ewald_param(self, ewald_param):
        self.ewald_param = ewald_param

    def data_check(self):
        pass

