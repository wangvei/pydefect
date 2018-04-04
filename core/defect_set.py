#!/usr/bin/env python
# -*- coding: utf-8 -*-

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.cell import Unitcell, PerfectSupercell, DefectSupercell
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.dft_results import UnitcellDftResults, SupercellDftResults

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

    def __init__(self, unitcell=None, perfect=None, defects=[]):
        self._unitcell = unitcell
        self._perfect = perfect
        self._defects = defects

    @classmethod
    def from_json_files(cls,
                        unitcell_json_file,
                        perfect_dft_results_json_file,
                        defect_json_files):
        """
        Constructs a class object from a set of json files
        """
        unitcell = Unitcell.json_load(unitcell_json_file)
        perfect = PerfectSupercell.json_load(perfect_dft_results_json_file)
        defects = [DefectSupercell(f[0], f[1], f[2]) for f in defect_json_files]

        cls(unitcell, perfect, defects)

    @staticmethod
    def make_dft_results_json_files(unitcell_directory_path,
                                    dielectric_const_directory_path,
                                    perfect_directory_path,
                                    defect_directory_paths,
                                    unitcell_json_name="unitcell.json",
                                    perfect_json_name="perfect.json",
                                    dft_results_json_name="dft_results.json"
                                    ):
        """
        """
        unitcell_dft_results = \
            UnitcellDftResults.from_vasp_files(unitcell_directory_path)

        unitcell_dft_results.set_dielectric_constants_from_outcar(
            dielectric_const_directory_path)

        unitcell_dft_results.to_json_file(
            unitcell_directory_path + "/" + unitcell_json_name)

        perfect_dft_results = \
            SupercellDftResults.from_vasp_files(perfect_directory_path)

        perfect_dft_results.to_json_file(
            perfect_directory_path + "/" + perfect_json_name)

        for d_path in defect_directory_paths:
            defect_dft_results = SupercellDftResults.from_vasp_files(d_path)
            defect_dft_results.to_json_file(d_path + "/" + dft_results_json_name)

    # @classmethod
    # def from_directory_paths(cls, unitcell_directory_path,
    #                          perfect_directory_path,
    #                          defect_directory_paths):

        # unitcell_json_file = unitcell_directory_path + "unitcell.json"
        # perfect_dft_results_json_file = \
        #     perfect_directory_path + "perfect.json"
        # defect_json_files = [[d + "defect_entry.json",
        #                       d + "dft_results.json",
        #                       d + "correction.json"]
        #                      for d in defect_directory_paths]

        # cls.from_json_files(unitcell_json_file, perfect_dft_results_json_file,
        #                     defect_json_files)

    # def data_check(self):
    #     pass
