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
    Integrate a set of defect-related properties.
    Args:
        unitcell (UnitCell):
        perfect (PerfectSupercell):
        defects ([DefectSupercell, ..]):
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
        Construct a class object from a set of json files
        defect_json_files = [[*Va_O1_0*] ,[*Va_O1_1*], ...]
                                  |
        [defect_entry.json, dft_results.json, correction.json]

        """
        unitcell = Unitcell.json_load(unitcell_json_file)
        perfect = PerfectSupercell.json_load(perfect_dft_results_json_file)
        defects = [DefectSupercell(f[0], f[1]) for f in defect_json_files]
#        defects = [DefectSupercell(f[0], f[1], f[2]) for f in defect_json_files]

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
        Construct dft_results json files from the vasp calculation results
        The vasp file names are fixed to CONTCAR, OUTCAR, and vasprun.xml
        in the directory_paths.
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
            defect_dft_results.to_json_file(
                d_path + "/" + dft_results_json_name)