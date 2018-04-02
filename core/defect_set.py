#!/usr/bin/env python

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.supercell import Perfect, Defect
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

    def set_unitcell_from_json(self):
        pass

    def set_unitcell_from_vasp(self):
        pass

    def set_perfect_from_json(self):
        pass

    def set_perfect_from_vasp(self):
        pass

    def set_defect_set_from_vasp(self):
        pass

    def set_defect_set_from_json(self):
        pass

    @classmethod
    def from_json(cls, unitcell_json, perfect_json, defect_json):
        return cls(Unitcell.json_load(unitcell_json),
                   Perfect.json_load(perfect_json),
                   Defect.json_load(defect_json))

    @property
    def ewald_param(self):
        return self._ewald_param

    @ewald_param.setter
    def ewald_param(self, ewald_param):
        self._ewald_param = ewald_param

    def data_check(self):
        pass

