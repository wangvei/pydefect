#!/usr/bin/env python

from abc import ABCMeta
import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.DFT_results import SupercellDftResults
from pydefect.input_maker.defect_entry import DefectEntry

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

    @property
    def ewald_param(self):
        return self._ewald_param

    @ewald_param.setter
    def ewald_param(self, ewald_param):
        self._ewald_param = ewald_param


class Perfect(Supercell):
    pass



    def test_ewald_param(self):
        ewald_param = 0.1
        self._MgO_perfect.ewald_param = ewald_param
        self.assertEqual(self._MgO_perfect.ewald_param, ewald_param)
