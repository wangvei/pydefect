# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.vasp_util.script.calc_band_center import calc_o2p_band_center

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "core", "MgO", "unitcell",
                        "dielectric_constants")


class Test(unittest.TestCase):

    def test(self):
        vasprun = os.path.join(test_dir, "vasprun.xml")
        calc_o2p_band_center(vasprun_file=vasprun)

