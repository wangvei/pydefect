# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.vasp_util.script.plot_dos import get_dos_plot

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


class PlotBandStructureTest(unittest.TestCase):

    def test(self):
        vasprun = os.path.join(test_dir, "vasprun.xml")
        a = get_dos_plot(vasprun, site=True)
        a.show()
