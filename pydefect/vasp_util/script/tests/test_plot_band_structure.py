# -*- coding: utf-8 -*-

import unittest

from pydefect.vasp_util.script.plot_band_structure import plot_band_structure

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class PlotBandStructureTest(unittest.TestCase):

    def test(self):
        plot_band_structure("KPOINTS", "vasprun-finish.xml")


if __name__ == '__main__':
    unittest.main()
