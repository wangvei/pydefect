# -*- coding: utf-8 -*-
import os

from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

from pydefect.core.interstitial_site import InterstitialSiteSet, InterstitialSite
from pydefect.input_maker.interstitial_site_maker import InterstitialSiteMaker

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


# __author__ = "Yu Kumagai"
# __maintainer__ = "Yu Kumagai"


class InterstitialSitesTest(PymatgenTest):

    def setUp(self):
        """ """
        self.interstitial_maker = \
            InterstitialSiteMaker(interstitial_set=InterstitialSiteSet())

    def test_add(self):
        filename = os.path.join(test_dir, "POSCAR-MgO64atoms")
        structure = Structure.from_file(filename)
        specie = "H"
        coords = [0.125, 0.125, 0.125]

        self.interstitial_maker.add_interstitial(structure, specie, coords,
                                                 "oct")

        print(self.interstitial_maker.interstitial_set[0])
