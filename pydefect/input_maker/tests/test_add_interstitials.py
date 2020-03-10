import os
import tempfile
from collections import OrderedDict
from copy import deepcopy
import filecmp
from pathlib import Path

from pydefect.core.interstitial_site import (
    InterstitialSiteSet, InterstitialSite)
from pydefect.input_maker.add_interstitials import add_interstitials
from pydefect.util.testing import PydefectTest

parent = Path(__file__).parent


class AddInterstitialsTest(PydefectTest):
    def test(self):
        add_interstitials(uc_coords=[0.25, 0.25, 0.25],
                          vicinage_radius=1.0,
                          dposcar=str(self.POSCARS_DIR / "POSCAR-MgO64atoms"),
                          interstitial_site_yaml="test_interstitial.yaml",
                          defect_in_file="defect_unittest.in")

    def tearDown(self) -> None:
        try:
            os.remove("test_interstitial.yaml")
        except FileNotFoundError:
            pass

