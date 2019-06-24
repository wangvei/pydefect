# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.core.irreducible_site import IrreducibleSite

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class PriorInfoTest(unittest.TestCase):

    def setUp(self):
        """ """
        irreducible_name = "O1"
        element = "O"
        first_index = 33
        last_index = 64
        representative_coords = [0, 0, 0]
        wyckoff = "m3m"
        site_symmetry = "Oh"
        coordination_distance = {"Mg": [2, 2, 2, 2, 2, 2]}

        self.O1_MgO = IrreducibleSite(irreducible_name,
                                      element,
                                      first_index,
                                      last_index,
                                      representative_coords,
                                      wyckoff,
                                      site_symmetry,
                                      coordination_distance)

    def test_dict(self):
        """ round trip test of to_dict and from_dict """
        d = self.O1_MgO.as_dict()
        prior_info_from_dict = IrreducibleSite.from_dict(d)
        self.assertTrue(d == prior_info_from_dict.as_dict())


