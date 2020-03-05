# -*- coding: utf-8 -*-

from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.util.testing import PydefectTest


class IrreducibleSiteTest(PydefectTest):

    def setUp(self):
        irreducible_name = "O1"
        element = "O"
        first_index = 32
        last_index = 63
        representative_coords = [0, 0, 0]
        wyckoff = "m3m"
        site_symmetry = "Oh"
        cutoff = 2.37
        coordination_distance = {"Mg": [2, 2, 2, 2, 2, 2]}

        self.mgo_o1 = IrreducibleSite(irreducible_name,
                                      element,
                                      first_index,
                                      last_index,
                                      representative_coords,
                                      wyckoff,
                                      site_symmetry,
                                      cutoff,
                                      coordination_distance)

    def test_dict(self):
        """ round trip test of to_dict and from_dict """
        d = self.mgo_o1.as_dict()
        irreducible_site_dict = IrreducibleSite.from_dict(d)
        self.assertTrue(d == irreducible_site_dict.as_dict())

    def test_msonable(self):
        self.assertMSONable(self.mgo_o1)
