# -*- coding: utf-8 -*-
# import os

from pymatgen.util.testing import PymatgenTest
from collections import OrderedDict

from pydefect.core.interstitial_site import InterstitialSites, InterstitialSite

# __author__ = "Yu Kumagai"
# __maintainer__ = "Yu Kumagai"


class InterstitialSitesTest(PymatgenTest):

    def setUp(self):
        """ """
        site_name = "i1"
        representative_coords = [0.25, 0.25, 0.25]
        wyckoff = "b"
        site_symmetry = "m-3m"
        symmetry_multiplicity = 8
        coordination_distances = {"Mg": [2.12] * 6}
        method = "Voronoi"

        i1 = InterstitialSite(
            site_name=site_name,
            representative_coords=representative_coords,
            wyckoff=wyckoff,
            site_symmetry=site_symmetry,
            symmetry_multiplicity=symmetry_multiplicity,
            coordination_distances=coordination_distances,
            method=method)

        site_name = "i2"
        representative_coords = [0.125, 0.125, 0.125]
        wyckoff = "b"
        site_symmetry = "m-3m"
        symmetry_multiplicity = 64
        coordination_distances = {"Mg": [2.12] * 6}
        method = "Voronoi"

        i2 = InterstitialSite(
            site_name=site_name,
            representative_coords=representative_coords,
            wyckoff=wyckoff,
            site_symmetry=site_symmetry,
            symmetry_multiplicity=symmetry_multiplicity,
            coordination_distances=coordination_distances,
            method=method)

        self.interstitial_sites = InterstitialSites([i1, i2])

    def test_dict(self):
        d = self.interstitial_sites.as_dict()
        from_dict = InterstitialSites.from_dict(d)
        print(from_dict.as_dict())
        self.assertTrue(d == from_dict.as_dict())

    def test_yaml(self):
        self.interstitial_sites.to_yaml_file()



    def test_from_interstitial_in(self):
        """
        """

        actual = InterstitialSites.from_interstitial_in()
        print(actual.as_dict())
        print(self.interstitial_sites.as_dict())
        self.assertTrue(actual.as_dict() == self.interstitial_sites.as_dict())

