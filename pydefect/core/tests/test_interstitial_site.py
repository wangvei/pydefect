# -*- coding: utf-8 -*-
import os
from collections import OrderedDict

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure

from pydefect.core.interstitial_site import InterstitialSiteSet, InterstitialSite

# __author__ = "Yu Kumagai"
# __maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")

class InterstitialSiteTest(PymatgenTest):
    def setUp(self):
        """ """

        site_name = "i1"
        representative_coords = [0.25, 0.25, 0.25]
        wyckoff = "b"
        site_symmetry = "m-3m"
        symmetry_multiplicity = 8
        coordination_distances = {"Mg": [2.12] * 6}
        method = "Voronoi"

        self.i = InterstitialSite(
                 site_name=site_name,
                 representative_coords=representative_coords,
                 wyckoff=wyckoff,
                 site_symmetry=site_symmetry,
                 symmetry_multiplicity=symmetry_multiplicity,
                 coordination_distances=coordination_distances,
                 method=method)

    def test_dict(self):
        d = self.i.as_dict()
        from_dict = InterstitialSite.from_dict(d)
        self.assertTrue(d == from_dict.as_dict())


class InterstitialSitesTest(PymatgenTest):

    def setUp(self):
        """ """
        structure = Structure.from_file("BPOSCAR-MgO")

        site_name_1 = "i1"
        representative_coords = [0.25, 0.25, 0.25]
        wyckoff = "b"
        site_symmetry = "m-3m"
        symmetry_multiplicity = 8
        coordination_distances = {"Mg": [2.12] * 6}
        method = "Voronoi"

        i1 = InterstitialSite(
            representative_coords=representative_coords,
            wyckoff=wyckoff,
            site_symmetry=site_symmetry,
            symmetry_multiplicity=symmetry_multiplicity,
            coordination_distances=coordination_distances,
            method=method)

        site_name_2 = "i2"
        representative_coords = [0.125, 0.125, 0.125]
        wyckoff = "b"
        site_symmetry = "m-3m"
        symmetry_multiplicity = 64
        coordination_distances = {"Mg": [2.12] * 6}
        method = "Voronoi"

        i2 = InterstitialSite(
            representative_coords=representative_coords,
            wyckoff=wyckoff,
            site_symmetry=site_symmetry,
            symmetry_multiplicity=symmetry_multiplicity,
            coordination_distances=coordination_distances,
            method=method)

        d = OrderedDict({site_name_1: i1, site_name_2: i2})

        self.interstitial_site_set = InterstitialSiteSet(structure, d)

    def test_dict(self):
        d = self.interstitial_site_set.as_dict()
        print(d)
#        from_dict = InterstitialSiteSet.from_dict(d)
#        self.assertTrue(d == from_dict.as_dict())

    def test_yaml(self):
        self.interstitial_site_set.site_set_to_yaml_file()

    def test_str(self):
        print(self.interstitial_site_set.interstitial_sites["i1"])

    def test_coords(self):
        actual = self.interstitial_site_set.coords
        expected = [[0.25, 0.25, 0.25], [0.125, 0.125, 0.125]]
        self.assertEqual(actual, expected)

    def test_site_names(self):
        actual = self.interstitial_site_set.site_names
        expected = ['i1', 'i2']
        self.assertEqual(actual, expected)

    def test_from_files(self):
        actual = InterstitialSiteSet.from_files(
            poscar="BPOSCAR-MgO", filename="interstitials.yaml").as_dict()
        expected = self.interstitial_site_set.as_dict()
        self.assertEqual(actual, expected)

    def test_add(self):
        coord = [0.175, 0.175, 0.175]
        self.interstitial_site_set.add_site(coord=coord)
        print(self.interstitial_site_set.interstitial_sites["i3"])
