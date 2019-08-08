# -*- coding: utf-8 -*-
import os
from collections import OrderedDict
import tempfile
import filecmp

from pymatgen.core.structure import Structure

from pydefect.util.testing import PydefectTest
from pydefect.core.interstitial_site import InterstitialSiteSet, \
    InterstitialSite

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InterstitialSiteTest(PydefectTest):
    def setUp(self):
        representative_coords = [0.25, 0.25, 0.25]
        wyckoff = "b"
        site_symmetry = "m-3m"
        symmetry_multiplicity = 8
        coordination_distances = {"Mg": [2.12] * 6}
        method = "manual"

        self.mgo_interstitial = \
            InterstitialSite(
                representative_coords=representative_coords,
                wyckoff=wyckoff,
                site_symmetry=site_symmetry,
                symmetry_multiplicity=symmetry_multiplicity,
                coordination_distances=coordination_distances,
                method=method)

    def test_repr(self):
        expected = """representative_coords: [0.25, 0.25, 0.25]
wyckoff: b
site_symmetry: m-3m
symmetry_multiplicity: 8
coordination_distances: {'Mg': [2.12, 2.12, 2.12, 2.12, 2.12, 2.12]}
method: manual"""
        self.assertEqual(expected, str(self.mgo_interstitial))

    def test_dict(self):
        d = self.mgo_interstitial.as_dict()
        rounded_d = InterstitialSite.from_dict(d).as_dict()
        self.assertEqual(d, rounded_d)


class InterstitialSiteSetTest(PydefectTest):

    def setUp(self):
        self.structure = self.get_structure_by_name("MgO64atoms")

        i1 = InterstitialSite(
            representative_coords=[0.125, 0.125, 0.125],
            wyckoff="b",
            site_symmetry="m-3m",
            symmetry_multiplicity=64,
            coordination_distances={"Mg": [2.12] * 4, "O": [2.12] * 4},
            method="manual")

        i2 = InterstitialSite(
            representative_coords=[0, 0.125, 0.125],
            wyckoff="b",
            site_symmetry="m-3m",
            symmetry_multiplicity=64,
            coordination_distances={"Mg": [2.12] * 6},
            method="Voronoi")

        d = OrderedDict({"i1": i1, "i2": i2})

        self.interstitial_site_set = InterstitialSiteSet(self.structure, d)

    def test_dict(self):
        d = self.interstitial_site_set.as_dict()
        repr_coords = d["interstitial_site_set"]["i1"]["representative_coords"]
        self.assertEqual(repr_coords, [0.125, 0.125, 0.125])
        rounded_d = InterstitialSiteSet.from_dict(d).as_dict()
        self.assertEqual(d, rounded_d)

    def test_yaml(self):
        tmp_file = tempfile.NamedTemporaryFile().name
        self.interstitial_site_set.site_set_to_yaml_file(tmp_file)
        self.assertTrue(filecmp.cmp("expected_interstitials.yaml", tmp_file))

    def test_str(self):
        actual = str(self.interstitial_site_set.interstitial_sites["i1"])
        expected = """representative_coords: [0.125, 0.125, 0.125]
wyckoff: b
site_symmetry: m-3m
symmetry_multiplicity: 64
coordination_distances: {'Mg': [2.12, 2.12, 2.12, 2.12], 'O': [2.12, 2.12, 2.12, 2.12]}
method: manual"""
        self.assertEqual(expected, actual)

    def test_coords(self):
        actual = self.interstitial_site_set.coords
        expected = [[0.125, 0.125, 0.125], [0, 0.125, 0.125]]
        self.assertEqual(expected, actual)

    def test_from_files(self):
        actual = InterstitialSiteSet.from_files(
            structure=self.structure, filename="expected_interstitials.yaml").as_dict()
        expected = self.interstitial_site_set.as_dict()
        self.assertEqual(expected, actual)

    def test_add(self):
        coord = [0.175, 0.175, 0.175]
        self.interstitial_site_set.add_site(coord=coord)
        print(self.interstitial_site_set.interstitial_sites["i3"])

    def test_add_from_charge_density(self):
        chgcar_name = self.get_filename("core/CHGCAR-MgO8atoms")
        self.interstitial_site_set.add_sites_from_charge_density(chgcar_name)

