# -*- coding: utf-8 -*-
import filecmp
import tempfile
from collections import OrderedDict
from copy import deepcopy
import filecmp

from pydefect.core.interstitial_site import (
    InterstitialSiteSet, InterstitialSite)
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InterstitialSiteTest(PydefectTest):
    def setUp(self):
        representative_coords = [0.25, 0.25, 0.25]
        wyckoff = "b"
        site_symmetry = "m-3m"
        multiplicity = 8
        coordination_distances = {"Mg": [2.12] * 6}
        method = "manual"

        self.mgo_interstitial = \
            InterstitialSite(
                representative_coords=representative_coords,
                wyckoff=wyckoff,
                site_symmetry=site_symmetry,
                multiplicity=multiplicity,
                coordination_distances=coordination_distances,
                method=method)

    def test_repr(self):
        expected = """representative_coords: [0.25, 0.25, 0.25]
wyckoff: b
site_symmetry: m-3m
multiplicity: 8
coordination_distances: {'Mg': [2.12, 2.12, 2.12, 2.12, 2.12, 2.12]}
method: manual"""
        self.assertEqual(expected, str(self.mgo_interstitial))

    def test_dict(self):
        d = self.mgo_interstitial.as_dict()
        rounded_d = InterstitialSite.from_dict(d).as_dict()
        self.assertEqual(d, rounded_d)

    def test_msonable(self):
        self.assertMSONable(self.mgo_interstitial)


class InterstitialSiteSetTest(PydefectTest):

    def setUp(self):
        """
        i1:
          representative_coords: [0.125, 0.125, 0.125]
          wyckoff: c
          site_symmetry: -43m
          multiplicity: 64
          coordination_distances:
            Mg: [1.84, 1.84, 1.84, 1.84]
            O: [1.84, 1.84, 1.84, 1.84]
          method: manual
        i2:
          representative_coords: [0.0, 0.125, 0.125]
          wyckoff: d
          site_symmetry: m.mm
          multiplicity: 192
          coordination_distances:
            Mg: [1.5, 1.5, 2.6, 2.6, 2.6, 2.6]
            O: [1.5, 1.5, 2.6, 2.6, 2.6, 2.6]
          method: manual
        """
        self.structure = self.get_structure_by_name("MgO64atoms")

        i1 = InterstitialSite(
            representative_coords=[0.125, 0.125, 0.125],
            wyckoff="c",
            site_symmetry="-43m",
            multiplicity=64,
            coordination_distances={"Mg": [1.84] * 4, "O": [1.84] * 4},
            method="manual")

        i2 = InterstitialSite(
            representative_coords=[0, 0.125, 0.125],
            wyckoff="d",
            site_symmetry="m.mm",
            multiplicity=192,
            coordination_distances={"Mg": [1.5] * 2 + [2.6] * 4,
                                    "O": [1.5] * 2 + [2.6] * 4},
            method="manual")

        d = OrderedDict({"i1": i1, "i2": i2})

        self.interstitial_site_set = InterstitialSiteSet(self.structure, d)

    def test_dict(self):
        d = self.interstitial_site_set.as_dict()
        rounded_d = InterstitialSiteSet.from_dict(d).as_dict()
        self.assertEqual(d, rounded_d)

    def test_yaml(self):
        tmp_file = tempfile.NamedTemporaryFile().name
        self.interstitial_site_set.site_set_to_yaml_file(tmp_file)
        self.assertTrue(filecmp.cmp("expected_interstitials.yaml", tmp_file))

    def test_str(self):
        actual = str(self.interstitial_site_set.interstitial_sites["i1"])
        expected = """representative_coords: [0.125, 0.125, 0.125]
wyckoff: c
site_symmetry: -43m
multiplicity: 64
coordination_distances: {'Mg': [1.84, 1.84, 1.84, 1.84], 'O': [1.84, 1.84, 1.84, 1.84]}
method: manual"""
        self.assertEqual(expected, actual)

    def test_coords(self):
        actual = self.interstitial_site_set.coords
        expected = [[0.125, 0.125, 0.125], [0, 0.125, 0.125]]
        self.assertEqual(expected, actual)

    def test_from_files(self):
        actual = InterstitialSiteSet.from_files(
            dposcar=self.structure,
            yaml_filename="expected_interstitials.yaml").as_dict()
        expected = self.interstitial_site_set.as_dict()
        self.assertEqual(expected, actual)

    def test_add(self):
        coords = [[0.125, 0, 0]]
        interstitial_set = deepcopy(self.interstitial_site_set)
        interstitial_set.add_sites(frac_coords=coords)
        expected = InterstitialSite(
            representative_coords=[0.125, 0, 0],
            wyckoff="c",
            site_symmetry="-43m",
            multiplicity=64,
            coordination_distances={"Mg": [1.05] + [2.35] * 4,
                                    "O": [1.05] + [2.35] * 4},
            method="manual")
        self.assertEqual(expected.as_dict(),
                         interstitial_set.interstitial_sites["i3"].as_dict())

    def test_msonable(self):
        self.assertMSONable(self.interstitial_site_set)