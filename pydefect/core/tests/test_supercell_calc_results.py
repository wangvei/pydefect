# -*- coding: utf-8 -*-
from copy import deepcopy
import tempfile

import numpy as np
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import (
     ProcarDefectProperty, SupercellCalcResults)
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.electronic_structure.core import Spin
from pydefect.util.testing import PydefectTest
from pydefect.util.tools import flatten_dict


class ProcarDefectPropertyTest(PydefectTest):
    def setUp(self) -> None:
        """ Va_O in the 2+ charge state in 64-atom supercells"""
        # TODO: Fix the hob_index toss123 and change related values.
        #       The true hob_index is 123 but is fine for unittest.
        hob_index = {Spin.up: 124, Spin.down: 124}
        procar = self.get_object_by_name(
            Procar, ["defects", "MgO", "Va_O1_2", "PROCAR"])
        vasprun = self.get_object_by_name(
            Vasprun, ["defects", "MgO", "Va_O1_2", "vasprun.xml"])
        eigenvalues = vasprun.eigenvalues
        structure = self.get_structure_by_name("MgO64atoms-Va_O1_2")
        neighboring_sites = [0, 4, 16, 17, 24, 26]

        self.prop = ProcarDefectProperty.analyze_procar(
            hob_index=hob_index,
            procar=procar,
            eigenvalues=eigenvalues,
            structure=structure,
            neighboring_sites=neighboring_sites)

    def test_band_edge_energies(self):
        expected = {  Spin.up: {'hob': {'top': 5.5148, 'bottom': 5.5148},
                                'lub': {'top': 8.6662, 'bottom': 8.6662}},
                    Spin.down: {'hob': {'top': 5.5148, 'bottom': 5.5148},
                                'lub': {'top': 8.6662, 'bottom': 8.6662}}}
        self.assertEqual(expected, self.prop.band_edge_energies)

    def test_orbital_character(self):
        expected = \
            {Spin.up: {'hob':   {'top': {'Mg': {'s': 0.018, 'p': 0.036,
                                                'd': 0.018, 'f': 0.0},
                                          'O': {'s': 0.018, 'p': 0.216,
                                                'd': 0.0,   'f': 0.0}},
                              'bottom': {'Mg': {'s': 0.018, 'p': 0.036,
                                                'd': 0.018, 'f': 0.0},
                                          'O': {'s': 0.018, 'p': 0.216,
                                                'd': 0.0,   'f': 0.0}}},
                       'lub':   {'top': {'Mg': {'s': 0.174, 'p': 0.006,
                                                'd': 0.0,   'f': 0.0},
                                          'O': {'s': 0.199, 'p': 0.114,
                                                'd': 0.0,   'f': 0.0}},
                              'bottom': {'Mg': {'s': 0.174, 'p': 0.006,
                                                'd': 0.0,   'f': 0.0},
                                          'O': {'s': 0.199, 'p': 0.114,
                                                'd': 0.0,   'f': 0.0}}}}}
        expected[Spin.down] = deepcopy(expected[Spin.up])
        for k1, k2, k3, k4, k5, v in flatten_dict(expected):
            self.assertAlmostEqual(
                v, self.prop.orbital_character[k1][k2][k3][k4][k5], 3)

    def test_participation_ratio(self):
        expected = {  Spin.up: {'hob': 0.235294, 'lub': 0.060852},
                    Spin.down: {'hob': 0.235294, 'lub': 0.060852}}
        for k1, k2, v in flatten_dict(expected):
            self.assertAlmostEqual(v, self.prop.participation_ratio[k1][k2], 5)


class SupercellDftResultsTest(PydefectTest):

    def setUp(self):
        """ Va_O in the 2+ charge state in 64-atom supercells"""
        self.mgO_perfect = \
            SupercellCalcResults.from_vasp_files(
                directory_path=self.DEFECTS_MGO_DIR / "perfect")

        filepath = ["defects", "MgO", "Va_O1_2", "defect_entry.json"]
        defect_entry = self.get_object_by_name(DefectEntry.load_json, filepath)

        self.mgo_va_o1_2 = SupercellCalcResults.from_vasp_files(
            directory_path=self.DEFECTS_MGO_DIR / "Va_O1_2",
            defect_entry=defect_entry)

    def test_from_vasp_files(self):
        # CAUTION: When constructing Structure object from Structure.from_file
        #          velocities are not stored, so equality check of Structure
        #          objects returns False. If the structure is converted via
        #          poscar file format, it may be solved.

        # energy
        expected = -399.85095628
        actual = self.mgo_va_o1_2.total_energy
        self.assertAlmostEqual(expected, actual, 5)

        # total_magnetization
        expected = 0.0
        actual = self.mgo_va_o1_2.total_magnetization
        self.assertAlmostEqual(expected, actual, 5)

        # eigenvalue: test only a single point
        expected = [-1.40215e+01, 1.0]
        actual = self.mgo_va_o1_2.eigenvalues[Spin.up][0][0]
        self.assertArrayAlmostEqual(expected, actual, 5)

    def test_dict(self):
        expected = self.mgo_va_o1_2.as_dict()
        actual = SupercellCalcResults.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self.mgo_va_o1_2.to_json_file(tmp_file.name)
        actual = SupercellCalcResults.load_json(tmp_file.name)
        np.testing.assert_equal(actual.eigenvalues[Spin.up],
                                self.mgo_va_o1_2.eigenvalues[Spin.up])

    def test_msonable(self):
        self.assertMSONable(self.mgO_perfect)
        self.assertMSONable(self.mgo_va_o1_2)
