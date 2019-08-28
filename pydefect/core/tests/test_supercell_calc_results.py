# -*- coding: utf-8 -*-
from copy import deepcopy
import tempfile

import numpy as np
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import (
    analyze_procar, SupercellCalcResults)
from pymatgen.io.vasp.outputs import Outcar, Vasprun, Procar
from pymatgen.electronic_structure.core import Spin
from pydefect.util.testing import PydefectTest
from pydefect.util.tools import all_combination

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class AnalyzeProcarTest(PydefectTest):
    def setUp(self) -> None:
        """ Va_O in the 2+ charge state in 64-atom supercells"""
        hob_index = {Spin.up: 124, Spin.down: 124}
        procar = self.get_object_by_name(
            Procar, ["defects", "MgO", "Va_O1_2", "PROCAR"])
        vasprun = self.get_object_by_name(
            Vasprun, ["defects", "MgO", "Va_O1_2", "vasprun.xml"])
        eigenvalues = vasprun.eigenvalues
        structure = self.get_structure_by_name("MgO64atoms-Va_O1_2")
        neighboring_sites = [0, 4, 16, 17, 24, 26]

        (self.band_edge_energies, self.orbital_character,
         self.participation_ratio) \
            = analyze_procar(hob_index=hob_index,
                             procar=procar,
                             eigenvalues=eigenvalues,
                             structure=structure,
                             neighboring_sites=neighboring_sites)

    def test_band_edge_energies(self):
        expected = {Spin.up: {'hob': {'top': 5.5148, 'bottom': 5.5148},
                              'lub': {'top': 8.6662, 'bottom': 8.6662}},
                    Spin.down: {'hob': {'top': 5.5148, 'bottom': 5.5148},
                                'lub': {'top': 8.6662, 'bottom': 8.6662}}}
        self.assertEqual(expected, self.band_edge_energies)

    def test_orbital_character(self):
        expected = {Spin.up: {'hob': {'top': {'Mg': {'s': 0.018, 'p': 0.036,
                                                     'd': 0.018, 'f': 0.0},
                                              'O': {'s': 0.018, 'p': 0.216,
                                                    'd': 0.0, 'f': 0.0}},
                                      'bottom': {'Mg': {'s': 0.018, 'p': 0.036,
                                                        'd': 0.018, 'f': 0.0},
                                                 'O': {'s': 0.018, 'p': 0.216,
                                                       'd': 0.0, 'f': 0.0}}},
                              'lub': {'top': {'Mg': {'s': 0.174, 'p': 0.006,
                                                     'd': 0.0, 'f': 0.0},
                                              'O': {'s': 0.199, 'p': 0.114,
                                                    'd': 0.0, 'f': 0.0}},
                                      'bottom': {'Mg': {'s': 0.174, 'p': 0.006,
                                                        'd': 0.0, 'f': 0.0},
                                                 'O': {'s': 0.199, 'p': 0.114,
                                                       'd': 0.0, 'f': 0.0}}}}}
        expected[Spin.down] = deepcopy(expected[Spin.up])
        for k1, k2, k3, k4, k5, v in all_combination(expected):
            self.assertAlmostEqual(
                v, self.orbital_character[k1][k2][k3][k4][k5], 3)

    def test_participation_ratio(self):
        expected = {Spin.up: {'hob': 0.235294, 'lub': 0.060852},
                    Spin.down: {'hob': 0.235294, 'lub': 0.060852}}
        for k1, k2, v in all_combination(expected):
            self.assertAlmostEqual(v, self.participation_ratio[k1][k2], 5)


class SupercellDftResultsTest(PydefectTest):

    def setUp(self):
        """ Va_O in the 2+ charge state in 64-atom supercells"""
        # self._defect_entry_perfect = \
        #     DefectEntry.load_json(os.path.join(
        #         test_dir, "MgO/defects/perfect", "defect_entry.json"))

        self._MgO_perfect = \
            SupercellCalcResults.from_vasp_files(
                directory_path=os.path.join(test_dir, "MgO/defects/perfect"),
                procar=True)
#                defect_entry=self._defect_entry_perfect)

        defect_entry_MgO_Va_O1_2 = \
            DefectEntry.load_json(os.path.join(
                test_dir, "MgO/defects/Va_O1_2", "defect_entry.json"))

        self._MgO_Va_O1_2 = \
            SupercellCalcResults.from_vasp_files(
                directory_path=os.path.join(test_dir, "MgO/defects/Va_O1_2"),
                procar=True,
                defect_entry=defect_entry_MgO_Va_O1_2)

    def test_nothing(self):
        pass
#        print(self._MgO_Va_O1_2.displacements)
        print(self._MgO_Va_O1_2.participation_ratio)
        print(self._MgO_Va_O1_2.orbital_character)
        print(self._MgO_Va_O1_2.band_edge_energies)
#        print(self._MgO_perfect.band_edge_energies)

    def test_from_vasp_files(self):
        # CAUTION: When constructing Structure object from Structure.from_file
        #          velocities are not stored.
        #          Therefore, equality check of Structure objects returns False.
        #          If the structure is converted via poscar file format, it may
        #          be solved.

        # energy
        expected = -93.64519527
        self.assertEqual(self._MgO_Va_O1_2.total_energy, expected)

        # total_magnetization
        expected = 3.4e-06
        self.assertEqual(self._MgO_Va_O1_2.total_magnetization, expected)

        # eigenvalue: test only a single point
        expected = [-1.43572e+01, 1.0]
        self.assertArrayEqual(self._MgO_Va_O1_2.eigenvalues[Spin.up][0][0],
                              expected)

        # electrostatic_potential
        expected = [-34.5752, -35.5301, -35.5344, -35.5357, -35.5346, -35.5326,
                    -35.5265, -34.5754, -70.027, -70.026, -70.0253, -70.0274,
                    -70.0271, -70.0272, -70.4853]
        self.assertEqual(self._MgO_Va_O1_2.electrostatic_potential, expected)

    def test_dict(self):
        """Round trip test"""
        d = self._MgO_Va_O1_2.as_dict()
        dd = SupercellCalcResults.from_dict(d).as_dict()
        self.assertEqual(d, dd)

    def test_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self._MgO_Va_O1_2.to_json_file(tmp_file.name)
        actual = SupercellCalcResults.load_json(tmp_file.name)
        np.testing.assert_equal(actual.eigenvalues[Spin.up],
                                self._MgO_Va_O1_2.eigenvalues[Spin.up])

