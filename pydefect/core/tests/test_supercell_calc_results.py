# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from collections import defaultdict

import numpy as np
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results \
    import SupercellCalcResults, defaultdict_to_dict
from pymatgen.electronic_structure.core import Spin
from pymatgen.util.testing import PymatgenTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefaultdictToDictTest(PymatgenTest):
    def setUp(self):
        """ """
        self.d = defaultdict(lambda: defaultdict(dict))
        self.d[1][2][4] = 10

    def test(self):
        print(self.d)
        print(defaultdict_to_dict(self.d))


class SupercellDftResultsTest(PymatgenTest):

    def setUp(self):
        """ """
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


if __name__ == "__main__":
    unittest.main()

