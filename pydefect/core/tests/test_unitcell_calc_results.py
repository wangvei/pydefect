# -*- coding: utf-8 -*-

from copy import deepcopy
import numpy as np
import os
import tempfile
import unittest

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.util.testing import PymatgenTest

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results \
    import SupercellCalcResults
from pydefect.util.structure_tools import defect_center, distances_from_point
from pydefect.core.unitcell_calc_results import UnitcellCalcResults

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class UnitcellDftResultsTest(PymatgenTest):

    def setUp(self):
        """ """
        self.unitcell = UnitcellCalcResults(band_edge=None,
                                            static_dielectric_tensor=None,
                                            ionic_dielectric_tensor=None,
                                            total_dos=None,
                                            volume=None)

        MgO_band_edges = [2.9978, 7.479]
        MgO_static_dielectric_tensor = [[3.166727, 0, 0],
                                        [0, 3.166727, 0],
                                        [0, 0, 3.166727]]
        MgO_ionic_dielectric_tensor = [[9.102401, 0, 0],
                                       [0, 9.102448, 0],
                                       [0, 0, 9.102542]]
        MgO_fictitious_dos = [[0] * 301] * 2
        MgO_fictitious_dos[1][300] = 23.3688

        MgO_volume = 19.1659131591
        self.MgO_unitcell = UnitcellCalcResults(
            band_edge=MgO_band_edges,
            static_dielectric_tensor=MgO_static_dielectric_tensor,
            ionic_dielectric_tensor=MgO_ionic_dielectric_tensor,
            total_dos=MgO_fictitious_dos,
            volume=MgO_volume)

    def test_set_static_dielectric_tensor(self):
        self.unitcell.static_dielectric_tensor = 3.166727
        self.assertArrayAlmostEqual(self.unitcell.static_dielectric_tensor,
                                    self.MgO_unitcell.static_dielectric_tensor)
        self.unitcell.ionic_dielectric_tensor = [9.102401, 9.102448, 9.102542]
        self.assertArrayAlmostEqual(self.unitcell.ionic_dielectric_tensor,
                                    self.MgO_unitcell.ionic_dielectric_tensor)
        # test upper triangle matrix form
        self.unitcell.ionic_dielectric_tensor = [1, 2, 3, 4, 5, 6]
        self.assertArrayAlmostEqual(self.unitcell.ionic_dielectric_tensor,
                                    [[1, 4, 5], [4, 2, 6], [5, 6, 3]])

    def test_band_edge(self):
        self.unitcell.set_band_edge_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/structure_optimization"))
        self.assertArrayAlmostEqual(self.unitcell.band_edge,
                                    self.MgO_unitcell.band_edge)

    def test_dielectric_constant(self):
        self.unitcell.set_static_dielectric_tensor_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/dielectric_constants"))
        self.unitcell.set_ionic_dielectric_tensor_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/dielectric_constants"))
        self.assertArrayAlmostEqual(self.unitcell.static_dielectric_tensor,
                                    self.MgO_unitcell.static_dielectric_tensor)
        self.assertArrayAlmostEqual(self.unitcell.ionic_dielectric_tensor,
                                    self.MgO_unitcell.ionic_dielectric_tensor)
        MgO_total_dielectric_tensor = [[12.269128, 0, 0],
                                       [0, 12.269175, 0],
                                       [0, 0, 12.269269]]
        self.assertArrayAlmostEqual(self.unitcell.total_dielectric_tensor,
                                    MgO_total_dielectric_tensor)

    def test_total_dos(self):
        self.unitcell.set_total_dos_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/structure_optimization"))
        # check length of dos
        self.assertEqual(len(self.unitcell.total_dos[0]), 301)
        # check 301st density of states
        self.assertAlmostEqual(self.unitcell.total_dos[1][300],
                               self.MgO_unitcell.total_dos[1][300])

    def test_volume(self):
        self.unitcell.set_volume_from_vasp(
            directory_path=os.path.join(test_dir,
                                        "MgO/unitcell/structure_optimization"))
        self.assertAlmostEqual(self.unitcell.volume, self.MgO_unitcell.volume)

    def test_dict_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self.MgO_unitcell.to_json_file(tmp_file.name)
        unitcell_from_json = UnitcellCalcResults.load_json(tmp_file.name)
        self.assertEqual(unitcell_from_json.as_dict,
                         self.MgO_unitcell.as_dict())

    def test_print(self):
        print(self.MgO_unitcell)


if __name__ == "__main__":
    unittest.main()

