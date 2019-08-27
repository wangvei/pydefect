# -*- coding: utf-8 -*-

from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class UnitcellDftResultsTest(PydefectTest):

    def setUp(self):
        self.unitcell_dir = self.TEST_FILES_DIR / "MgAl2Se4_unitcell"

        self.unitcell = UnitcellCalcResults(band_edge=None,
                                            static_dielectric_tensor=None,
                                            ionic_dielectric_tensor=None,
                                            total_dos=None,
                                            volume=None)

    def test(self):
        self.unitcell.set_static_dielectric_tensor_from_vasp(
            directory_path=self.unitcell_dir, outcar_name="dielectric_OUTCAR")
        expected = [[6.072649, 0, 0], [0, 6.072649, 0], [0, 0, 5.643136]]
        actual = self.unitcell.static_dielectric_tensor
        self.assertArrayAlmostEqual(expected, actual, 5)

        self.unitcell.set_ionic_dielectric_tensor_from_vasp(
            directory_path=self.unitcell_dir, outcar_name="dielectric_OUTCAR")
        expected = [[3.769646, 0, 0], [0, 3.769646, 0], [0, 0, 1.239479]]
        actual = self.unitcell.ionic_dielectric_tensor
        self.assertArrayAlmostEqual(expected, actual, 5)

        self.unitcell.set_volume_from_vasp(
            directory_path=self.unitcell_dir, contcar_name="CONTCAR")
        expected = 168.72742284380834
        actual = self.unitcell.volume
        self.assertEqual(expected, actual)

        self.unitcell.set_band_edge_from_vasp(
            directory_path=self.unitcell_dir, vasprun_name="band_vasprun.xml")
        expected = [3.127, 4.197]
        actual = self.unitcell.band_edge
        self.assertArrayAlmostEqual(expected, actual, 3)

        self.unitcell.set_total_dos_from_vasp(
            directory_path=self.unitcell_dir, vasprun_name="dos_vasprun.xml")

        self.assertMSONable(self.unitcell)
