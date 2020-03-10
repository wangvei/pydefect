# -*- coding: utf-8 -*-

from argparse import Namespace
import os
from pathlib import Path
from unittest.mock import patch
from shutil import copyfile, rmtree
from pydefect.util.testing import PydefectTest

from pydefect.cli.main_functions import (
    unitcell_calc_results, initial_setting
)


parent_dir = Path(__file__).parent


class UnitcellCalcResultsTest(PydefectTest):
    @patch("pydefect.core.unitcell_calc_results.UnitcellCalcResults.load_json")
    def test_print(self, mock):
        args_print = Namespace(print=True, json_file="unitcell.json")
        unitcell_calc_results(args_print)
        mock.assert_called_once_with(filename="unitcell.json")

    def test_load_json(self):
        d = Path("band_edge_dir")
        d.mkdir()
        outcar = self.TEST_FILES_DIR / "MgAl2Se4_unitcell" / "band_OUTCAR"
        vasprun = self.TEST_FILES_DIR / "MgAl2Se4_unitcell" / "band_vasprun.xml"
        copyfile(outcar, d / "OUTCAR")
        copyfile(vasprun, d / "vasprun.xml")

        d = Path("diele")
        d.mkdir()
        outcar = self.TEST_FILES_DIR / "MgAl2Se4_unitcell" / "dielectric_OUTCAR"
        copyfile(outcar, d / "OUTCAR")

        d = Path("dos")
        d.mkdir()
        vasprun = self.TEST_FILES_DIR / "MgAl2Se4_unitcell" / "dos_vasprun.xml"
        copyfile(vasprun, d / "vasprun.xml")

        args = Namespace(print=False,
                         static_diele=None,
                         ionic_diele=None,
                         band_edge_dir="band_edge_dir",
                         static_diele_dir="diele",
                         ionic_diele_dir="diele",
                         total_dos_dir="dos",
                         outcar="OUTCAR",
                         vasprun="vasprun.xml",
                         json_file="unitcell.json")
        unitcell_calc_results(args)

    def tearDown(self) -> None:
        for f in ["band_edge_dir", "diele", "dos"]:
            try:
                rmtree(f)
            except FileNotFoundError:
                pass

        try:
            Path("unitcell.json").unlink()
        except FileNotFoundError:
            pass


class InitialSettingTest(PydefectTest):
    def setUp(self) -> None:
        self.kwargs = {
            "dopants": ["Ga", "In"],
            "antisite": False,
            "en_diff": 4.1,
            "included": None,
            "excluded": None,
            "displacement_distance": 0.1,
            "symprec": 0.1,
            "angle_tolerance": 0.1,
            "interstitials": ["all"],
            "complex_defect_names": [],
        }

        self.kwargs_supercells = {
            "isotropy_criterion": 0.1,
            "min_num_atoms": 40,
            "conventional_base": False,
            "most_isotropic": True,
            "rhombohedral_angle": 80}

    def test_print(self):
        initial_setting(Namespace(print_dopant="H"))

    def test_matrix(self):
        args = Namespace(print_dopant=None,
                         matrix=[2, 2, 2],
                         poscar="POSCAR",
                         **self.kwargs)
        copyfile(self.POSCARS_DIR / "POSCAR-MgO", "POSCAR")
        self.assertIsNone(initial_setting(args))

    def test_no_supercell(self):
        args = Namespace(print_dopant=None,
                         matrix=None,
                         poscar="POSCAR",
                         supercell_set=False,
                         max_num_atoms=40,
                         **self.kwargs,
                         **self.kwargs_supercells)
        copyfile(self.POSCARS_DIR / "POSCAR-MgO", "POSCAR")
        self.assertFalse(initial_setting(args))

    def test_normal(self):
        args = Namespace(print_dopant=None,
                         matrix=None,
                         poscar="POSCAR",
                         supercell_set=False,
                         max_num_atoms=100,
                         **self.kwargs,
                         **self.kwargs_supercells)
        copyfile(self.POSCARS_DIR / "POSCAR-MgO", "POSCAR")
        self.assertIsNone(initial_setting(args))

    def test_set(self):
        args = Namespace(print_dopant=None,
                         matrix=None,
                         poscar="POSCAR",
                         supercell_set=True,
                         max_num_atoms=100,
                         **self.kwargs,
                         **self.kwargs_supercells)
        copyfile(self.POSCARS_DIR / "POSCAR-MgO", "POSCAR")
        self.assertIsNone(initial_setting(args))

    def tearDown(self) -> None:
        for f in ["defect.in", "DPOSCAR", "POSCAR"]:
            try:
                Path(f).unlink()
            except FileNotFoundError:
                pass

        try:
            rmtree("p3x3x3_54_0.0")
        except FileNotFoundError:
            pass



    # @patch("pydefect.input_maker.defect_initial_setting.DefectInitialSetting.from_basic_settings")
    # @patch("pydefect.input_maker.supercell_maker.Supercell")
    # @patch("pymatgen.core.structure.Structure.from_file")
    # def test_matrix(self, mock_structure, mock_supercell, mock_initial_setting):
    #     initial_setting(self.args_matrix)
    #     mock_structure.assert_called_once_with(self.kwargs["poscar"])
    #     mock_supercell.assert_called_once_with(
    #         structure=mock_structure,
    #         trans_mat=[10],
    #         check_unitcell=True,
    #         symprec=self.kwargs["symprec"],
    #         angle_tolerance=self.kwargs["angle_tolerance"])

        # mock_initial_setting.assert_called_once_with(
        #     structure=mock_structure,
        #     transformation_matrix=mock_structure.trans_mat_tolist,
        #     cell_multiplicity=mock_structure.multiplicity,
        #     **self.kwargs)






