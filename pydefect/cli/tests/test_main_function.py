# -*- coding: utf-8 -*-

from argparse import Namespace
from unittest.mock import patch

from pydefect.util.testing import PydefectTest

from pydefect.cli.main_functions import (
    unitcell_calc_results, initial_setting
)


class Test(PydefectTest):
    def setUp(self) -> None:
        self.args_print = Namespace(print=True,
                                    json_file="unitcell.json")

        # self.args_num = Namespace(number=110, poscar="POSCAR")
        # self.kwargs = {"elements": ["Mg", "O"],
        #                "e_above_hull": 1.5,
        #                "molecules": True}
        # self.args_elements = Namespace(**self.kwargs)
        # self.args_none = Namespace()

    @patch("pydefect.core.unitcell_calc_results.UnitcellCalcResults.load_json")
    def test_print(self, mock):
        unitcell_calc_results(self.args_print)
        mock.assert_called_once_with(filename="unitcell.json")


class InitialSettingTest(PydefectTest):
    def setUp(self) -> None:
        self.args_print_dopant = Namespace(print_dopant="H")

        self.kwargs = {
            "dopants": ["Ga", "In"],
            "antisite": False,
            "en_diff": 4.1,
            "included": ["b"],
            "excluded": ["c"],
            "displacement_distance": 5.1,
            "symprec": 6.1,
            "angle_tolerance": 7.1,
            "interstitials": ["i1", "i2"],
            "complex_defect_names": ["d", "e"],
        }

        self.args_matrix = Namespace(print_dopant=None,
                                     matrix=[10],
                                     poscar="a",
                                     **self.kwargs)

        # self.args_matrix = Namespace(
        #     isotropy_criterion=0.1,
        #     min_num_atoms=1,
        #     max_num_atoms=2,
        #     conventional_base=False,
        #     most_isotropic=True,
        #     rhombohedral_angle=3.1,
        #     supercell_set=True,
        #     print_dopant="Cs",

    def test_print(self):
        initial_setting(self.args_print_dopant)

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






