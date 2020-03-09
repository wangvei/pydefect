# -*- coding: utf-8 -*-
import os
from pathlib import Path

from pydefect.util.testing import PydefectTest
from pydefect.corrections.vertical_transition_energy_correction import \
    VerticalTransitionEnergyCorrection
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.defect_entry import DefectEntry

parent = Path(__file__).parent


class VerticalTransitionEnergyCorrectionTest(PydefectTest):
    def setUp(self) -> None:
        unitcell = UnitcellCalcResults.load_json(parent / "MgO_unitcell.json")
        dielectric_tensor = unitcell.total_dielectric_tensor
        static_dielectric_tensor = unitcell.static_dielectric_tensor
        initial_efnv_cor = ExtendedFnvCorrection.load_json(
            parent / "MgO_Va_O_1_correction.json")
        initial_calc_results = SupercellCalcResults.load_json(
            parent / "MgO_Va_O_1_dft_results.json")
        final_defect_entry = DefectEntry.load_json(
            parent / "MgO_Va_O_1-added_defect_entry.json")
        final_calc_results = SupercellCalcResults.load_json(
            parent / "MgO_Va_O_1-added_dft_results.json")

        self.vertical_correction = \
            VerticalTransitionEnergyCorrection.from_files(
                dielectric_tensor,
                static_dielectric_tensor,
                initial_efnv_cor,
                initial_calc_results,
                final_defect_entry,
                final_calc_results,
                8)

    def test(self):
        print(self.vertical_correction)
        self.vertical_correction.plot_potential()
        os.remove("vertical_correction.pdf")

