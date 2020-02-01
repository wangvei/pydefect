from pydefect.util.testing import PydefectTest
from pymatgen.core.structure import Structure
from pydefect.corrections.vertical_transition_energy_correction import \
    VerticalTransitionEnergyCorrection
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class VerticalTransitionEnergyCorrectionTest(PydefectTest):
    def setUp(self) -> None:
        unitcell = UnitcellCalcResults.load_json("MgO_unitcell.json")
        initial_efnv_cor = \
            ExtendedFnvCorrection.load_json("MgO_Va_O_1_correction.json")
        initial_calc_results = \
            SupercellCalcResults.load_json("MgO_Va_O_1_dft_results.json")
        final_defect_entry = \
            DefectEntry.load_json("MgO_Va_O_1-added_defect_entry.json")
        final_calc_results = \
            SupercellCalcResults.load_json("MgO_Va_O_1-added_dft_results.json")

        self.vertical_correction = \
            VerticalTransitionEnergyCorrection.from_files(unitcell,
                                                          initial_efnv_cor,
                                                          initial_calc_results,
                                                          final_defect_entry,
                                                          final_calc_results,
                                                          8)

    def test(self):
        print(self.vertical_correction.correction_energy)
        print(self.vertical_correction.pc_correction_energy)
        print(self.vertical_correction.alignment_correction_energy)

