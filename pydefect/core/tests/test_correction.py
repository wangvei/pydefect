import unittest
import os
from pydefect.core.correction import Ewald, Correction, CorrectionMethod
from pydefect.core.dft_results import SupercellDftResults, UnitcellDftResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 20, 2018"

# TODO: write test directory.
test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core", "MgO")


dirname_unitcell = test_dir + "/unitcell/structure_optimization"
dirname_dielectric = test_dir + "/unitcell/dielectric_constants"

dirname_perfect = test_dir + "/defects/perfect"
# interstitial directory

vacancy_expected_alignment = "NEED_TO_SPECIFY!"
dirname_vacancy = test_dir + "/defects/Va_O1_2/"
vac_defect_entry_json = dirname_vacancy + "/defect_entry.json"


class EwaldTest(unittest.TestCase):

    def setUp(self):
        unitcell = UnitcellDftResults()
        unitcell.set_static_dielectric_tensor_from_vasp(dirname_dielectric)
        unitcell.set_ionic_dielectric_tensor_from_vasp(dirname_dielectric)
        perfect = SupercellDftResults.from_vasp_files(dirname_perfect)
        self._structure = perfect.final_structure
        self._dielectric_tensor = unitcell.total_dielectric_tensor

    def test_optimize(self):
        ewald =\
            Ewald.from_optimization(self._structure, self._dielectric_tensor)
        print(ewald.ewald_param)
        print(ewald.num_real_lattice)
        print(ewald.num_reciprocal_lattice)
        print(ewald.dielectric_tensor)
        print(ewald.num_real_lattice)


class CorrectionTest(unittest.TestCase):

    def setUp(self):
        self._unitcell = UnitcellDftResults()
        self._unitcell.set_static_dielectric_tensor_from_vasp(
            dirname_dielectric)
        self._unitcell.set_ionic_dielectric_tensor_from_vasp(dirname_dielectric)
        self._perfect = SupercellDftResults.from_vasp_files(dirname_perfect)
        self._vacancy_entry = DefectEntry.json_load(vac_defect_entry_json)
        self._vacancy = SupercellDftResults.from_vasp_files(dirname_vacancy)
        self._structure = self._perfect.final_structure
        self._dielectric_tensor = self._unitcell.total_dielectric_tensor
        self._ewald = Ewald(
            lattice_matrix=self._structure.lattice.matrix,
            dielectric_tensor=self._unitcell.total_dielectric_tensor,
            ewald_param=0.3314065604136566,
            num_real_lattice=11051,
            num_reciprocal_lattice=10692)
        # self._ewald = \
        #     Ewald.from_optimization(self._structure, self._dielectric_tensor)
        print("setUp completed")

    def test_compute_extended_fnv(self):
        self._vacancy_correction = \
            Correction.compute_correction(self._vacancy_entry,
                                          self._vacancy,
                                          self._perfect,
                                          self._unitcell)
        # TODO: write expected value.
        # vacancy
        # self.assertAlmostEqual(actual_vacancy,
        #                        vacancy_expected_alignment)

    # def test_yaml(self):
    #     ewald = Ewald()
    #     correction = Correction(CorrectionMethod.extended_fnv,
    #                             )


if __name__ == "__main__":
    unittest.main()
