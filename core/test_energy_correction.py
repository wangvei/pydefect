import unittest
from pydefect.core.energy_correction import compute_alignment_by_extended_fnv
from pydefect.core.unitcell import Unitcell
from pydefect.core.supercell import Perfect, Defect

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 20, 2018"

# TODO: write test directory.
UNITCELL_DIRECTORY = "NEED_TO_SPECIFY!"
SUPERCELL_DIRECTORY = "NEED_TO_SPECIFY!"
# interstitial directory

VACANCY_JSON = "NEED_TO_SPECIFY!"
VACANCY_EXPECTED_ALIGNMENT = "NEED_TO_SPECIFY!"

SUBSTITUTIONAL_JSON = "NEED_TO_SPECIFY!"
SUBSTITUTIONAL_EXPECTED_ALIGNMENT = "NEED_TO_SPECIFY!"

INTERSTITIAL_JSON = "NEED_TO_SPECIFY!"
INTERSTITIAL_EXPECTED_ALIGNMENT = "NEED_TO_SPECIFY!"

EWALD_PARAM = "NEED_TO_SPECIFY!"


class DefectCorrectionTest(unittest.TestCase):

    def setUp(self):
        self._unitcell = Unitcell(UNITCELL_DIRECTORY)
        self._supercell = Perfect(SUPERCELL_DIRECTORY)
        self._vacancy = Defect.json_load(VACANCY_JSON)
        self._substitutional = Defect.json_load(SUBSTITUTIONAL_JSON)
        self._interstitial = Defect.json_load(INTERSTITIAL_JSON)
        self._ewald_param = "NEED_TO_SPECIFY!"

    def test_correct_energy_extended_fnv(self):
        # TODO: write expected value.
        # vacancy
        actual_vacancy = compute_alignment_by_extended_fnv(
            self._unitcell.dielectric_tensor,
            self._ewald_param,
            self._supercell.structure,
            self._supercell.electrostatic_potential,
            self._vacancy.structure,
            self._vacancy.electrostatic_potential,
            self._vacancy.defect_coords,
            self._vacancy.atomic_position_without_defect,
            self._vacancy.distance_from_defect,
            self._vacancy.charge)
        self.assertAlmostEqual(actual_vacancy,
                               VACANCY_EXPECTED_ALIGNMENT)

        # substitutional
        actual_substitutional = compute_alignment_by_extended_fnv(
            self._unitcell.dielectric_tensor,
            self._ewald_param,
            self._supercell.structure,
            self._supercell.electrostatic_potential,
            self._substitutional.structure,
            self._substitutional.electrostatic_potential,
            self._substitutional.defect_coords,
            self._substitutional.atomic_position_without_defect,
            self._substitutional.distance_from_defect,
            self._substitutional.charge)
        self.assertAlmostEqual(actual_substitutional,
                               SUBSTITUTIONAL_EXPECTED_ALIGNMENT)

        # interstitial
        actual_interstitial = compute_alignment_by_extended_fnv(
            self._unitcell.dielectric_tensor,
            self._ewald_param,
            self._supercell.structure,
            self._supercell.electrostatic_potential,
            self._interstitial.structure,
            self._interstitial.electrostatic_potential,
            self._interstitial.defect_coords,
            self._interstitial.atomic_position_without_defect,
            self._interstitial.distance_from_defect,
            self._interstitial.charge)
        self.assertAlmostEqual(actual_interstitial,
                               INTERSTITIAL_EXPECTED_ALIGNMENT)


if __name__ == "__main__":
    unittest.main()