import unittest
from pydefect.taka_my_classes.perfect import Unitcell
from pydefect.input_generator.defect import Defect
from pydefect.taka_my_classes.energy_correction import correct_energy

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
VACANCY_ENERGY = "NEED_TO_SPECIFY!"
VACANCY_ELECTROSTATIC_POTENTIAL = "NEED_TO_SPECIFY!"
VACANCY_FINAL_STRUCTURE = "NEED_TO_SPECIFY!"
VACANCY_EXPECTED_ENERGY = "NEED_TO_SPECIFY!"

SUBSTITUTIONAL_JSON = "NEED_TO_SPECIFY!"
SUBSTITUTIONAL_ENERGY = "NEED_TO_SPECIFY!"
SUBSTITUTIONAL_ELECTROSTATIC_POTENTIAL = "NEED_TO_SPECIFY!"
SUBSTITUTIONAL_FINAL_STRUCTURE = "NEED_TO_SPECIFY!"
SUBSTITUTIONAL_EXPECTED_ENERGY = "NEED_TO_SPECIFY!"

INTERSTITIAL_JSON = "NEED_TO_SPECIFY!"
INTERSTITIAL_ENERGY = "NEED_TO_SPECIFY!"
INTERSTITIAL_ELECTROSTATIC_POTENTIAL = "NEED_TO_SPECIFY!"
INTERSTITIAL_FINAL_STRUCTURE = "NEED_TO_SPECIFY!"
INTERSTITIAL_EXPECTED_ENERGY = "NEED_TO_SPECIFY!"


class DefectCorrectionTest(unittest.TestCase):

    def setUp(self):
        self._unitcell = Unitcell(UNITCELL_DIRECTORY)
        self._supercell = Unitcell(SUPERCELL_DIRECTORY)
        self._vacancy_input = Defect.json_load(VACANCY_JSON)
        self._substitutional_input = Defect.json_load(SUBSTITUTIONAL_JSON)
        self._interstitial_input = Defect.json_load(INTERSTITIAL_JSON)

    def test_correct_energy_extended_fnv(self):
        # TODO: write expected value.
        # vacancy
        actual_vacancy = correct_energy(self._unitcell,
                                        self._supercell,
                                        self._vacancy_input,
                                        VACANCY_ENERGY,
                                        VACANCY_ELECTROSTATIC_POTENTIAL,
                                        VACANCY_FINAL_STRUCTURE,
                                        method = "ExtendedFNV")
        self.assertAlmostEqual(actual_vacancy,
                               VACANCY_EXPECTED_ENERGY)

        # substitutional
        actual_substitutional \
            = correct_energy(self._unitcell,
                             self._supercell,
                             self._substitutional_input,
                             SUBSTITUTIONAL_ENERGY,
                             SUBSTITUTIONAL_ELECTROSTATIC_POTENTIAL,
                             SUBSTITUTIONAL_FINAL_STRUCTURE,
                             method = "ExtendedFNV")
        self.assertAlmostEqual(actual_substitutional,
                               SUBSTITUTIONAL_EXPECTED_ENERGY)

        # interstitial
        actual_interstitial \
            = correct_energy(self._unitcell,
                             self._supercell,
                             self._interstitial_input,
                             INTERSTITIAL_ENERGY,
                             INTERSTITIAL_ELECTROSTATIC_POTENTIAL,
                             INTERSTITIAL_FINAL_STRUCTURE,
                             method = "ExtendedFNV")
        self.assertAlmostEqual(actual_interstitial,
                               INTERSTITIAL_EXPECTED_ENERGY)


if __name__ == "__main__":
    unittest.main()
