import unittest
from pydefect.taka_my_classes.perfect import Perfect
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
PERFECT_DIR = "NEED_TO_SPECIFY!"
DEFECT_JSON = "NEED_TO_SPECIFY!"
DEFECT_ENERGY = "NEED_TO_SPECIFY!"
DEFECT_ELECTROSTATIC_POTENTIAL = "NEED_TO_SPECIFY!"
DEFECT_FINAL_STRUCTURE = "NEED_TO_SPECIFY!"


class DefectCorrectionTest(unittest.TestCase):

    def setUp(self):
        self._perfect = Perfect(PERFECT_DIR)
        self._defect_input = Defect.json_load(DEFECT_JSON)
        pass

    def test_correct_energy(self):
        # TODO: write expected value.
        expected = "NEED_TO_SPECIFY!"
        actual = correct_energy(self._perfect, self._defect_input,
                                DEFECT_ENERGY,
                                DEFECT_ELECTROSTATIC_POTENTIAL,
                                DEFECT_FINAL_STRUCTURE,
                                method = "ExtendedFNV")
        self.assertAlmostEqual(actual, expected)

if __name__ == "__main__":
    unittest.main()
