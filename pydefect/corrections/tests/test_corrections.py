from pydefect.util.testing import PydefectTest

from pydefect.corrections.corrections import ManualCorrection

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class NoCorrectionTest(PydefectTest):
    def setUp(self):
        self._correction = ManualCorrection(manual_correction_energy=1.5)

    def test_dict(self):
        d = self._correction.as_dict()
        actual = ManualCorrection.from_dict(d).as_dict()
        expected = d
        self.assertEqual(actual, expected)


