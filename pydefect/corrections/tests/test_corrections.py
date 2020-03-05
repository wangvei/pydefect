# -*- coding: utf-8 -*-
from pydefect.util.testing import PydefectTest

from pydefect.corrections.corrections import ManualCorrection


class ManualCorrectionTest(PydefectTest):
    def setUp(self):
        self.correction = ManualCorrection(manual_correction_energy=1.5)

    def test_dict(self):
        d = self.correction.as_dict()
        actual = ManualCorrection.from_dict(d).as_dict()
        expected = d
        self.assertEqual(actual, expected)

    def test_msonable(self):
        self.assertMSONable(self.correction)
