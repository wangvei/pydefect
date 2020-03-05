# -*- coding: utf-8 -*-
from pydefect.core.defect_name import DefectName
from pydefect.util.testing import PydefectTest


class DefectNameTest(PydefectTest):

    def setUp(self):
        self.va_o_1 = DefectName("Va_O1", 1)
        self.mg_i1_m1 = DefectName.from_str("Mg_i1_1")

    def test_to_str(self):
        actual = str(self.va_o_1)
        expected = "Va_O1_1"
        self.assertEqual(actual, expected)

    def test_is_name_matched(self):
        self.assertTrue(self.va_o_1.is_name_matched("Va_O[0-9]_1"))
        self.assertFalse(
            self.va_o_1.is_name_matched(["Va_N[0-9]_1", "Va_Mg[0-9]"]))

    def test_msonable(self):
        self.assertMSONable(self.va_o_1)
