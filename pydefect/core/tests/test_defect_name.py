# -*- coding: utf-8 -*-
import unittest

from pydefect.core.defect_name import DefectName

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class DefectNameTest(unittest.TestCase):

    def setUp(self):
        """
        """
        self.va_o_1 = DefectName("Va_O1", 1)
        self.mg_i1_m1 = DefectName.from_str("Mg_i1_1")

    def test_to_str(self):
        actual = str(self.va_o_1)
        expected = "Va_O1_1"
        self.assertEqual(actual, expected)

    def test_is_name_matched(self):
        self.assertTrue(self.va_o_1.is_name_matched("Va_O[0-9]_1"))
        self.assertFalse(self.va_o_1.is_name_matched(["Va_N[0-9]_1", "Va_Mg[0-9]"]))
