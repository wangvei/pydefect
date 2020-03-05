# -*- coding: utf-8 -*-

from pydefect.analysis.defect import Defect
from pydefect.analysis.defect_structure import (
    DefectStructure, defect_structure_matcher)
from pydefect.util.testing import PydefectTest


class DefectStructureTest(PydefectTest):

    def setUp(self):
        filename = ["defects", "MgO", "Va_O1_0", "defect.json"]
        d = self.get_object_by_name(Defect.load_json, filename)
        self.va_o_0 = DefectStructure.from_defect(d)

        filename = ["defects", "MgO", "Va_O1_2", "defect.json"]
        d = self.get_object_by_name(Defect.load_json, filename)
        self.va_o_2 = DefectStructure.from_defect(d)

        filename = ["defects", "MgO", "Mg_i1_0", "defect.json"]
        d = self.get_object_by_name(Defect.load_json, filename)
        self.mg_i_0 = DefectStructure.from_defect(d)

        filename = ["defects", "MgO", "Mg_i1_1", "defect.json"]
        d = self.get_object_by_name(Defect.load_json, filename)
        self.mg_i_1 = DefectStructure.from_defect(d)

    def test2(self):
        print(self.va_o_2.final_structure)
        print(self.va_o_2.final_local_structure)
        print(self.va_o_2.show_displacements())
#        print(self.va_o_0.comparator(defect_local_structure=s))

    def test3(self):
        print(self.mg_i_0.show_displacements())
        print(self.mg_i_1.show_displacements())
        print(self.mg_i_0.comparator(self.mg_i_1))
        print(self.va_o_2.comparator(self.va_o_0))

    def test(self):
        group = defect_structure_matcher(
            [self.va_o_0, self.va_o_2, self.mg_i_1, self.mg_i_0])
        print(group)
