# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.analysis.defect_structure import DefectStructure
from pydefect.corrections.corrections import ExtendedFnvCorrection
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.defect_entry import DefectEntry
from pydefect.analysis.defect import Defect

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectStructureTest(unittest.TestCase):

    def setUp(self):
        """ """
        d = os.path.join(test_dir, "MgO/defects", "Va_O1_0")
        defect_entry = \
            DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
        dft_results = \
            SupercellCalcResults.load_json(
                os.path.join(d, "dft_results.json"))
        correction = \
            ExtendedFnvCorrection.load_json(os.path.join(d, "correction.json"))

        self.va_o_0 = DefectStructure.from_defect(
            Defect(defect_entry=defect_entry,
                   dft_results=dft_results,
                   correction=correction))

        d = os.path.join(test_dir, "MgO/defects", "Va_O1_2")
        defect_entry = \
            DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
        dft_results = \
            SupercellCalcResults.load_json(
                os.path.join(d, "dft_results.json"))
        correction = \
            ExtendedFnvCorrection.load_json(os.path.join(d, "correction.json"))

        self.va_o_2 = DefectStructure.from_defect(
            Defect(defect_entry=defect_entry,
                   dft_results=dft_results,
                   correction=correction))

#    def test(self):
#        print(self.va_o_0.initial_local_structure)
#        print(self.va_o_2.show_displacements)

    def test2(self):
        print(self.va_o_0.final_structure)
        s = self.va_o_2.final_local_structure
        print(s)
        print(self.va_o_0.final_local_structure)
#        print(self.va_o_0.comparator(defect_local_structure=s))


if __name__ == "__main__":
    unittest.main()
