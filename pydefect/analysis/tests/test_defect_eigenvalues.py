# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.analysis.defect_eigenvalues import DefectEigenvalue
from pydefect.corrections.corrections import ExtendedFnvCorrection
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_entry import DefectEntry
from pydefect.analysis.defect_energies import Defect

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectEigenvalueTest(unittest.TestCase):

    def setUp(self):
        """ """
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        unitcell = UnitcellCalcResults.load_json(unitcell_file)
        perfect_file = os.path.join(test_dir,
                                    "MgO/defects/perfect/dft_results.json")
        perfect = SupercellCalcResults.load_json(perfect_file)

        d = os.path.join(test_dir, "MgO/defects/Mg_i1_1")
#        d = os.path.join(test_dir, "MgO/defects/Va_O1_2")
        defect_entry = \
            DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
        dft_results = \
            SupercellCalcResults.load_json(
                os.path.join(d, "dft_results.json"))
        correction = \
            ExtendedFnvCorrection.load_json(os.path.join(d, "correction.json"))

        defect = Defect(defect_entry=defect_entry,
                        dft_results=dft_results,
                        correction=correction)

        self.defect_eigenvalues = DefectEigenvalue.from_files(unitcell=unitcell,
                                                              perfect=perfect,
                                                              defect=defect)

    def test_diagnose_shallow(self):
        print(self.defect_eigenvalues.diagnose_shallow_states())

    def test(self):
        print(self.defect_eigenvalues.vbm)

    def test_plot(self):
        self.defect_eigenvalues.plot()


if __name__ == "__main__":
    unittest.main()

