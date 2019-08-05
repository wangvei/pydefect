# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from chempotdiag.chem_pot_diag import ChemPotDiag

from pymatgen.util.testing import PymatgenTest

from pydefect.analysis.defect_energies import convert_str_in_dict
from pydefect.analysis.defect import Defect
from pydefect.corrections.corrections import ExtendedFnvCorrection
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core", "MgO", "defects")


class DefectEnergiesTest(PymatgenTest):

    def setUp(self):
        """ """
        perfect_file = os.path.join(test_dir, "perfect", "dft_results.json")
        perfect = SupercellCalcResults.load_json(perfect_file)

        defect_dir = os.path.join(test_dir, "Va_O1_2")
        defect_entry_file = os.path.join(defect_dir, "defect_entry.json")
        defect_entry = DefectEntry.load_json(defect_entry_file)
        dft_results_file = os.path.join(defect_dir, "dft_results.json")
        dft_results = SupercellCalcResults.load_json(dft_results_file)
        correction_file = os.path.join(defect_dir, "correction.json")
        correction = ExtendedFnvCorrection.load_json(correction_file)

        self.defect = Defect.from_objects(defect_entry=defect_entry,
                                          dft_results=dft_results,
                                          perfect_dft_results=perfect,
                                          correction=correction)

    def test_dict(self):
        """ round trip test of to_json and from_json """
        defect_from_dict = Defect.from_dict(self.defect.as_dict())
        print(defect_from_dict.as_dict())
        self.assertEqual(defect_from_dict.as_dict(), self.defect.as_dict())

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.defect.to_json_file(tmp_file.name)
        defect_from_json = Defect.load_json(tmp_file.name)
        print(defect_from_json.as_dict())
        print(self.defect.as_dict())
        self.assertEqual(defect_from_json.as_dict(), self.defect.as_dict())


    # def test_print(self):
    #     print(self.defect_energies)

    # def test_energies(self):
    #     d = self.defect_energies.as_dict()
    #     de = DefectEnergies.from_dict(d)
    #     dd = de.as_dict()
    #     self.assertEqual(d, dd)
    #     print(d)
    #     print(dd)


if __name__ == "__main__":
    unittest.main()

