# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.analysis.defect import Defect
from pydefect.analysis.defect_eigenvalues import DefectEigenvalue
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection
from pydefect.util.testing import PydefectTest


class DefectEigenvalueTest(PydefectTest):

    def setUp(self):
        """ """
        filename = ["defects", "MgO", "unitcell.json"]
        unitcell = self.get_object_by_name(
            UnitcellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "perfect", "dft_results.json"]
        perfect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "Mg_i1_2", "dft_results.json"]
        dft_results = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "Mg_i1_2", "defect_entry.json"]
        defect_entry = self.get_object_by_name(
            DefectEntry.load_json, filename)

        filename = ["defects", "MgO", "Mg_i1_2", "correction.json"]
        correction = self.get_object_by_name(
            ExtendedFnvCorrection.load_json, filename)

        defect = Defect.from_objects(defect_entry=defect_entry,
                                     dft_results=dft_results,
                                     perfect_dft_results=perfect,
                                     correction=correction)

        self.defect_eigenvalues = DefectEigenvalue.from_files(unitcell=unitcell,
                                                              defect=defect)

    def test(self):
        print(self.defect_eigenvalues.vbm)

    def test_plot(self):
        self.defect_eigenvalues.plot()

