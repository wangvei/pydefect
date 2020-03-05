# -*- coding: utf-8 -*-

from pydefect.analysis.defect import (
    BandEdgeState, too_close_atom_pairs, Defect)
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.corrections.efnv_corrections import (
    ExtendedFnvCorrection)
from pydefect.util.testing import PydefectTest


class TooCloseAtomPairsTest(PydefectTest):
    def setUp(self) -> None:
        self.kzn4p3_struct = self.get_structure_by_name("KZn4P3")
        self.kzn4p3_struct.append("H", coords=[0.48, 0.48, 0.48],
                                  coords_are_cartesian=False)

    def test(self):
        self.assertTrue(too_close_atom_pairs(self.kzn4p3_struct))


class BandEdgeStateTest(PydefectTest):
    def setUp(self):
        self.donor_phs = BandEdgeState.donor_phs
        self.no_in_gap = BandEdgeState.no_in_gap

    def test_to_str(self):
        actual = str(self.donor_phs)
        expected = "Donor PHS"
        self.assertEqual(actual, expected)

    def test_shallow(self):
        self.assertTrue(self.donor_phs.is_shallow)
        self.assertFalse(self.no_in_gap.is_shallow)


class DiagnoseBandEdgesTest(PydefectTest):
    pass


class DefectTest(PydefectTest):
    def setUp(self) -> None:
        filename = ["defects", "MgO", "Va_O1_2", "dft_results.json"]
        dft_results = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "Va_O1_2", "defect_entry.json"]
        defect_entry = self.get_object_by_name(
            DefectEntry.load_json, filename)

        filename = ["defects", "MgO", "perfect", "dft_results.json"]
        perfect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        filename = ["defects", "MgO", "Va_O1_2", "correction.json"]
        correction = self.get_object_by_name(
            ExtendedFnvCorrection.load_json, filename)

        self.defect = Defect.from_objects(defect_entry=defect_entry,
                                          dft_results=dft_results,
                                          perfect_dft_results=perfect,
                                          correction=correction)

    def test_msonable(self):
        self.assertMSONable(self.defect)

    def test_dict(self):
        expected = self.defect.as_dict()
        actual = Defect.from_dict(expected).as_dict()
        self.assertEqual(expected, actual)


