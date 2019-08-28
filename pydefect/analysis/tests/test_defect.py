# -*- coding: utf-8 -*-

from pydefect.analysis.defect import (
    BandEdgeState)

from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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
    pass
    # def setUp(self):
    #     """ """
    #     perfect_file = os.path.join(test_dir, "perfect", "dft_results.json")
    #     perfect = SupercellCalcResults.load_json(perfect_file)

        # defect_dir = os.path.join(test_dir, "Va_O1_2")
        # defect_entry_file = os.path.join(defect_dir, "defect_entry.json")
        # defect_entry = DefectEntry.load_json(defect_entry_file)
        # dft_results_file = os.path.join(defect_dir, "dft_results.json")
        # dft_results = SupercellCalcResults.load_json(dft_results_file)
        # correction_file = os.path.join(defect_dir, "correction.json")
        # correction = ExtendedFnvCorrection.load_json(correction_file)

        # self.defect = Defect.from_objects(defect_entry=defect_entry,
        #                                   dft_results=dft_results,
        #                                   perfect_dft_results=perfect,
        #                                   correction=correction)

    # def test_dict(self):
    #     """ round trip test of to_json and from_json """
    #     defect_from_dict = Defect.from_dict(self.defect.as_dict())
    #     print(defect_from_dict.as_dict())
    #     self.assertEqual(defect_from_dict.as_dict(), self.defect.as_dict())

#     def test_json(self):
#         """ round trip test of to_json and from_json """
# #        tmp_file = os.path.join(test_dir, "Va_O1_2")
#         tmp_file = tempfile.NamedTemporaryFile()
#         self.defect.to_json_file(tmp_file.name)
#         defect_from_json = Defect.load_json(tmp_file.name)
#         print(defect_from_json.as_dict())
#         print(self.defect.as_dict())
#         self.assertEqual(defect_from_json.as_dict(), self.defect.as_dict())


#     # def test_print(self):
#     #     print(self.defect_energies)

#     # def test_energies(self):
#     #     d = self.defect_energies.as_dict()
#     #     de = DefectEnergies.from_dict(d)
#     #     dd = de.as_dict()
#     #     self.assertEqual(d, dd)
#     #     print(d)
#     #     print(dd)

