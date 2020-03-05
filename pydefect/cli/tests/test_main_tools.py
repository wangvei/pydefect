# -*- coding: utf-8 -*-

from pydefect.core.defect_entry import DefectEntry
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection
from pydefect.cli.main_tools import generate_objects_from_json_files
from pydefect.util.testing import PydefectTest


class GenerateObjectsFromJsonFilesTest(PydefectTest):
    def test_success(self) -> None:
        directory = self.DEFECTS_MGO_DIR / "Va_O1_0"
        filenames = ["defect_entry.json", "correction.json"]
        classes = [DefectEntry, ExtendedFnvCorrection]
        self.objects = generate_objects_from_json_files(
            directory, filenames, classes)
        self.assertTrue(isinstance(self.objects[0], DefectEntry))
        self.assertTrue(isinstance(self.objects[1], ExtendedFnvCorrection))
