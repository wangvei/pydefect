# -*- coding: utf-8 -*-
from pathlib import Path

from pydefect.core.defect_entry import DefectEntry
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection
from pydefect.util.main_tools import (
    get_default_args, generate_objects_from_json_files)
from vise.util.main_tools import potcar_str2dict, list2dict, get_user_settings, \
    dict2list
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class GetUserSettingsTest(PydefectTest):
    def setUp(self) -> None:
        self.user_setting = get_user_settings()

    def test_dict(self):
        actual = self.user_setting['ldauu']
        expected = {'Mn': 5, 'O': 3}
        self.assertEqual(expected, actual)

    def test_bool(self):
        actual = self.user_setting['perfect_incar_setting']['ALGO']
        self.assertEqual(False, actual)

    def test_float(self):
        actual = self.user_setting['symprec']
        self.assertEqual(0.01, actual)

    def test_string(self):
        actual = self.user_setting['xc']
        self.assertEqual("hse", actual)

    def test_string_list(self):
        actual = self.user_setting['potcar_set']
        self.assertEqual("Mg_pv O_h", actual)

    def test_fullpath(self):
        actual = self.user_setting['volume_dir']
        cwd = Path.cwd()
        expected = cwd / "test" / "unitcell" / "structure_opt"
        self.assertEqual(str(expected), actual)

    def test_fail(self):
        with self.assertRaises(ValueError):
            get_user_settings("pydefect_fail.yaml")


class GetDefaultArgsTest(PydefectTest):
    def setUp(self) -> None:
        def test_func(a, b=1, c=True):
            pass

        class test_cls:
            def __init__(self, d, e=2, f=None):
                pass
            def test_method(self, g, h=3.0):
                pass

        self.default_func = get_default_args(test_func)
        self.default_cls = get_default_args(test_cls)
        self.default_method = get_default_args(test_cls.test_method)

    def test_func(self):
        expected = {"b": 1, "c": True}
        self.assertEqual(expected, self.default_func)

    def test_cls(self):
        expected = {"e": 2, "f": None}
        self.assertEqual(expected, self.default_cls)

    def test_method(self):
        expected = {"h": 3.0}
        self.assertEqual(expected, self.default_method)


class PotcarStr2DictTest(PydefectTest):
    def test_potcar_str2dict(self):
        actual = potcar_str2dict("Mg_pv O_h Se")
        expected = {"Mg": "Mg_pv", "O": "O_h", "Se": "Se"}
        self.assertEqual(expected, actual)


class Dict2ListTest(PydefectTest):
    def test_dict2list(self):
        actual = dict2list({"a": 1, "b": "2 3 4", "c": True})
        expected = ["a", "1", "b", "2", "3", "4", "c", "True"]
        self.assertEqual(expected, actual)


class List2DictTest(PydefectTest):
    def setUp(self) -> None:
        self.key_candidates = ["ENCUT", "MAGMOM", "LWAVE"]

    def test_dict2list(self):
        flattened_list = ["ENCUT", "500", "MAGMOM", "4", "4.0", "LWAVE", "F"]
        actual = list2dict(flattened_list, self.key_candidates)
        expected = {"ENCUT": 500, "MAGMOM": [4, 4.0], "LWAVE": False}
        self.assertEqual(expected, actual)

    def test_fail(self):
        flattened_list = ["ENCAT", "500"]
        with self.assertRaises(ValueError):
            list2dict(flattened_list, self.key_candidates)

    def test_fail2(self):
        flattened_list = ["ENCUT", "500", "MAGMOM"]
        with self.assertRaises(ValueError):
            list2dict(flattened_list, self.key_candidates)


class GenerateObjectsFromJsonFilesTest(PydefectTest):
    def test_success(self) -> None:
        directory = self.DEFECTS_MGO_DIR / "Va_O1_0"
        filenames = ["defect_entry.json", "correction.json"]
        classes = [DefectEntry, ExtendedFnvCorrection]
        self.objects = generate_objects_from_json_files(
            directory, filenames, classes)
        self.assertTrue(isinstance(self.objects[0], DefectEntry))
        self.assertTrue(isinstance(self.objects[1], ExtendedFnvCorrection))
