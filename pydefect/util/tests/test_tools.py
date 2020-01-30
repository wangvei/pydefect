# -*- coding: utf-8 -*-
from collections import defaultdict

from pydefect.util.tools import (
    defaultdict_to_dict, flatten_dict, mod_defaultdict)
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class DefaultdictToDictTest(PydefectTest):
    def setUp(self):
        """ """
        self.d = {"a": defaultdict(dict)}
        self.d["a"]["b"] = {"c": 1}

    def test(self):
        print(self.d)
        print(defaultdict_to_dict(self.d))


class FlattenDictTest(PydefectTest):
    def setUp(self) -> None:
        d = {"a": defaultdict(dict)}
        d["a"]["b"] = {"c": 1}
        self.d = defaultdict_to_dict(d)

    def test(self):
        self.assertEqual([['a', {'b': {'c': 1}}]], flatten_dict(self.d, 1))
        self.assertEqual([['a', 'b', {'c': 1}]], flatten_dict(self.d, 2))
        self.assertEqual([['a', 'b', 'c', 1]], flatten_dict(self.d))


class ModDefaultDictTest(PydefectTest):
    def test(self):
        d = mod_defaultdict(depth=2)
        self.assertEqual(d[0][0], None)
        with self.assertRaises(TypeError):
            d[0][0][0]
