# -*- coding: utf-8 -*-
from collections import defaultdict

from pydefect.util.tools import defaultdict_to_dict
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class DefaultdictToDictTest(PydefectTest):
    def setUp(self):
        """ """
        self.d = defaultdict(lambda: defaultdict(dict))
        self.d[1][2][4] = 10

    def test(self):
        print(self.d)
        print(defaultdict_to_dict(self.d))