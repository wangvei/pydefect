# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import det
from pydefect.database.symmetry import (
    tm_from_standard_to_primitive, tm_from_primitive_to_standard)
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class TmFromPrimitiveToStandardTest(PydefectTest):
    def test(self):
        expected = np.array([[ 1,  0,  0],
                             [ 0,  1,  1],
                             [ 0, -1,  1]], dtype=int)
        actual = tm_from_primitive_to_standard("A")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 1, -1,  0],
                             [ 1,  1,  0],
                             [ 0,  0,  1]], dtype=int)
        actual = tm_from_primitive_to_standard("C")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 1,  0,  1],
                             [-1,  1,  1],
                             [ 0, -1,  1]], dtype=int)
        actual = tm_from_primitive_to_standard("R")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 0,  1,  1],
                             [ 1,  0,  1],
                             [ 1,  1,  0]], dtype=int)
        actual = tm_from_primitive_to_standard("I")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[-1,  1,  1],
                             [ 1, -1,  1],
                             [ 1,  1, -1]], dtype=int)
        actual = tm_from_primitive_to_standard("F")
        self.assertArrayEqual(expected, actual)


