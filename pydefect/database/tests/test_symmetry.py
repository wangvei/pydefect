# -*- coding: utf-8 -*-
import numpy as np

from pydefect.database.symmetry import transmat_primitive2standard
from pydefect.util.testing import PydefectTest


class TransmatPrimitive2StandardTest(PydefectTest):
    def test(self):
        expected = np.array([[ 1,  0,  0],
                             [ 0,  1,  1],
                             [ 0, -1,  1]], dtype=int)
        actual = transmat_primitive2standard("A")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 1, -1,  0],
                             [ 1,  1,  0],
                             [ 0,  0,  1]], dtype=int)
        actual = transmat_primitive2standard("C")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 1, -1,  0],
                             [ 0,  1, -1],
                             [ 1,  1,  1]], dtype=int)
        actual = transmat_primitive2standard("R")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[ 0,  1,  1],
                             [ 1,  0,  1],
                             [ 1,  1,  0]], dtype=int)
        actual = transmat_primitive2standard("I")
        self.assertArrayEqual(expected, actual)

        expected = np.array([[-1,  1,  1],
                             [ 1, -1,  1],
                             [ 1,  1, -1]], dtype=int)
        actual = transmat_primitive2standard("F")
        self.assertArrayEqual(expected, actual)


