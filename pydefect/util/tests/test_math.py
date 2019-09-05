# -*- coding: utf-8 -*-
import numpy as np
from pydefect.util.math import normalized_random_3d_vector
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class NormalizedRandom3dVectorTest(PydefectTest):

    def setUp(self):
        self.v = normalized_random_3d_vector()

    def test(self):
        print("normalized_random_3d_vector: ", self.v)
        print("norm: ", np.linalg.norm(self.v))
        self.assertAlmostEqual(np.linalg.norm(self.v), 1.0)


