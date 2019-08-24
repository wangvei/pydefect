# -*- coding: utf-8 -*-
import numpy as np
from pydefect.util.math import normalized_random_3d_vector, random_vector
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


class RandomVectorTest(PydefectTest):

    def setUp(self):
        self.distance = 3.0
        normalized_v = normalized_random_3d_vector()
        self.v = random_vector(normalized_v, self.distance)

    def test(self):
        print("random_3d_vector: ", self.v)
        print("displacement_distance: ", self.distance)
        print("norm: ", np.linalg.norm(self.v))
        self.assertLessEqual(np.linalg.norm(self.v), self.distance)
