# -*- coding: utf-8 -*-
import numpy as np
import os
import unittest

from pydefect.util.math_tools import normalized_random_3d_vector, random_vector

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class NormalizedRandom3dVectorTest(unittest.TestCase):

    def setUp(self):
        self.v = normalized_random_3d_vector()

    def test(self):
        print("normalized_random_3d_vector: ", self.v)
        print("norm: ", np.linalg.norm(self.v))
        self.assertAlmostEqual(np.linalg.norm(self.v), 1.0)


class RandomVectorTest(unittest.TestCase):

    def setUp(self):
        self.distance = 3.0
        normalized_v = normalized_random_3d_vector()
        self.v = random_vector(normalized_v, self.distance)

    def test(self):
        print("random_3d_vector: ", self.v)
        print("displacement_distance: ", self.distance)
        print("norm: ", np.linalg.norm(self.v))
        self.assertLessEqual(np.linalg.norm(self.v), self.distance)