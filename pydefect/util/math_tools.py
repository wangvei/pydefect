# -*- coding: utf-8 -*-
import numpy as np

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def normalized_random_3d_vector():
    """
    Generates a random 3d unit vector with a uniform spherical distribution.
    stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """
    phi = np.random.uniform(0, np.pi * 2)
    cos_theta = np.random.uniform(-1, 1)
    theta = np.arccos(cos_theta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def random_vector(normed_vector: np.array,
                  distance: float):
    """
    Returns a vector scaled by displacement_distance * x, where 0.1 < x < 1.
    The finite minimum displacement is needed for detecting perturbed sites.

    Args:
        normed_vector (3x1 array): Normed 3d vector.
        distance (float): displacement_distance
    """
    return normed_vector * distance * (np.random.random() * 0.9 + 0.1)
