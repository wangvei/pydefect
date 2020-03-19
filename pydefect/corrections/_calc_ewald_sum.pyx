# -*- coding: utf-8 -*-

from functools import reduce
from typing import Tuple
from cython cimport boundscheck, wraparound


import numpy as np
from numpy import dot, pi, exp, cos
from scipy.special import erfc

@boundscheck(False)
@wraparound(False)
def calc_ewald_sum(float[:, :] dielectric_tensor,
                   double[:, :] real_lattice_set,
                   double[:, :] reciprocal_lattice_set,
                   float mod_ewald_param,
                   float root_det_epsilon,
                   float volume):
    """Return real and reciprocal Ewald summations at given parameters"""

    cdef float[:, :] epsilon_inv
    cdef float root_r_inv_epsilon_r, g_epsilon_g
    cdef float real_sum = 0.0, real_part, reciprocal_sum = 0.0, reciprocal_part

    epsilon_inv = np.linalg.inv(dielectric_tensor)
    # Skip the potential caused by the defect itself
    for v in real_lattice_set:
        root_r_inv_epsilon_r = np.sqrt(reduce(dot, [v.T, epsilon_inv, v]))
        real_sum += \
            erfc(mod_ewald_param * root_r_inv_epsilon_r) / root_r_inv_epsilon_r
    real_part = real_sum / (4 * pi * root_det_epsilon)

    # Ewald reciprocal part
    # sum exp(-g * epsilon * g / (4 * ewald ** 2)) / g * epsilon * g [1/A]
    for g in reciprocal_lattice_set:
        g_epsilon_g = reduce(dot, [g.T, dielectric_tensor, g])
        reciprocal_sum += \
            (exp(- g_epsilon_g / 4.0 / mod_ewald_param ** 2)
             / g_epsilon_g * cos(dot(g, np.zeros(3))))  # [A^2]
    reciprocal_part = reciprocal_sum / volume

    return real_part, reciprocal_part