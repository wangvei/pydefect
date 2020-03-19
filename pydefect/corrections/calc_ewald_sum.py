from functools import reduce
from typing import Tuple

import numpy as np
from numpy import dot, pi, exp, cos
from scipy.special import erfc


# TODO: use cython
def calc_ewald_sum(dielectric_tensor: np.ndarray,
                   real_lattice_set: np.ndarray,
                   reciprocal_lattice_set: np.ndarray,
                   mod_ewald_param: float,
                   root_det_epsilon: float,
                   volume: float,
                   ) -> Tuple[float, float]:
    """Return real and reciprocal Ewald summations at given parameters"""

    epsilon_inv = np.linalg.inv(dielectric_tensor)
    real_sum = 0
    # Skip the potential caused by the defect itself
    for v in real_lattice_set:
        root_r_inv_epsilon_r = np.sqrt(reduce(dot, [v.T, epsilon_inv, v]))
        real_sum += \
            erfc(mod_ewald_param * root_r_inv_epsilon_r) / root_r_inv_epsilon_r
    real_part = real_sum / (4 * pi * root_det_epsilon)

    # Ewald reciprocal part
    # sum exp(-g * epsilon * g / (4 * ewald ** 2)) / g * epsilon * g [1/A]
    reciprocal_sum = 0
    for g in reciprocal_lattice_set:
        g_epsilon_g = reduce(dot, [g.T, dielectric_tensor, g])
        reciprocal_sum += \
            (exp(- g_epsilon_g / 4.0 / mod_ewald_param ** 2)
             / g_epsilon_g * cos(dot(g, np.zeros(3))))  # [A^2]
    reciprocal_part = reciprocal_sum / volume

    return real_part, reciprocal_part
