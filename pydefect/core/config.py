# -*- coding: utf-8 -*-
"""Defaults that determine condition of automatic defect calculations. """

# electronegativity difference for antisites and substitutional impurities
ELECTRONEGATIVITY_DIFFERENCE = 1.0
# Maximum displacement displacement_distance
DISPLACEMENT_DISTANCE = 0.2
SYMMETRY_TOLERANCE = 0.01
# The following must be used after structure optimization anytime.
DEFECT_SYMMETRY_TOLERANCE = 0.07
# For seeking inequivalent interstitial sites.
INTERSTITIAL_SYMPREC = 0.15

ANGLE_TOL = 5

DEFECT_KPT_DENSITY = 2

# Factor multiplied with the minimum distance, of which distance determine
# the coordination environment and in which atoms are perturbed. 1.3 is less
# than sqrt(2), meaning the 2NNs are excluded in e.g., the rocksalt structure.
CUTOFF_FACTOR = 1.3

COLOR = [
    "xkcd:blue",
    "xkcd:brown",
    "xkcd:crimson",
    "xkcd:darkgreen",
    "xkcd:gold",
    "xkcd:magenta",
    "xkcd:orange",
    "xkcd:darkblue",
    "xkcd:navy",
    "xkcd:red",
    "xkcd:olive",
    "xkcd:black",
    "xkcd:indigo"
]


