# -*- coding: utf-8 -*-
__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


# Following defaults determine the condition of automatic defect calculations.
# electronegativity difference for antisites and substitutional impurities
ELECTRONEGATIVITY_DIFFERENCE = 1.0
# Maximum displacement displacement_distance
DISPLACEMENT_DISTANCE = 0.2
SYMMETRY_TOLERANCE = 0.01
# The following must be used after structure optimization anytime.
DEFECT_SYMMETRY_TOLERANCE = 0.03

ANGLE_TOL = 5

PERFECT_KPT_DENSITY = 3
DEFECT_KPT_DENSITY = 2
ENCUT_FACTOR_STR_OPT = 1.3

# Factor multiplied with the minimum distance, of which distance determine
# the coordination environment and in which atoms are perturbed. 1.4 is less
# than sqrt(2), meaning the 2NNs are excluded in e.g., the rocksalt structure.
CUTOFF_FACTOR = 1.4

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


