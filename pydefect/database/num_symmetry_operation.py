# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def num_symmetry_operation(point_group):
    """ Return number of symmetry operations from Hermann–Mauguin notation. """
    d = {"1":     1,
         "-1":    2,
         "2":     2,
         "m":     2,
         "2/m":   4,
         "222":   4,
         "2mm":   4,
         "mm2":   4,
         "mmm":   8,
         "2/mmm": 8,
         "4":     4,
         "-4":    4,
         "4/m":   8,
         "422":   8,
         "4mm":   8,
         "-4m2":  8,
         "-42m":  8,
         "4/mmm": 16,
         "3":     3,
         "-3":    6,
         "32":    6,
         "3m":    6,
         "-3m":   12,
         "6":     6,
         "-6":    6,
         "6/m":   12,
         "622":   12,
         "6mm":   12,
         "-6m2":  12,
         "6/mmm": 24,
         "23":    12,
         "m3":    24,
         "432":   24,
         "-43m":  24,
         "m-3m":  48}

    point_group = "".join([s for s in point_group if s != "."])

    return d[point_group]
