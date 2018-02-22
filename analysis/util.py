#!/bin/env python3                                                                                                                       

import numpy as np
from copy import deepcopy
from itertools import product


def calc_min_distance_and_its_v2coord(v1, v2, axis):
    candidate_list = []
    v = np.dot(axis, v2 - v1)
    for index in product((-1, 0, 1), repeat=3):
        index = np.array(index)
        delta_v = np.dot(axis, index)
        distance = np.linalg.norm(delta_v + v)
        candidate_list.append((distance, v2 + index))
    return min(candidate_list, key = lambda t: t[0])


def make_local_structure_for_visualization(structure, defect_pos_frac,
                                           distances_from_defect, rate):
    threshold = min(distances_from_defect) * rate
    structure_visualize = deepcopy(structure) # poscar for visualization
    remove_list = []
    for i, d in enumerate(distances_from_defect):
        if d >= threshold:
            remove_list.append(i)
        else:
            translate = -np.array(defect_pos_frac) + np.array([0.5, 0.5, 0.5])
            structure_visualize.translate_sites(i, translate)
    structure_visualize.remove_sites(remove_list)
    return structure_visualize


