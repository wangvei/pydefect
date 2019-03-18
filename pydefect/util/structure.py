# -*- coding: utf-8 -*-

from collections import defaultdict

import numpy as np
import seekpath
import spglib

from pymatgen import Structure
from pymatgen.core.periodic_table import Element

from pydefect.database.atom import symbols_to_atom, charge
from pydefect.util.math import normalized_random_3d_vector, random_vector

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "April 4, 2018"


def perturb_neighbors(structure, center, cutoff, distance):
    """
    Return the structure with randomly perturbed atoms near the input point in
    structure.

    Args:
        structure (Structure):
            pmg Structure class object
        center (list):
            Fractional coordinates of a central position.
        cutoff (float):
            Radius of a sphere in which atoms are perturbed.
        distance (float):
            Max distance for the perturbation.
    """
    if len(center) == 3:
        from copy import deepcopy
        perturbed_structure = deepcopy(structure)
        cartesian_coords = structure.lattice.get_cartesian_coords(center)
        neighbors = structure.get_sites_in_sphere(
            cartesian_coords, cutoff, include_index=True)
    else:
        raise ValueError

    sites = []
    # Since translate_sites accepts only one vector, we need to iterate this.
    for i in neighbors:
        vector = random_vector(normalized_random_3d_vector(), distance)
        site = i[2]
        sites.append(site)
        perturbed_structure.translate_sites(site, vector, frac_coords=False)

    return perturbed_structure, sites


