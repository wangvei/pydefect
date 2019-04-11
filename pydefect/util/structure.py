# -*- coding: utf-8 -*-

from collections import defaultdict
from itertools import product

from math import acos, pi
import numpy as np
import seekpath
import spglib
from typing import Union

from pymatgen import Structure
from pymatgen.core.periodic_table import Element

from pydefect.database.atom import symbols_to_atom, charge
from pydefect.util.math import normalized_random_3d_vector, random_vector
from pydefect.util.error_classes import ImproperInputStructureError

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "April 4, 2018"


def perturb_neighboring_atoms(structure, center, cutoff, distance):
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


def get_displacements(final_structure: Structure,
                      initial_structure: Structure,
                      center: Union[list, tuple],
                      anchor_atom_index: int):
    """
    Return the list of displacements.

    Args:
        final_structure (Structure):
            Relaxed Structure
        initial_structure (Structure):
            Initial Structure
        center (list):
            Fractional coordinates of a central position.
    """

    if len(final_structure) != len(initial_structure):
        raise ImproperInputStructureError("The number of atoms are not the"
                                          "same between two input structures.")
    elif final_structure.lattice != initial_structure.lattice:
        raise ImproperInputStructureError("The Lattice constants are different"
                                          "between two input structures.")

    center = list(center)

    offset_coords = final_structure[anchor_atom_index].frac_coords - \
                    initial_structure[anchor_atom_index].frac_coords

    displacements = []
    for f, i in zip(final_structure, initial_structure):
        disp_coords = (f.frac_coords - offset_coords) - i.frac_coords
        displacements.\
            append(initial_structure.lattice.get_cartesian_coords(disp_coords))

    initial_distances = distance_list(initial_structure, center)
    final_distances = distance_list(final_structure, center + offset_coords)

    angles = []
    for i, d, f in zip(initial_distances, displacements, final_distances):
        angles = acos((i * i + d * d - f * f) / (2 * i * d)) / pi * 180

    return initial_distances, final_distances, displacements, angles


def defect_center_from_coords(inserted_atom_coords, removed_atom_coords,
                              structure):
    """
    """
    defect_coords = inserted_atom_coords + removed_atom_coords
    # First defect_coords is used as a base point under periodic boundary
    # condition. Here, we are aware of the case when two defect positions
    # are, e.g., [0.01, 0.01, 0.01] and [0.99, 0.99, 0.99].
    base = defect_coords[0]
    shortest_defect_coords = []

    for dc in defect_coords:
        diff = min_distance_under_pbc(
            np.array(base), np.array(dc), structure.lattice.matrix)[1]
        shortest_defect_coords.append(dc + diff)

    # np.array([[0, 0.1, 0.2], [0.3, 0.4, 0.5]]).transpose() =
    # np.array([[0, 0.3], [0.1, 0.4], [0.2, 0.5]])
    return [np.mean(i) for i in np.array(shortest_defect_coords).transpose()]


def defect_center(defect_entry, structure=None):
    """
    Return a fractional coordinates of the defect center that is calculated
    by averaging the coordinates of vacancies and interstitials.
    When len(defect_coords) == 1, simply returns defect_coords[0].
    First defect_coords is used as a base point when the periodic boundary
    condition is considered.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        defect_entry (DefectEntry):
            Related DefectEntry class object
    """
    # If structure is not given, initial_structure of defect_entry is used.
    if structure is None:
        structure = defect_entry.initial_structure

    inserted_atom_coords = list([structure.frac_coords[k]
                                 for k in defect_entry.inserted_atoms])
    removed_atom_coords = list(defect_entry.removed_atoms.values())

    return defect_center_from_coords(inserted_atom_coords, removed_atom_coords,
                                     structure)


def min_distance_under_pbc(frac1, frac2, lattice_parameters):
    """
    Return the shortest distance between two points in fractional coordinates
    under periodic boundary condition.

    Args:
       frac1 (1x3 np.array):
           1st fractional coordinates
       frac2 (1x3 np.array):
           2nd fractional coordinates
       lattice_parameters (3x3 numpy array):
           a, b, c lattice vectors
    """
    candidate = []
    diff = np.dot(lattice_parameters, frac2 - frac1)

    # (-1, -1, -1), (-1, -1, 0), ..., (1, 1, 1)
    for index in product((-1, 0, 1), repeat=3):
        index = np.array(index)
        delta_diff = np.dot(lattice_parameters, index)
        distance = np.linalg.norm(delta_diff + diff)
        candidate.append([distance, index])

    return min(candidate, key=lambda x: x[0])


def distance_list(structure, coords):
    """
    Return a list of the shortest distances between a point and its images
    under periodic boundary condition.
    Args:
       structure (Structure):
           pmg structure class object
       coords (1x3 numpy array):
           Fractional coordinates
    """
    lattice_parameters = structure.lattice.matrix

    return [min_distance_under_pbc(host, coords, lattice_parameters)[0]
            for host in structure.frac_coords]


def distances_from_point(structure, defect_entry):
    """
    Returns a list of distances at atomic sites from a defect center defined
    by defect_entry. Note that in the case of an interstitial-type defect,
    zero is also set to the interstitial site.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        defect_entry (DefectEntry):
            DefectEntry class object considered
    """
    return distance_list(structure, defect_center(defect_entry, structure))