# -*- coding: utf-8 -*-
from itertools import product, combinations, chain
from math import acos, pi, floor
import numpy as np
from tqdm import tqdm
from typing import Union

from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.util.math import normalized_random_3d_vector, random_vector
from pydefect.core.error_classes import ImproperInputStructureError
from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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
            Max displacement_distance for the perturbation.
    """
    if len(center) == 3:
        perturbed_structure = structure.copy()
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
    Return the shortest displacement_distance between two points in fractional coordinates
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


def fold_positions(structure):
    """
    Fold atomic positions with fractional coords less than 0 or larger than 1
    into from 0 to 1.
    For example, coords of site changes from [-0.3, 1.9, 0.5] to [0.7, 0.9, 0.5]

    Args:
        structure(Structure):

    Returns:
        Structure
    """
    for i, site in enumerate(structure):
        modification_vector = [-floor(v) for v in site.frac_coords]
        structure.translate_sites(i, modification_vector)
    return structure


def fold_positions_in_poscar(poscar):
    """

    Args:
        poscar (Poscar):

    Returns:
        Poscar:

    """
    s = poscar.structure
    fold_positions(s)
    return Poscar(s)


def atomic_distances(lattice: Lattice,
                     points: list
                     ) -> np.array:
    """ return a list of distances between the given points.
    Args:
        lattice (Lattice):
        points (list):
    """
    distances = []
    for a, b in combinations(points, 2):
        distances.append(lattice.get_distance_and_image(a, b)[0])
    return np.sort(np.array(distances))


def check_distances(lattice: Lattice,
                    points: list,
                    orig_distances: np.array,
                    symprec: float,
                    ) -> bool:
    """ return a list of distances between the given points.
    Args:
        lattice (Lattice):
        points (list):
        orig_distances (np.array):
        symprec (float):
    """
    distances = []
    for a, b in combinations(points, 2):
        distance = lattice.get_distance_and_image(a, b)[0]
        if any(abs(distance - orig_distances) < 0.0001):
            distances.append(distance)
        else:
            return False

    if all(abs(np.sort(np.array(distances)) - orig_distances) < symprec):
        return True
    else:
        return False


def create_saturated_interstitial_structure(structure: Structure,
                                            inserted_atom_coords: list,
                                            dist_tol: float = 0.1):
    """ generates the sublattice for it based on the structure's space group.
    Originally copied from pymatgen.analysis.defects.core

    Args:
        structure (Structure):
        inserted_atom_coords (list):
        dist_tol (float):
            changing distance tolerance of saturated structure,
            allowing for possibly overlapping sites
            but ensuring space group is maintained

    Returns:
        Structure object decorated with interstitial site equivalents and
        list of inserted atom indices.
    """
    sga = SpacegroupAnalyzer(structure)
    sg_ops = sga.get_symmetry_operations()

    saturated_defect_structure = structure.copy()
    saturated_defect_structure.DISTANCE_TOLERANCE = dist_tol

    lattice = saturated_defect_structure.lattice

    from pymatgen.core.periodic_table import DummySpecie
    from pymatgen.core.sites import PeriodicSite

    inserted_atom_indices = []
    for i in inserted_atom_coords:
        already_exist = False
        for j, site in enumerate(saturated_defect_structure):
            if lattice.get_distance_and_image(i, site.frac_coords)[0] < dist_tol:
                inserted_atom_indices.append(j)
                already_exist = True

        if already_exist:
            continue

        for sgo in sg_ops:
            if already_exist is False:
                inserted_atom_indices.append(len(saturated_defect_structure))

            new_interstitial_coords = sgo.operate(i[:])
            poss_new_site = PeriodicSite(
                    DummySpecie(),
                    new_interstitial_coords,
                    saturated_defect_structure.lattice)

            already_exist = True

            try:
                # will raise value error if site already exists in structure
                saturated_defect_structure.append(
                            poss_new_site.specie, poss_new_site.frac_coords,
                            validate_proximity=True)
            except ValueError:
                pass

    # do final space group analysis to make sure symmetry not lowered by
    # saturating defect structure
    # saturated_sga = SpacegroupAnalyzer(saturated_defect_struct)
    # print(saturated_defect_struct)

    # if saturated_sga.get_space_group_number() != sga.get_space_group_number():
    #     raise ValueError("Warning! Interstitial sublattice generation "
    #                      "has changed space group symmetry. I recommend "
    #                      "reducing dist_tol and trying again...")

    return saturated_defect_structure, inserted_atom_indices


def count_equivalent_clusters(perfect_structure: Structure,
                              inserted_atom_coords: list,
                              removed_atom_indices: list,
                              symprec: float = SYMMETRY_TOLERANCE,
                              angle_tolerance: float = ANGLE_TOL):
    """Return a dict of atomic distances

    Args:
        perfect_structure (Structure): Supercell is assumed to big enough.
        inserted_atom_coords (list):
        removed_atom_indices (list):
            Needs to begin from 0.
        symprec (float):
        angle_tolerance (float):

    Returns:
        count (int):

    """
    saturated_structure, inserted_atom_indices = \
        create_saturated_interstitial_structure(
            structure=perfect_structure.copy(),
            inserted_atom_coords=inserted_atom_coords)

    # make host site groups for removed_atoms
    sga = SpacegroupAnalyzer(saturated_structure, symprec, angle_tolerance)
    symmetrized_saturated_structure = sga.get_symmetrized_structure()
    equiv_indices = symmetrized_saturated_structure.equivalent_indices
    equiv_num_atoms = []

    for i in equiv_indices:
        atoms_from_i = \
            [j for j in removed_atom_indices + inserted_atom_indices if j in i]
        equiv_num_atoms.append(len(atoms_from_i))

    orig_atom_indices = removed_atom_indices + inserted_atom_indices
    orig_atom_frac_coords = \
        [saturated_structure[i].frac_coords for i in orig_atom_indices]
    orig_distances = atomic_distances(perfect_structure.lattice,
                                      orig_atom_frac_coords)

    groups = []
    for i, num_atoms in enumerate(equiv_num_atoms):
        groups.append(list(combinations(equiv_indices[i], num_atoms)))

    count = 0
    for i in tqdm(list(product(*groups))):
        frac_coords = \
            [saturated_structure[i].frac_coords for i in list(chain(*i))]

        if check_distances(perfect_structure.lattice, frac_coords,
                           orig_distances, symprec):
            count += 1

    return count

