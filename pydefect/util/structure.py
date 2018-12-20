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


def get_symmetry_dataset(structure, symprec):
    """
    Args:
        structure (Structure):
            Pymatgen Structure class object
    """
    cell = structure_to_spglib_cell(structure)
    return spglib.get_symmetry_dataset(cell, symprec=symprec)


def get_point_group_from_dataset(sym_dataset, coords, lattice, symprec):
    """
    Args:
        sym_dataset (dict):
            spglib get_symmetry_dataset.
        coords (list):
            Fractional coordinates.
        lattice (numpy.array):
            3x3 numpy array
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
    """
    full_rotations = sym_dataset["rotations"]
    translations = sym_dataset["translations"]
    rotations = get_rotations(coords, lattice, full_rotations, translations,
                              symprec)
    return get_point_group_from_rotations(rotations)


def get_point_group_from_rotations(rotations):
    ptg = spglib.get_pointgroup(rotations)
    return ptg[0].strip(), ptg[2]


def get_rotations(coords, lattice, rotations, translations, symprec):
    """
    Args:
        coords (list):
            Cartesian coordinates.
        lattice (numpy.array):
            3x3 numpy array
        rotations (dict):
            list of 3x3 rotation matrix.
        translations (dict):
            list of 3 translation column.
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
    """
    site_symmetries = []

    for r, t in zip(rotations, translations):
        rot_pos = np.dot(coords, r.T) + t
        diff = coords - rot_pos
        diff -= np.rint(diff)
        diff = np.dot(diff, lattice)
        if np.linalg.norm(diff) < symprec:
            site_symmetries.append(r)

    return np.array(site_symmetries, dtype='intc')


def structure_to_spglib_cell(structure):
    """
    Return a *cell* tuple parsed by spglib that is composed of lattice
    parameters, atomic positions in fractional coordinates, and corresponding
    atomic numbers.
    Args:
        structure (Structure):
            Pymatgen Structure class object
    """
    lattice = list(structure.lattice.matrix)
    positions = structure.frac_coords.tolist()
    atomic_numbers = [i.specie.number for i in structure.sites]
    return lattice, positions, atomic_numbers


def spglib_cell_to_structure(cell):
    """
    Return a pymatgen Structure class object from spglib cell tuple.
    Args:
        cell (3 tuple):
            Lattice parameters, atomic positions in fractional coordinates,
            and corresponding atom numbers
    """
    species = [symbols_to_atom[i] for i in cell[2]]
    return Structure(cell[0], species, cell[1])


def find_spglib_standard_conventional(structure, symprec=1e-02):
    """
    Return a standard conventional unit cell.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
    """
    cell = structure_to_spglib_cell(structure)
    return spglib_cell_to_structure(
        spglib.standardize_cell(cell, to_primitive=False, no_idealize=False,
                                symprec=symprec))


def find_spglib_standard_primitive(structure, symprec=1e-02):
    """
    Return a primitive unit cell.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
    """
    cell = structure_to_spglib_cell(structure)
    return spglib_cell_to_structure(spglib.find_primitive(cell, symprec))


def find_hpkot_primitive(structure, symprec=1e-02, angle_tolerance=-1.0):
    """
    Return a hpkot primitive unit cell.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analysis.
    """
    cell = structure_to_spglib_cell(structure)
    res = seekpath.get_explicit_k_path(structure=cell, symprec=symprec,
                                       angle_tolerance=angle_tolerance)

    return seekpath_to_hpkot_structure(res)


def structure_to_seekpath(structure, time_reversal=True, ref_distance=0.025,
                          recipe='hpkot', threshold=1.e-7, symprec=1e-02,
                          angle_tolerance=-1.0):
    """
    Return the full information for the band path of the given Structure class
    object generated by seekpath.
    Args:
        structure (Structure):
            Pymatgen Structure class object
        time_reversal (bool):
            If the time reversal symmetry exists
        ref_distance (float):
            Distance for the k-point mesh.
        threshold (float):
            To use to verify if we are in edge case (see seekpath).
        symprec (float):
            Distance tolerance in cartesian coordinates Unit is compatible with
            the cell.
        angle_tolerance (float):
            Angle tolerance used for symmetry analysis.
    """
    cell = structure_to_spglib_cell(structure)
    res = seekpath.get_explicit_k_path(cell,
                                       with_time_reversal=time_reversal,
                                       reference_distance=ref_distance,
                                       recipe=recipe,
                                       threshold=threshold,
                                       symprec=symprec,
                                       angle_tolerance=angle_tolerance)

    # If numpy.allclose is too strict in pymatgen.core.lattice __eq__,
    # make almost_equal
    if structure.lattice == seekpath_to_hpkot_structure(res).lattice:
        return res
    else:
        raise NotStandardizedPrimitiveError(
            "The given structure is not standardized primitive cell.")


def seekpath_to_hpkot_structure(res):
    """
    Return a pymatgen Structure class object from seekpath res dictionary.
    Args:
        res (dict):
            seekpath res dictionary.
    """
    lattice = res["primitive_lattice"]
    element_types = res["primitive_types"]
    species = [symbols_to_atom[i] for i in element_types]
    positions = res["primitive_positions"]
    return Structure(lattice, species, positions)


def get_coordination_environment(structure, index, factor=1.2):
    """
    Args:
        structure (Structure):
            pmg Structure class object
        index (int):
            The atomic index
        factor (float):
            Multiple number of the distance of the sum of the averaged ionic
            radii that is used to detect the coordination.

    Return:
        coords (dict):
            values are tuples of Element object and distance.
            example  {'O': [(PeriodicSite: O (0.0000, 0.0000, -2.1234)
            [-0.5000, -0.5000, 0.5000], 2.123447), ...
    """
    element = structure.sites[index].specie
    c = charge[element.symbol]
    ionic_radius = element.ionic_radii[c]
    coords = []

    for e in structure.types_of_specie:
        c2 = charge[e.symbol]

        if c * c2 >= 0:
            continue

        another_ionic_radius = e.ionic_radii[c2]

        cutoff = (ionic_radius + another_ionic_radius) * factor
        neighbors = structure.get_neighbors(structure.sites[index], cutoff)
        for site in neighbors:
            if site[0].specie == e:
                coords.append(site)

    return coords


def get_coordination_distances(structure, index, factor=1.2):
    """
    Return:
        coords (dict):
            example {"Mg": [1.92, 1.73], "Al": [2.01, 2.11]}
    """
    coordination_environment = \
        get_coordination_environment(structure, index, factor)
    coordination_distances = defaultdict(list)
    for ce in coordination_environment:
        coordination_distances[ce[0].species_string].append(ce[1])

    return coordination_distances


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
        structure.translate_sites(site, vector, frac_coords=False)

    return structure, sites


class NotStandardizedPrimitiveError(Exception):
    pass
