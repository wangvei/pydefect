# -*- coding: utf-8 -*-
from itertools import product, combinations, chain
from math import acos, floor, degrees
import numpy as np
from tqdm import tqdm
from typing import Union, List

from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import DummySpecie, Specie
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import pbc_shortest_vectors

from obadb.util.structure_handler import get_symmetry_dataset, get_rotations

import spglib

from pydefect.util.logger import get_logger
from pydefect.util.math_tools import normalized_random_3d_vector, random_vector
from pydefect.core.error_classes import StructureError
from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def perturb_neighboring_atoms(structure: Structure,
                              center: List[float],
                              cutoff: float,
                              distance: float,
                              inserted_atom_indices: List[int]):
    """ Return the structure with randomly perturbed atoms near the center

    Args:
        structure (Structure):
            pmg Structure class object
        center (list):
            Fractional coordinates of a central position.
        cutoff (float):
            Radius of a sphere in which atoms are perturbed.
        distance (float):
            Max displacement_distance for the perturbation.
        inserted_atom_indices (list):
            Inserted atom indices, which will not be perturbed.
    """
    perturbed_structure = structure.copy()
    cartesian_coords = structure.lattice.get_cartesian_coords(center)
    neighbors = structure.get_sites_in_sphere(
        cartesian_coords, cutoff, include_index=True)
    if not neighbors:
        logger.warning("No neighbors withing the cutoff {}.".format(cutoff))

    sites = []
    # Since translate_sites accepts only one vector, we need to iterate this.
    for i in neighbors:
        if i[2] in inserted_atom_indices:
            continue
        vector = random_vector(normalized_random_3d_vector(), distance)
        site_index = i[2]
        sites.append(site_index)
        perturbed_structure.translate_sites(site_index, vector,
                                            frac_coords=False)

    return perturbed_structure, sites


def get_displacements(final_structure: Structure,
                      initial_structure: Structure,
                      center: Union[list, int],
                      anchor_atom_index: int = None):
    """ Return the list of displacements.

    Args:
        final_structure (Structure):
            Relaxed Structure
        initial_structure (Structure):
            Initial Structure
        center (list):
            Fractional coordinates of a central position or an atom index
        anchor_atom_index (int):
            Atom index that is assumed not to be moved during the structure
            optimization, which is usually the farthest atom from a defect.
    """
    if len(final_structure) != len(initial_structure):
        raise StructureError("The number of atoms are different between two "
                             "input structures.")
    elif final_structure.lattice != initial_structure.lattice:
        raise StructureError("The lattice constants are different between two "
                             "input structures.")

    center = np.array(center)

    if anchor_atom_index:
        offset_coords = final_structure[anchor_atom_index].frac_coords - \
                        initial_structure[anchor_atom_index].frac_coords
    else:
        offset_coords = [0.0, 0.0, 0.0]
    offset_coords = np.array(offset_coords)

    disp_vectors = []
    disp_norms = []
    for f, i in zip(final_structure, initial_structure):
        # disp_vector is in cartesian coordinates.
        disp_vector, d2 = pbc_shortest_vectors(initial_structure.lattice,
                                               i.frac_coords,
                                               f.frac_coords - offset_coords,
                                               return_d2=True)
        disp_vectors.append(disp_vector[0][0])
        disp_norms.append(d2[0][0] ** 0.5)

    initial_distances = distance_list(initial_structure, center)
    final_distances = distance_list(final_structure, center + offset_coords)

    angles = []
    for i, d, f in zip(initial_distances, disp_norms, final_distances):
        if i < 1e-5 or d < 1e-5:
            angles.append(None)
            continue
        x = (i * i + d * d - f * f) / (2 * i * d)
        if 1 < x < 1 + 1e-8:
            x = 1
        if -1 - 1e-8 < x < -1:
            x = -1
        angles.append(degrees(acos(x)))

    # angles are nan when the displacements are zero or diverged.
    return initial_distances, final_distances, disp_vectors, disp_norms, angles


def defect_center_from_coords(defect_coords: list,
                              structure: Structure) -> list:
    """Return defect center from given coordinates.

    Args:
        defect_coords (list):
        structure (Structure):
    """
    # First defect_coords is used as a base point under periodic boundary
    # condition. Here, we are aware of the case when two defect positions
    # are, e.g., [0.01, 0.01, 0.01] and [0.99, 0.99, 0.99].
    base = defect_coords[0]
    shortest_defect_coords = []

    for dc in defect_coords:
        diff = structure.lattice.get_distance_and_image(base, dc)[1]
        shortest_defect_coords.append(dc + diff)

    # np.array([[0, 0.1, 0.2], [0.3, 0.4, 0.5]]).transpose() =
    # np.array([[0, 0.3], [0.1, 0.4], [0.2, 0.5]])
    return [np.mean(i) for i in np.array(shortest_defect_coords).transpose()]


def distance_list(structure, coords):
    """ Return a list of the shortest distances between a point and its images

    Args:
       structure (Structure):
           pmg structure class object
       coords (1x3 numpy array):
           Fractional coordinates
    """
    return [structure.lattice.get_distance_and_image(host, coords)[0]
            for host in structure.frac_coords]


def distances_from_defect_center(structure, defect_entry):
    """ Returns a list of distances at atomic sites from defect center

    Note that in the case of an interstitial-type defect, zero is also set for
    the interstitial site itself.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        defect_entry (DefectEntry):
            DefectEntry class object considered
    """
    return distance_list(structure, defect_entry.defect_center)


def fold_positions(structure):
    """ Fold atomic positions in fractional coords into from 0 to 1.

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
    """ Same as fold_positions but for poscar

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
    return np.array(distances)


def are_distances_same(lattice: Lattice,
                       points: list,
                       compared_distances: np.array,
                       symprec: float,
                       ) -> bool:
    """ whether the calculated inter-points distances are the same

    Args:
        lattice (Lattice):
        points (list):
        compared_distances (np.array):
        symprec (float):
    """
    distances = []
    for a, b in combinations(points, 2):
        distance = lattice.get_distance_and_image(a, b)[0]
        if any(abs(distance - compared_distances) < 0.0001):
            distances.append(distance)
        else:
            return False

    return all(abs(np.sort(np.array(distances))
                   - np.sort(compared_distances)) < symprec)


def create_saturated_interstitial_structure(structure: Structure,
                                            inserted_atom_coords: list,
                                            species: list = None,
                                            dist_tol: float = 0.1,
                                            symprec: float = SYMMETRY_TOLERANCE,
                                            angle_tolerance: float = ANGLE_TOL):
    """ generates the sublattice for it based on the structure's space group.

    Args:
        structure (Structure):
        inserted_atom_coords (list):
            List of 3x1 fractional coords or 3x1 fractional coords
        species (list): If not provided, "H" are inserted.
        dist_tol (float):
            changing distance tolerance of saturated structure,
            allowing for possibly overlapping sites
            but ensuring space group is maintained
        symprec (float):
        angle_tolerance (float):

    Returns:
        saturated_defect_structure (Structure):
            Saturated structure decorated with equivalent interstitial sites.
        inserted_atom_indices (list):
            The representative inserted atom indices.
        symmetry_dataset (dict):
            Spglib style symmetry dataset with the following propeties.
                number: International space group number
                international: International symbol
                hall: Hall symbol
                transformation_matrix: Transformation matrix from lattice of
                input cell to Bravais lattice L^bravais = L^original * Tmat
                origin shift: Origin shift in the setting of "Bravais lattice"
                rotations, translations: Rotation matrices and translation
                vectors. Space group operations are obtained by
                [(r,t) for r, t in zip(rotations, translations)]
                wyckoffs: Wyckoff letters
    """
    if species and len(inserted_atom_coords) != len(species):
        raise ValueError("Numbers of inserted_atom_coords and species differ.")
    if isinstance(inserted_atom_coords[0], float):
        inserted_atom_coords = [inserted_atom_coords[:]]

    sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
    symmops = sga.get_symmetry_operations()
    symmetry_dataset = sga.get_symmetry_dataset()
    saturated_defect_structure = structure.copy()
    saturated_defect_structure.DISTANCE_TOLERANCE = dist_tol

    inserted_atom_indices = []
    for i, coord in enumerate(inserted_atom_coords):
        # Whether the inserted atoms already exist.
        try:
            inserted_atom_index = \
                neighboring_atom_indices(
                    saturated_defect_structure, coord, dist_tol)[0]
        except IndexError:
            inserted_atom_index = len(saturated_defect_structure)

            for symmop in symmops:
                specie = species[i] if species else DummySpecie()
                new_interstitial_coords = symmop.operate(coord[:])
                try:
                    saturated_defect_structure.append(specie,
                                                      new_interstitial_coords,
                                                      validate_proximity=True)
                except ValueError:
                    pass

        inserted_atom_indices.append(inserted_atom_index)

    return saturated_defect_structure, inserted_atom_indices, symmetry_dataset


def neighboring_atom_indices(structure, coord, dist_tol):

    neighboring_indices = []
    for j, site in enumerate(structure):
        distance = structure.lattice.get_distance_and_image(coord,
                                                            site.frac_coords)[0]
        if distance < dist_tol:
            neighboring_indices.append(j)

    return neighboring_indices


def count_equivalent_clusters(perfect_structure: Structure,
                              inserted_atom_coords: list,
                              removed_atom_indices: list,
                              displacement_distance: float,
                              symprec: float = SYMMETRY_TOLERANCE,
                              angle_tolerance: float = ANGLE_TOL):
    """Return a dict of atomic distances

    NOTE: This can be calculated using mathematics more elegantly.
          Symmetry of complex defect is direct product of each defect symmetry?

    Args:
        perfect_structure (Structure):
            Supercell is assumed to big enough.
        inserted_atom_coords (list):
        removed_atom_indices (list):
            Needs to begin from 0.
        displacement_distance (float)
        symprec (float):
        angle_tolerance (float):

    Returns:
        count (int):

    """
    inserted_atom_coords = list(inserted_atom_coords)
    lattice = perfect_structure.lattice

    # If the removed atoms are substituted, the removed sites shouldn't be
    # considered as vacant sites for counting multiplicity.
    vacant_atom_indices = []
    for removed_atom_index in removed_atom_indices:
        removed_atom_coord = perfect_structure[removed_atom_index].frac_coords
        for inserted_atom_coord in inserted_atom_coords:

            distance = lattice.get_distance_and_image(
                inserted_atom_coord, removed_atom_coord)[0]

            if distance < displacement_distance:
                break
        else:
            vacant_atom_indices.append(removed_atom_index)

    if inserted_atom_coords:
        saturated_structure, inserted_atom_indices, _ = \
            create_saturated_interstitial_structure(
                structure=perfect_structure.copy(),
                inserted_atom_coords=inserted_atom_coords)
    else:
        saturated_structure = perfect_structure
        inserted_atom_indices = []

    sga = SpacegroupAnalyzer(saturated_structure, symprec, angle_tolerance)
    symmetrized_saturated_structure = sga.get_symmetrized_structure()
    # symmetrically equivalent sites are grouped.
    equiv_indices = symmetrized_saturated_structure.equivalent_indices

    defect_atom_indices = vacant_atom_indices + inserted_atom_indices

    # count the symmetrically equivalent defects
    equiv_num_defects = []
    for i in equiv_indices:
        atoms_from_i = [j for j in defect_atom_indices if j in i]
        equiv_num_defects.append(len(atoms_from_i))

    # Here, interatomic distances of original cluster are calculated.
    orig_atom_frac_coords = \
        [saturated_structure[i].frac_coords for i in defect_atom_indices]
    orig_distances = atomic_distances(lattice, orig_atom_frac_coords)

    # Calculate the combination of each equivalent defect site.
    groups = []
    for i, num_atoms in enumerate(equiv_num_defects):
        groups.append(list(combinations(equiv_indices[i], num_atoms)))

    # Consider all the possible defect site combinations.
    count = 0
    for i in tqdm(list(product(*groups))):
        frac_coords = \
            [saturated_structure[i].frac_coords for i in list(chain(*i))]

        if are_distances_same(lattice, frac_coords, orig_distances, symprec):
            count += 1

    return count


# def count_equivalent_clusters2(perfect_structure: Structure,
#                                inserted_atom_coords: list,
#                                removed_atom_indices: list,
#                                displacement_distance: float,
#                                symprec: float = SYMMETRY_TOLERANCE,
#                                angle_tolerance: float = ANGLE_TOL):
#     """Return a dict of atomic distances

# NOTE: This can be calculated using mathematics more elegantly.
#       Symmetry of complex defect is direct product of each defect symmetry?

# Args:
#     perfect_structure (Structure):
#         Supercell is assumed to big enough.
#     inserted_atom_coords (list):
#     removed_atom_indices (list):
#         Needs to begin from 0.
#     displacement_distance (float)
#     symprec (float):
#     angle_tolerance (float):

# Returns:
#     count (int):

# """
# inserted_atom_coords = list(inserted_atom_coords)

# # If the removed atoms are substituted, the removed sites shouldn't be
# # considered as vacant sites for counting multiplicity.
# vacant_atom_indices = []
# for removed_atom_index in removed_atom_indices:
#     for inserted_atom_coord in inserted_atom_coords:

# distance_and_image = \
#     perfect_structure.lattice.get_distance_and_image(
#         inserted_atom_coord,
#         perfect_structure[removed_atom_index].frac_coords)

#     distance = distance_and_image[0]
#     if distance < displacement_distance:
#         break
# else:
#     vacant_atom_indices.append(removed_atom_index)

# saturated_structure, inserted_atom_indices = \
#     create_saturated_interstitial_structure(
#         structure=perfect_structure.copy(),
#         inserted_atom_coords=inserted_atom_coords)

# sga = SpacegroupAnalyzer(saturated_structure, symprec, angle_tolerance)
# symmetrized_saturated_structure = sga.get_symmetrized_structure()
# # symmetrically equivalent sites are grouped.
# equiv_indices = symmetrized_saturated_structure.equivalent_indices
# print(equiv_indices)

# sym_dataset = get_symmetry_dataset(symmetrized_saturated_structure, symprec,
#                                    angle_tolerance)

# lattice = symmetrized_saturated_structure.lattice
# defect_atom_indices = vacant_atom_indices + inserted_atom_indices


# # count the symmetrically equivalent defects
# count = 0
# for sublattice_equiv_indices in equiv_indices:
#     sublattice_defect_indices = \
#         [j for j in defect_atom_indices if j in sublattice_equiv_indices]
#     if not sublattice_defect_indices:
#         continue
#     else:
#         sublattice_defect_coords = \
#             np.array([symmetrized_saturated_structure[i].frac_coords
#                       for i in sublattice_defect_indices])

# # site symmetry
# symmetry_operation_number = \
#     get_point_group_op_number(sym_dataset, sublattice_defect_coords[0],
#                               lattice.matrix, symprec)

# structure = Structure(lattice=lattice,
#                       species=["H"] * len(sublattice_defect_coords),
#                       coords=sublattice_defect_coords)

# equivalent_symmetry_reduction = \
#     len(SpacegroupAnalyzer(structure, symprec, angle_tolerance).
#         get_space_group_operations())

# print("aaaaaaaaaaaaaaaa")
# print(sublattice_defect_coords)
# print(symmetry_operation_number)
# print(equivalent_symmetry_reduction)


# count += len(sublattice_equiv_indices) * symmetry_operation_number / equivalent_symmetry_reduction

# return count


def get_point_group_op_number(sym_dataset: dict,
                              coords: list,
                              lattice: np.array,
                              symprec: float = SYMMETRY_TOLERANCE):
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
    return len(rotations)


def get_symmetry_multiplicity(sym_dataset: dict,
                              coords: list,
                              lattice: np.array,
                              symprec: float = SYMMETRY_TOLERANCE):

    return int(len(sym_dataset["rotations"]) /
               get_point_group_op_number(sym_dataset, coords, lattice, symprec))


def first_appearance_index(structure: Structure,
                           specie: Union[str, Specie]) -> int:
    """Return first index where the specie appears. Return 0 if it doesn't exist

    Used for inserting an element to a Structure, so if the element does not
    exist, 0 is returned such that the new specie is inserted to
    the first place by structure.insert(inserted_index, specie, coords).

    :param structure:
    :param specie:
    :return:
    """
    if specie in structure.symbol_set:
        return min(structure.indices_from_symbol(specie))
    else:
        return 0


class ModSpacegroupAnalyzer(SpacegroupAnalyzer):
    def get_refined_structure(self, is_sorted=False):
        """
        Get the refined structure based on the symmetry. The refined
        structure is a *conventional* cell setting with atoms moved to the
        expected symmetry positions.

        Returns:
            Refined structure.
        """
        # Atomic positions have to be specified by scaled positions for spglib.
        lattice, scaled_positions, numbers \
            = spglib.refine_cell(self._cell, self._symprec, self._angle_tol)

        species = [self._unique_species[i - 1] for i in numbers]
        s = Structure(lattice, species, scaled_positions)

        return s.get_sorted_structure() if is_sorted else s

    def get_anchored_refined_structure(self, index):
        """
        Get the refined structure based on detected symmetry. The refined
        structure is a *conventional* cell setting with atoms moved to the
        expected symmetry positions.

        Returns:
            Refined structure.
        """
        # Atomic positions have to be specified by scaled positions for spglib.
        lattice, scaled_positions, numbers \
            = spglib.refine_cell(self._cell, self._symprec, self._angle_tol)

        species = [self._unique_species[i - 1] for i in numbers]
        s = Structure(lattice, species, scaled_positions)
        return s.get_sorted_structure(), numbers


