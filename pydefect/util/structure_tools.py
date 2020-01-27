# -*- coding: utf-8 -*-
from collections import defaultdict
from typing import Union, List, Tuple, Optional

import numpy as np
from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.error_classes import StructureError
from pydefect.database.symmetry import num_symmetry_operation
from pydefect.util.logger import get_logger
from pydefect.util.math import normalized_random_3d_vector
from pymatgen.core.periodic_table import DummySpecie, Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import pbc_shortest_vectors

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def perturb_neighboring_atoms(structure: Structure,
                              center: List[float],
                              cutoff: float,
                              distance: float,
                              inserted_atom_indices: List[int]
                              ) -> Tuple[Structure, list]:
    """ Return the structure with randomly perturbed atoms near the center

    Args:
        structure (Structure):
            pmg Structure class object
        center (list):
            Fractional coordinates of a central position.
        cutoff (float):
            Radius of a sphere in which atoms are perturbed.
        distance (float):
            Max displacement distance for the perturbation.
        inserted_atom_indices (list):
            Inserted atom indices, which will not be perturbed.

    Return:
        Tuple of perturbed Structure and list of perturbed atom indices.
    """
    perturbed_structure = structure.copy()
    cartesian_coords = structure.lattice.get_cartesian_coords(center)
    # neighbor is (PeriodicSite, distance, index)
    neighbors = structure.get_sites_in_sphere(
        pt=cartesian_coords, r=cutoff, include_index=True)
    if not neighbors:
        logger.warning(f"No neighbors withing the cutoff {cutoff}.")

    sites = []
    for i in neighbors:
        # not perturb inserted atoms
        if i[2] in inserted_atom_indices:
            continue
        # Since translate_sites accepts only one vector, iterate this.
        vector = normalized_random_3d_vector() * distance * np.random.random()
        site_index = i[2]
        sites.append(site_index)
        perturbed_structure.translate_sites(
            indices=site_index, vector=vector, frac_coords=False)

    return perturbed_structure, sites


def get_min_distance(structure: Structure) -> float:
    """Return the minimum distance in the inter-atomic distances

    Calculate the minimum distance from the pair distances between an atom and
    another atom in the structure, and it does not include the distance between
    an atom and its image.

    Args:
        structure (Structure): Input structure:

    Return:
        distance
    """
    distance_matrix = structure.distance_matrix
    mask = np.ones(distance_matrix.shape, dtype=bool)
    np.fill_diagonal(mask, 0)

    return distance_matrix[mask].min()


def get_displacements(final_structure: Structure,
                      initial_structure: Structure,
                      defect_center: Union[list, int],
                      anchor_atom_index: int = None) -> dict:
    """ Return information related to atomic displacements.

    Args:
        final_structure (Structure):
            Relaxed final Structure
        initial_structure (Structure):
            Initial Structure
        defect_center (list / int):
            Fractional coordinates of a central position or an atom index.
        anchor_atom_index (int):
            Atom index that is assumed not to be moved during the structure
            optimization, which is usually the farthest atom from a defect.

    Variables:
        drift_frac_coords (np.array):
            A vector showing how the farthest atom drifts during the structure
            optimization.

    Return (dict):
        Keys:
        + initial_distances:
            Distances from a defect in the initial structure.
        + final_distances:
            Distances from a defect in the final structure.
        + displacement_vectors:
            Displacement vectors of atoms from initial to final structures.
        + displacement_norms:
            Norms of displacement vectors.
        + initial_vectors:
            Vectors from a defect position to atoms in the initial supercell.
        + final_vectors:
            Vectors from a defect position to atoms in the final supercell.
        + defect_migration_distance:
            Distance the defect migrates defined only for interstitials,
            antisites, and substituted defects.
    """
    if len(final_structure) != len(initial_structure):
        raise StructureError("The number of atoms are different between two "
                             "input structures.")
    elif final_structure.lattice != initial_structure.lattice:
        logger.warning("The lattice constants are different between two input "
                       "structures. Anchoring the farthest atom is switched "
                       "off as it bears erroneous result.")
        anchor_atom_index = None

    if anchor_atom_index:
        drift_frac_coords = final_structure[anchor_atom_index].frac_coords - \
                             initial_structure[anchor_atom_index].frac_coords
    else:
        drift_frac_coords = np.zeros(3)

    # Except for interstitials, it is assumed that the defect center is not
    # moved w.r.t. the anchoring atom.
    if isinstance(defect_center, int):
        initial_center = initial_structure[defect_center].frac_coords
        final_center = final_structure[defect_center].frac_coords
        if anchor_atom_index:
            final_fcoords = final_structure[defect_center].frac_coords
            initial_fcoords = (initial_structure[defect_center].frac_coords +
                               drift_frac_coords)
            # Note: Defect migration distance is estimated based on the initial
            #       lattice constants. Now, it is fine as anchor_atom_index is
            #       set only when the lattice constants are unchanged.
            _, d2 = pbc_shortest_vectors(lattice=initial_structure.lattice,
                                         fcoords1=final_fcoords,
                                         fcoords2=initial_fcoords,
                                         return_d2=True)
            defect_migration_distance = d2[0][0] ** 0.5

        else:
            defect_migration_distance = None

    else:
        initial_center = np.array(defect_center)
        final_center = np.array(defect_center) - drift_frac_coords
        defect_migration_distance = None

    displacement_vectors = []
    displacement_norms = []
    for final_site, initial_site in zip(final_structure, initial_structure):

        # displacement_vectors are cartesian coordinates.
        # Return of pbc_shortest_vectors is a tuple of these two lists
        # * Array of displacement vectors from fcoords1 to fcoords2.
        # * Squared distances
        # Fo both, First index is fcoords1 index, and second is fcoords2 index
        displacement_vector, d2 = \
            pbc_shortest_vectors(lattice=initial_structure.lattice,
                                 fcoords1=initial_site.frac_coords,
                                 fcoords2=(final_site.frac_coords
                                           - drift_frac_coords),
                                 return_d2=True)
        displacement_vectors.append(list(displacement_vector[0][0]))
        displacement_norms.append(d2[0][0] ** 0.5)

    initial_vectors, initial_distances_2 = \
        pbc_shortest_vectors(lattice=initial_structure.lattice,
                             fcoords1=initial_center,
                             fcoords2=initial_structure.frac_coords,
                             return_d2=True)

    initial_distances = [d2 ** 0.5 for d2 in initial_distances_2[0]]

    final_vectors, final_distances_2 = \
        pbc_shortest_vectors(lattice=final_structure.lattice,
                             fcoords1=final_center,
                             fcoords2=final_structure.frac_coords,
                             return_d2=True)

    final_distances = [d2 ** 0.5 for d2 in final_distances_2[0]]

    # angles are nan when the displacements are zero or diverged.
    # [0] is needed for the variables with double loops.
    return {"initial_distances":          initial_distances,
            "final_distances":            final_distances,
            "displacement_vectors":       displacement_vectors,
            "displacement_norms":         displacement_norms,
            "initial_vectors":            initial_vectors[0].tolist(),
            "final_vectors":              final_vectors[0].tolist(),
            "defect_migration_distance":  defect_migration_distance}


def defect_center_from_coords(defect_coords: list,
                              structure: Structure) -> list:
    """Return defect center in fractional coordinates from given coordinates.

    Args:
        defect_coords (list):
        structure (Structure):

    Return:
        Averaged fractional coordinates along x-, y-, and z-directions .
    """
    # First defect_coords is used as a base point under periodic boundary
    # condition. Here, we are concerned about the case when two defect positions
    # are, e.g., [0.01, 0.01, 0.01] and [0.99, 0.99, 0.99].
    base = defect_coords[0]
    shortest_defect_coords = []

    for dc in defect_coords:
        diff = structure.lattice.get_distance_and_image(base, dc)[1]
        shortest_defect_coords.append(dc + diff)

    # np.array([[0, 0.1, 0.2], [0.3, 0.4, 0.5]]).transpose() =
    # np.array([[0, 0.3], [0.1, 0.4], [0.2, 0.5]])
    return [np.mean(i) for i in np.array(shortest_defect_coords).transpose()]


def distance_list(structure: Structure,
                  coords: np.array,
                  remove_self: bool = True) -> list:
    """Return a list of the shortest distances between a point and atoms

    Args:
       structure (Structure):
           pmg structure class object
       coords (1x3 numpy array):
           Fractional coordinates of a single point.
       remove_self (bool):
            Whether to remove self site distance, i.e. distance = 0.0

    Return:
        List of the shortest distances from coords to atoms
    """
    distance = []
    for frac_coord in structure.frac_coords:
        # return of get_distance_and_image is (distance, jimage).
        dist = structure.lattice.get_distance_and_image(frac_coord, coords)[0]
        if dist > 1e-5 or remove_self is False:
            distance.append(dist)

    return distance


def min_distance_from_coords(structure: Structure, coords: np.array) -> float:
    """Return the minimum distance

    Args:
       structure (Structure):
           pmg structure class object
       coords (1x3 numpy array):
           Fractional coordinates of a single point.

    Return:
        distance
    """
    return min(distance_list(structure, coords, remove_self=True))


def create_saturated_interstitial_structure(
        structure: Structure,
        inserted_atom_coords: list,
        dist_tol: float = 0.1,
        symprec: float = SYMMETRY_TOLERANCE,
        angle_tolerance: float = ANGLE_TOL) -> tuple:
    """ generates the sublattice for it based on the structure's space group.

    Original idea comes from pymatgen.
    pymatgen.analysis.defects.core.create_saturated_interstitial_structure

    Args:
        structure (Structure):
        inserted_atom_coords (list):
            List of 3x1 fractional coords
        dist_tol (float):
            changing distance tolerance of saturated structure,
            allowing for possibly overlapping sites
            but ensuring space group is maintained
        symprec (float):
        angle_tolerance (float):

    Returns:
        saturated_defect_structure (Structure):
            Saturated structure decorated with equivalent interstitial sites.
        atom_indices (list):
            The representative inserted atom indices.
        are_inserted (list):
            Show whether the inserted_atom_coords are inserted.
        symmetry_dataset (dict):
            Spglib style symmetry dataset with the various properties.
                + number:
                    International space group number
                + international:
                    International symbol
                + hall:
                    Hall symbol
                + transformation_matrix:
                    Transformation matrix from lattice of input cell to Bravais
                    lattice L^bravais = L^original * Tmat
                + origin shift:
                    Origin shift in the setting of "Bravais lattice" rotations,
                    translations: Rotation matrices and translation vectors.
                    Space group operations are obtained by
                    [(r,t) for r, t in zip(rotations, translations)]
                + wyckoffs:
                    Wyckoff letters
    """
    sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
    symmops = sga.get_symmetry_operations()
    saturated_structure = structure.copy()
    saturated_structure.DISTANCE_TOLERANCE = dist_tol

    atom_indices = []
    are_inserted = []
    for i, coord in enumerate(inserted_atom_coords):
        # Check whether the inserted atoms already exist.
        neighboring_atom_indices, distances = \
            get_neighboring_atom_indices(saturated_structure,
                                         coord, dist_tol)

        if neighboring_atom_indices:
            are_inserted.append(False)
            # If the inserted atom locates near other atoms within dist_tol,
            # first neighbor atom is set as inserted_atom_index.
            inserted_atom_index = neighboring_atom_indices[0]
            specie = saturated_structure[inserted_atom_index].specie

            logger.warning(f"Inserted position is too close to {specie}.\n  "
                           f"The distance is {distances[0]:5.3f} A.")

        else:
            are_inserted.append(True)
            # Get the last atom index to be inserted
            inserted_atom_index = len(saturated_structure)
            specie = DummySpecie()

            for symmop in symmops:
                new_interstitial_coords = symmop.operate(coord[:])
                try:
                    # validate_proximity (bool):
                    # Whether to check if inserted site is too close to an
                    # existing site. Defaults to False. If it is caught,
                    # ValueError is raised. For criterion, DISTANCE_TOLERANCE
                    # set to Structure is used. Here, this is used such that
                    # high-symmetric sites are reduced due to degeneracy.
                    saturated_structure.append(specie, new_interstitial_coords,
                                               validate_proximity=True)
                except ValueError:
                    pass

        atom_indices.append(inserted_atom_index)

    return saturated_structure, atom_indices, are_inserted


def get_neighboring_atom_indices(structure: Structure,
                                 coords: list,
                                 cutoff: float) -> Tuple[list, list]:
    """ Return the neighboring atom indices within the cutoff distance.
   
    Args:
        structure (Structure): Input structure.
        coords (list): Fractional coordinates of a centering point
        cutoff (float):
            Radius of a sphere in which atoms are considered as neighbors

    Return:
         Tuple of (neighboring indices, distances)
    """
    neighboring_indices = []
    distances = []
    for j, site in enumerate(structure):
        # returned tuple of (distance, periodic lattice translations)
        distance = structure.lattice.get_distance_and_image(coords,
                                                            site.frac_coords)
        if distance[0] < cutoff:
            neighboring_indices.append(j)
            distances.append(distance[0])

    return neighboring_indices, distances


def num_equivalent_clusters(structure: Structure,
                            inserted_atom_coords: Optional[list],
                            removed_atom_indices: Optional[list],
                            symprec: float = SYMMETRY_TOLERANCE,
                            angle_tolerance: float = ANGLE_TOL
                            ) -> Tuple[int, str]:
    """Calculate number of equivalent clusters in the structure.

    Args:
        structure (Structure):
            Supercell is assumed to big enough.
        inserted_atom_coords (list):
        removed_atom_indices (list):
            Needs to begin from 0.
        symprec (float):
        angle_tolerance (float):
            Angle tolerance in degree used for identifying the space group.

    Returns:
        Tuple of (num_equivalent_clusters (int), point_group (str))
    """
    inserted_atom_coords = inserted_atom_coords or []
    removed_atom_indices = removed_atom_indices or []

    sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
    num_symmop = len(sga.get_symmetry_operations())

    structure_with_cluster = structure.copy()
    for i in inserted_atom_coords:
        structure_with_cluster.append(DummySpecie(), i)
    structure_with_cluster.remove_sites(removed_atom_indices)

    sga_with_cluster = \
        SpacegroupAnalyzer(structure_with_cluster, symprec, angle_tolerance)
    sym_dataset = sga_with_cluster.get_symmetry_dataset()
    point_group = sym_dataset["pointgroup"]

    return int(num_symmop / num_symmetry_operation(point_group)), point_group


def first_appearing_index(structure: Structure,
                          element: Union[str, Element]) -> int:
    """Return first index where the specie appears. Return 0 if not exist

    If the element does not exist, 0 is returned.

    Args:
        structure: Input structure
        element: String of element or Element
    Return:
        Int of first index
    """
    if isinstance(element, Element):
        element = str(element)

    if element in structure.symbol_set:
        return min(structure.indices_from_symbol(element))
    else:
        return 0


def get_coordination_distances(structure: Structure,
                               atom_index: int,
                               cutoff: float) -> dict:
    """Calculated the coordination environment at the given atomic index.

    Args:
        structure (Structure):
            pmg Structure class object
        atom_index (int):
            The atomic index
        cutoff (float):
            Radius of a sphere in which atoms are considered as neighbors

    Return:
        Dict composed of specie name keys and distance list.
        E.g., {"Mg": [1.82, 1.82, 1.82, 1.82], ..}
    """
    neighbors = structure.get_neighbors(
        structure.sites[atom_index], cutoff, include_index=True)

    coordination_distances = defaultdict(list)

    for n in neighbors:
        # remove oxidation state from species_string, e.g., Mg2+ -> Mg
        specie = ''.join([i for i in str(n._species) if i.isalpha()])
        # float is needed as the numpy.float is not compatible with the
        # combination of OrderedDict and yaml. Otherwise, the interstitial.yaml
        # shows ugly printing.
        coordination_distances[specie].append(round(float(n[1]), 2))

    for distances in coordination_distances.values():
        distances.sort()

    return dict(coordination_distances)
