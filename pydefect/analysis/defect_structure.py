# -*- coding: utf-8 -*-

from monty.json import MSONable

from pymatgen.core.structure import Structure


def get_neighbors_with_distance(structure, coords, neighbor_tolerance=1.1,
                                initial_cutoff=2, is_frac=True):
    """
    Args:
        structure (Structure):
        coords (3x1 array):
        neighbor_tolerance (float):
            get sites whose displacement_distance < min_distance * neighbor_tolerance
        initial_cutoff (float): initial cutoff displacement_distance when searching neighbors
        is_frac (bool): Is coordinate fractional, not Cartesian
    Returns:
        (list of (Site, displacement_distance))
    """

    # conversion to Cartesian
    if is_frac:
        cart_coords = structure.lattice.get_cartesian_coords(coords)
    else:
        cart_coords = coords

    self_threshold = 1e-2

    def exclude_self(distance_list):
        return [site_dist for site_dist in distance_list
                if site_dist[1] > self_threshold]

    # search neighbor
    cutoff = initial_cutoff
    sites_dist = []
    while not exclude_self(sites_dist):
        cutoff = cutoff * 1.2
        sites_dist = structure.get_sites_in_sphere(cart_coords, cutoff)

    min_distance = min(distance for _, distance in exclude_self(sites_dist))
    max_distance = min_distance * neighbor_tolerance
    candidates = structure.get_sites_in_sphere(cart_coords, max_distance)
    sites_dist = [site_dist for site_dist in candidates
                  if site_dist[1] >= min_distance]

    return sites_dist


def get_neighbors(structure, coords, neighbor_tolerance,
                  initial_cutoff=3, is_frac=True, include_self=False):
    """
    Args:
        structure (Structure):
        coords (3x1 array):
        neighbor_tolerance (float):
            get sites whose displacement_distance < min_distance * neighbor_tolerance
        initial_cutoff (float):
            initial cutoff displacement_distance when searching neighbors
        is_frac (bool):
            Whether the coordinate is in fractional instead of cartesian
        include_self (bool):
            Whether the sites whose displacement_distance < 1e-2 are included.
    Returns:
        (list of Site)
    """
    return [site for site, _ in get_neighbors_with_distance(
        structure, coords, neighbor_tolerance,
        initial_cutoff=initial_cutoff, is_frac=is_frac)]


class DefectStructure(MSONable):
    """
    A class related to eigenvalues in a defect calculation.
    """
    def __init__(self,
                 name: str,
                 volume: float,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 final_structure: Structure,
                 displacements: list,
                 initial_site_symmetry: str,
                 final_site_symmetry: str,
                 initial_symmetry_multiplicity: int,
                 final_symmetry_multiplicity: int,
                 perturbed_sites: list,
                 num_equiv_sites: int):
        """
        Args:
            name (str):
                Name of a defect
            volume (float):
                Volume of the supercell.
            initial_structure (Structure):
                Structure with a defect before the structure optimization.
            perturbed_initial_structure (Structure):
                Initial structure with perturbation of neighboring atoms.
            final_structure (Structure):
                pmg Structure class object. Usually relaxed structures
            displacements (list):
                Displacements
            initial_site_symmetry (str):
                Initial site symmetry such as D4h.
            final_site_symmetry (str):
                Final site symmetry.
            initial_symmetry_multiplicity (int):
                Symmetry multiplicity in the initial structure.
            initial_symmetry_multiplicity (int):
                Symmetry multiplicity in the final structure.
            perturbed_sites (list):
                Indices of the perturbed site for reducing the symmetry
            num_equiv_sites (int):
                Number of equivalent sites in the given supercell.

        """
        self.name = name
        self.volume = volume
        self.initial_structure = initial_structure
        self.perturbed_initial_structure = perturbed_initial_structure
        self.final_structure = final_structure
        self.initial_site_symmetry = initial_site_symmetry
        self.final_site_symmetry = final_site_symmetry
        self.displacements = displacements
        self.initial_symmetry_multiplicity = initial_symmetry_multiplicity
        self.final_symmetry_multiplicity = final_symmetry_multiplicity
        self.perturbed_sites = perturbed_sites
        self.num_equiv_sites = num_equiv_sites



