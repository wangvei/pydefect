# -*- coding: utf-8 -*-

from copy import deepcopy
from monty.json import MSONable

from pymatgen.core.structure import Structure

from pydefect.analysis.defect import Defect
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


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
        return [d for d in distance_list if d[1] > self_threshold]

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
    """ A class related to local structures around defect """

    def __init__(self,
                 name: str,
                 volume: float,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 final_structure: Structure,
                 displacements: list,
                 initial_site_symmetry: str,
                 final_site_symmetry: str,
                 neighboring_sites: list):
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
            neighboring_sites (list):
                Indices of the perturbed site for reducing the symmetry

        """
        self.name = name
        self.volume = volume

        self.initial_structure = deepcopy(initial_structure)
        self.perturbed_initial_structure = deepcopy(perturbed_initial_structure)
        self.final_structure = deepcopy(final_structure)

        self.initial_site_symmetry = initial_site_symmetry
        self.final_site_symmetry = final_site_symmetry

        self.displacements = displacements

        self.neighboring_sites = neighboring_sites

        removed_sites = []
        for i in range(len(self.final_structure)):
            if i not in self.neighboring_sites:
                removed_sites.append(i)

        self.initial_local_structure = deepcopy(initial_structure)
        self.initial_local_structure.remove_sites(removed_sites)
        self.final_local_structure = deepcopy(final_structure)
        self.final_local_structure.remove_sites(removed_sites)

    @classmethod
    def from_defect(cls, defect: Defect):
        """ Construct class obj from Defect object.

        Args:
            defect (Defect):
                namedtuple("Defect",
                           "defect_entry", "dft_results", "correction")
        """

        return cls(name=defect.defect_entry.name,
                   volume=defect.dft_results.volume,
                   initial_structure=defect.defect_entry.initial_structure,
                   perturbed_initial_structure=
                   defect.defect_entry.perturbed_initial_structure,
                   final_structure=defect.dft_results.final_structure,
                   initial_site_symmetry=
                   defect.defect_entry.initial_site_symmetry,
                   final_site_symmetry=defect.dft_results.site_symmetry,
                   displacements=defect.dft_results.displacements,
                   neighboring_sites=defect.defect_entry.neighboring_sites)

    def __repr__(self):
        outs = ["DefectStructure Summary"]
        return " ".join(outs)

    def show_displacements(self, all_atoms=False):
        lines = ["element initial   final   disp  angle            ",
                 " name   dist(A)  dist(A)  dist  (deg)   direction"]
        if all_atoms:
            candidate = range(len(self.final_structure))
        else:
            candidate = self.neighboring_sites
        for s in candidate:

            element = str(self.final_structure[s].specie)
            try:
                angle = int(round(self.displacements[4][s], 0))
                if angle < 60:
                    direction = "inward"
                elif angle > 120:
                    direction = "outward"
                else:
                    direction = "tangent"
                angle = str(angle)
            except TypeError:
                angle = "None"
                direction = "None"

            lines.append(" {:>5}    {:3.2f}     {:3.2f}   {:4.2f}  {:>4} {:>12}".format(
                element,
                round(self.displacements[0][s], 3),
                round(self.displacements[1][s], 3),
                round(self.displacements[3][s], 3), angle, direction))
        return "\n".join(lines)

    def comparator(self,
                   defect_local_structure: Structure,
                   ltol: float,
                   stol: float):
        return self.final_local_structure.matches(
            other=defect_local_structure, ltol=ltol, stol=stol,
            primitive_cell=False, scale=False)
