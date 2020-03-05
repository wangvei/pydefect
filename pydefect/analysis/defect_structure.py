# -*- coding: utf-8 -*-

from copy import deepcopy
from itertools import groupby
from typing import Union, List
from operator import attrgetter

from monty.json import MSONable
from pydefect.analysis.defect import Defect
from pydefect.util.logger import get_logger
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure import Structure

logger = get_logger(__name__)

"""
This module provides a class used to analyze the atomic structures around 
defects.
"""


class DefectStructure(MSONable):
    """ A class related to local structures around defect

    Attributes:
        initial_local_structure (Structure)
            Initial local structure with only neighboring atoms where the
            neighbors are determined at the final structure.
        final_local_structure (Structure)
            Same as initial_local_structure but for the final structure.
    """

    def __init__(self,
                 name: str,
                 charge: int,
                 final_volume: float,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 final_structure: Structure,
                 initial_site_symmetry: str,
                 final_site_symmetry: str,
                 displacements: dict,
                 defect_center: Union[int, list],
                 defect_center_coords: list,
                 neighboring_sites: list):
        """
        Args:
            name (str):
                Name of a defect
            final_volume (float):
                Volume of the supercell.
            initial_structure (Structure):
                Structure with a defect before the structure optimization.
            perturbed_initial_structure (Structure):
                Initial structure with perturbation of neighboring atoms.
            final_structure (Structure):
                pmg Structure class object. Usually relaxed structures
            initial_site_symmetry (str):
                Initial site symmetry such as D4h.
            final_site_symmetry (str):
                Final site symmetry.
            displacements (dict):
                see get_displacements for keys
                Keys:
                + initial_distances:
                    Distances from a defect in the initial structure.
                + final_distances:
                    Distances from a defect in the final structure.
                + displacement_vectors:
                    Displacement vectors of atoms from initial to final
                    structures.
                + displacement_norms:
                    Norms of displacement vectors.
                + initial_vectors:
                    Vectors from a defect position to atoms in the initial
                    supercell.
                + final_vectors:
                    Vectors from a defect position to atoms in the final
                    supercell.
                + defect_migration_distance:
                    Distance the defect migrates defined only for interstitials,
                    antisites, and substituted defects.
            defect_center (int/list):
                Show the defect center. When the defect center is an
                atomic position, atom index number is set. Contrary, when
                the defect center is defined by the initial position, such as
                vacancies and complex defects, the fractional coordinates are
                set.
            defect_center_coords (list):
                Fractional coordinates of defect center.
            neighboring_sites (list):
                Atomic indices of the neighboring sites.
                If is_defect_center_atom is True, neighboring_sites_after_relax
                is used. Otherwise, neighboring_sites in DefectEntry is used
                because the defect position is not well defined after the
                structure is relaxed.

        """
        self.name = name
        self.charge = charge
        self.final_volume = final_volume

        self.initial_structure = deepcopy(initial_structure)
        self.perturbed_initial_structure = \
            deepcopy(perturbed_initial_structure)
        self.final_structure = deepcopy(final_structure)

        self.initial_site_symmetry = initial_site_symmetry
        self.final_site_symmetry = final_site_symmetry

        self.displacements = displacements
        self.defect_center = defect_center
        self.defect_center_coords = defect_center_coords
        self.neighboring_sites = neighboring_sites

        # distant sites from the defect in the *final_structure*.
        distant_sites = [i for i in range(len(self.final_structure))
                         if i not in self.neighboring_sites]

        self.initial_local_structure = deepcopy(initial_structure)
        self.initial_local_structure.remove_sites(distant_sites)
        self.final_local_structure = deepcopy(final_structure)
        self.final_local_structure.remove_sites(distant_sites)

    @classmethod
    def from_defect(cls, defect: Defect):
        """ Retrieve from Defect class object"""
        return cls(name=defect.name,
                   charge=defect.charge,
                   final_volume=defect.final_volume,
                   initial_structure=defect.initial_structure,
                   perturbed_initial_structure=
                   defect.perturbed_initial_structure,
                   final_structure=defect.final_structure,
                   initial_site_symmetry=defect.initial_symmetry,
                   final_site_symmetry=defect.final_symmetry,
                   displacements=defect.displacements,
                   defect_center=defect.defect_center,
                   defect_center_coords=defect.defect_center_coords,
                   neighboring_sites=defect.neighboring_sites)

    def show_displacements(self, all_atoms: bool = False) -> str:
        is_defect_center_atom = isinstance(self.defect_center, int)
        defect_center_str = [round(i, 3) for i in self.defect_center_coords]
        lines = [f"Is defect center atomic position?: {is_defect_center_atom}",
                 f"Defect center position: {defect_center_str}"]

        if isinstance(self.displacements["defect_migration_distance"], float):
            migration_distance = \
                round(self.displacements["defect_migration_distance"], 3)
        else:
            migration_distance = None

        if migration_distance:
            lines.append(f"Defect traveling distance: {migration_distance}")

        lines.extend(
            [f"Site symmetry: "
             f"{self.final_site_symmetry} <- {self.initial_site_symmetry}",
             f"    element  final <-initial   disp",
             f"index name   dist(A)  dist(A)  dist"
             f"   coordination (final) <- coordination (initial)"])

        if all_atoms:
            candidate = range(len(self.final_structure))
        else:
            candidate = self.neighboring_sites

        for i in candidate:
            element = str(self.final_structure[i].specie)
            i_d = round(self.displacements["initial_distances"][i], 3)
            f_d = round(self.displacements["final_distances"][i], 3)
            d_n = round(self.displacements["displacement_norms"][i], 3)
            f_v = " ".join(["{:6.2f}".format(round(i, 2))
                            for i in self.displacements["final_vectors"][i]])
            i_v = " ".join(["{:6.2f}".format(round(i, 2))
                            for i in self.displacements["initial_vectors"][i]])

            lines.append(f"{i:>4} {element:>5}   {f_d:5.2f} <- {i_d:5.2f}    "
                         f"{d_n:4.2f}  {f_v}  <- {i_v}")

        return "\n".join(lines)

    def comparator(self,
                   defect_structure,
                   stol: float = 0.03,
                   **kwargs) -> bool:

        compared_structure = defect_structure.final_local_structure
        return self.final_local_structure.matches(other=compared_structure,
                                                  stol=stol,
                                                  primitive_cell=False,
                                                  scale=False,
                                                  **kwargs)


def defect_structure_matcher(d_list: List[DefectStructure],
                             site_tolerance: float = 0.05) -> dict:
    """ A list of DefectStructure is grouped by structural equality.

    Args:
        d_list ([DefectStructure]): List of DefectStructure to be grouped
        site_tolerance: Defined as the fraction of the average free length per
            atom := ( V / Nsites ) ** (1/3) Default in pymatgen is 0.3.
            See docstrings of StructureMatcher in
            pymatgen.analysis.structure_matcher.
    Returns:
        A list of lists of matched structures
        Assumption: if s1 == s2 but s1 != s3, than s2 and s3 will be put
        in different groups without comparison.
    """
    sm = StructureMatcher(stol=site_tolerance,
                          primitive_cell=False,
                          scale=False)
    group = {}
    # Pre-grouped by name
    d_list.sort(key=attrgetter("name"))
    for name, g in groupby(d_list, key=attrgetter("name")):
        unmatched = list(g)
        group[name] = []
        while len(unmatched) > 0:
            d = unmatched.pop(0)
            matches = [d.charge]
            ref = d.final_structure
            inds = filter(lambda i: sm.fit(ref, unmatched[i].final_structure),
                          list(range(len(unmatched))))
            # filter returns generator, so needs to change to list
            inds = list(inds)
            matches.extend([unmatched[i].charge for i in inds])
            unmatched = [unmatched[i] for i in range(len(unmatched))
                         if i not in inds]
            group[name].append(matches)

    return group

