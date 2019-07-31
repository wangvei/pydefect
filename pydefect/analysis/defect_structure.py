# -*- coding: utf-8 -*-

from copy import deepcopy
from monty.json import MSONable
from typing import Union

from pymatgen.core.structure import Structure

from pydefect.analysis.defect import Defect
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class DefectStructure(MSONable):
    """ A class related to local structures around defect """

    def __init__(self,
                 name: str,
                 volume: float,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 final_structure: Structure,
                 initial_site_symmetry: str,
                 final_site_symmetry: str,
                 displacements: list,
                 defect_center: Union[int, list],
                 defect_center_coords: list,
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
            initial_site_symmetry (str):
                Initial site symmetry such as D4h.
            final_site_symmetry (str):
                Final site symmetry.
            displacements (list):
                Displacements composed of 5 quantities.
                [initial_distances, final_distances, displacement_vectors,
                displacement_distances, angles_wrt_the_defect_site]
            defect_center (int/list):
                Show the defect center. When the defect center is an
                atomic position, atom index number is set. Contrary, when
                the defect center is defined by the initial position, such as
                vacancies and complex defects, the fractional coordinates are
                set.
            defect_center_coords (list):
            neighboring_sites (list):
                Atomic indices of the neighboring sites.
                If is_defect_center_atom is True, neighboring_sites_after_relax
                is used. Otherwise, neighboring_sites in DefectEntry is used
                because the defect position is not defined after the structure
                is relaxed.

        """
        self.name = name
        self.volume = volume

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

        removed_sites = [i for i in range(len(self.final_structure))
                         if i not in self.neighboring_sites]

        self.initial_local_structure = deepcopy(initial_structure)
        self.initial_local_structure.remove_sites(removed_sites)
        self.final_local_structure = deepcopy(final_structure)
        self.final_local_structure.remove_sites(removed_sites)

    @classmethod
    def from_defect(cls, defect: Defect):

        defect_center = defect.dft_results.defect_center
        if isinstance(defect_center, int):
            neighboring_sites = \
                defect.dft_results.neighboring_sites_after_relax
            defect_center_coords = \
                defect.dft_results.final_structure[defect_center].frac_coords
        else:
            neighboring_sites = defect.defect_entry.neighboring_sites
            defect_center_coords = defect.defect_entry.defect_center_coords

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
                   defect_center=defect_center,
                   defect_center_coords=defect_center_coords,
                   neighboring_sites=neighboring_sites)

    def __repr__(self):
        outs = ["DefectStructure Summary"]
        return " ".join(outs)

    def show_displacements(self, all_atoms: bool = False):
        is_defect_center_atom = \
            True if isinstance(self.defect_center, int) else False

        defect_center_coords = [round(i, 3) for i in self.defect_center_coords]
#        defect_index = self.defect_center

        lines = [f"Is defect center atomic position?: {is_defect_center_atom}",
                 f"Defect center position: {defect_center_coords}",
                 f"Site symmetry: "
                 f"{self.final_site_symmetry} <- {self.initial_site_symmetry}",
                 f"    element  final <-initial   disp  angle            ",
                 f"index name   dist(A)  dist(A)  dist  (deg)   direction"
                 f"   coordination (final) <- coordination (initial)"]
        if all_atoms:
            candidate = range(len(self.final_structure))
        else:
            candidate = self.neighboring_sites

        for s in candidate:
            element = str(self.final_structure[s].specie)
            print(self.final_structure[s].frac_coords)
            initial_distance = round(self.displacements[0][s], 3)
            final_distance = round(self.displacements[1][s], 3)
            displacement_distance = round(self.displacements[3][s], 3)
            final_vector = " ".join(["{:6.2f}".format(round(i, 2))
                                     for i in self.displacements[6][s]])
            initial_vector = " ".join(["{:6.2f}".format(round(i, 2))
                                       for i in self.displacements[5][s]])

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

            lines.append("{:>4} {:>5}   {:5.2f} <- {:5.2f}   {:4.2f}  {:>4} "
                         "{:>12}   {} <- {}".format(s, element, final_distance, initial_distance,
                                         displacement_distance,
                                         angle, direction, final_vector, initial_vector))
        return "\n".join(lines)

    def comparator(self,
                   defect_local_structure: Structure,
                   ltol: float,
                   stol: float):
        return self.final_local_structure.matches(
            other=defect_local_structure, ltol=ltol, stol=stol,
            primitive_cell=False, scale=False)
