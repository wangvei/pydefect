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
                 displacements: dict,
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
            displacements (dict):
                see get_displacements for keys
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
        return cls(name=defect.name,
                   volume=defect.volume,
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

    def __repr__(self):
        outs = ["DefectStructure Summary"]
        return " ".join(outs)

    def show_displacements(self, all_atoms: bool = False):
        is_defect_center_atom = isinstance(self.defect_center, int)
        defect_center_str = [round(i, 3) for i in self.defect_center_coords]
        lines = [f"Is defect center atomic position?: {is_defect_center_atom}",
                 f"Defect center position: {defect_center_str}"]

        migration_distance = \
            round(self.displacements["defect_migration_distance"], 3)
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
                   defect_local_structure: Structure,
                   ltol: float,
                   stol: float):
        return self.final_local_structure.matches(
            other=defect_local_structure, ltol=ltol, stol=stol,
            primitive_cell=False, scale=False)
