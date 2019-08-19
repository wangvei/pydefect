# -*- coding: utf-8 -*-
from copy import deepcopy
from enum import Enum, unique
import json
import numpy as np
import os
import ruamel.yaml as yaml
from typing import Optional, Tuple
from pathlib import Path
import shutil

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.error_classes import StructureError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import (
    count_equivalent_clusters, defect_center_from_coords, distance_list)
from pydefect.util.vasp_util import element_diff_from_structures
from pydefect.core.config import (
    DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL, CUTOFF_RADIUS)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

# When inheriting Enum class, MSONable does not work well.
@unique
class DefectType(Enum):
    vacancy = "vacancy"
    substituted = "substituted"
    interstitial = "interstitial"
    complex = "complex"

    def __repr__(self):
        return self.value

    # This is a must.
    def __str__(self):
        return self.value

    @classmethod
    def from_string(cls, s: str):
        for m in cls:
            if m.value == s or m.name == s:
                return m
        raise AttributeError("Defect type " + str(s) + " is not proper.\n" +
                             "Supported info:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.name for e in cls])

    @property
    def is_defect_center_atom(self):
        return self in [DefectType.substituted, DefectType.interstitial]


def determine_defect_type(inserted_atoms: list,
                          removed_atoms: list) -> DefectType:
    d = DefectType.complex

    if len(removed_atoms) == 1:
        if len(inserted_atoms) == 1 and \
                removed_atoms[0]["coords"] == inserted_atoms[0]["coords"]:
            d = DefectType.substituted
        elif len(inserted_atoms) == 0:
            d = DefectType.vacancy
    elif len(removed_atoms) == 0 and len(inserted_atoms) == 1:
        d = DefectType.interstitial

    return d


class DefectEntry(MSONable):
    """ Holds information related to initial setting of a single defect. """

    def __init__(self,
                 name: str,
                 defect_type: DefectType,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 removed_atoms: list,
                 inserted_atoms: list,
                 changes_of_num_elements: dict,
                 charge: int,
                 initial_site_symmetry: str,
                 cutoff: float,
                 neighboring_sites: list,
                 annotation: Optional[str] = None,
                 num_equiv_sites: Optional[int] = None):
        """
        Args:
            name (str):
                Name of a defect without charge.
            defect_type (DefectType):
                Defect type defined in DefectType enumeration.
            initial_structure (Structure):
                Structure with a defect before the structure optimization.
            perturbed_initial_structure (Structure):
                Initial structure with perturbation of neighboring atoms.
            removed_atoms (list):
                List of dict with the following values.
                + "element" (str):
                + "index" (int): Removed index in the original supercell.
                + "coords" (list):
            inserted_atoms (list):
                List of dict with the following values.
                + "element" (str):
                + "index" (int): Removed index in the defect supercell.
                + "coords" (list):
            changes_of_num_elements (dict):
                Keys: Element names
                Values: Change of the numbers of elements wrt perfect supercell.
            charge (int):
                Defect charge state. Charge is also added to the structure.
            initial_site_symmetry (str):
                Initial site symmetry such as D4h.
            cutoff (float):
                Cutoff Radius determining the neighboring sites.
            neighboring_sites (list):
                Atomic indices of the neighboring sites within cutoff radius.
            annotation (str):
                Annotation used for distinguishing defects with same
                net atom exchange but different local structures.
            num_equiv_sites (int):
                Number of equivalent sites in the given structure.
        """
        self.name = name
        self.defect_type = defect_type
        self.initial_structure = initial_structure
        self.perturbed_initial_structure = perturbed_initial_structure
        self.removed_atoms = deepcopy(removed_atoms)
        self.inserted_atoms = deepcopy(inserted_atoms)
        self.changes_of_num_elements = deepcopy(changes_of_num_elements)
        self.charge = charge
        self.initial_site_symmetry = initial_site_symmetry
        self.cutoff = cutoff
        if not neighboring_sites:
            raise StructureError(
                "No neighboring site detected, so increase cutoff radius.")
        self.neighboring_sites = neighboring_sites[:]
        self.annotation = annotation
        self.num_equiv_sites = num_equiv_sites

    def __repr__(self):
        annotation = "" if self.annotation is None else self.annotation

        outs = ["name: " + str(self.name),
                "defect type: " + str(self.defect_type),
                "charge: " + str(self.charge),
                "annotation: " + annotation,
                "perturbed initial structure: \n" +
                str(self.perturbed_initial_structure),
                "initial site symmetry: " + str(self.initial_site_symmetry),
                "removed_atoms: " + str(self.removed_atoms),
                "inserted atoms: " + str(self.inserted_atoms),
                "changes of num element: " + str(self.changes_of_num_elements),
                "cut off radius: " + str(self.cutoff),
                "neighboring sites: " + str(self.neighboring_sites),
                "num_equiv sites: " + str(self.num_equiv_sites)]
        return "\n".join(outs)

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_dict(cls, d: dict):
        # The keys need to be converted to integers.
        defect_type = DefectType.from_string(d["defect_type"])

#        removed_atoms = {int(k): v for k, v in d["removed_atoms"].items()}
#        inserted_atoms = {int(k): v for k, v in d["inserted_atoms"].items()}

        initial_structure = d["initial_structure"]
        if isinstance(initial_structure, dict):
            initial_structure = Structure.from_dict(initial_structure)

        perturbed_initial_structure = d["perturbed_initial_structure"]
        if isinstance(perturbed_initial_structure, dict):
            perturbed_initial_structure = \
                Structure.from_dict(perturbed_initial_structure)

        return cls(
            name=d["name"],
            defect_type=defect_type,
            initial_structure=initial_structure,
            perturbed_initial_structure=perturbed_initial_structure,
            removed_atoms=d["removed_atoms"],
            inserted_atoms=d["inserted_atoms"],
            changes_of_num_elements=d["changes_of_num_elements"],
            charge=d["charge"],
            initial_site_symmetry=d["initial_site_symmetry"],
            cutoff=d["cutoff"],
            neighboring_sites=d["neighboring_sites"],
            annotation=d["annotation"],
            num_equiv_sites=d["num_equiv_sites"])

    def as_dict(self) -> dict:
        defect_type = str(self.defect_type)

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "name": self.name,
             "defect_type": defect_type,
             "initial_structure": self.initial_structure,
             "perturbed_initial_structure": self.perturbed_initial_structure,
             "removed_atoms": self.removed_atoms,
             "inserted_atoms": self.inserted_atoms,
             "changes_of_num_elements": self.changes_of_num_elements,
             "charge": self.charge,
             "initial_site_symmetry": self.initial_site_symmetry,
             "cutoff": self.cutoff,
             "neighboring_sites": self.neighboring_sites,
             "annotation": self.annotation,
             "num_equiv_sites": self.num_equiv_sites}

        return d

    @classmethod
    def from_yaml(cls,
                  yaml_filename: Optional[str] = None,
                  disp_dist: float = 0.2,
                  symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                  angle_tolerance: float = ANGLE_TOL,
                  cutoff: float = CUTOFF_RADIUS,
                  calc_num_equiv_site: bool = True,
                  defect_name: str = None):
        """Construct the DefectEntry object from perfect and defective POSCARs.

        Note1: displacement_distance needs to be the same as the twice of max
               displacement distance used for reducing the symmetry.
        Note2: Only unrelaxed but perturbed defective POSCAR structure is
               accepted as the inserted atoms are assumed not to be moved from
               their original positions.

        filename (str):
            yaml filename.
        displacement_distance (float):
            Tolerance to judge whether the atoms are the same between the given
            defective and perfect supercell structures.
        symprec (float):
            Length tolerance in Angstrom used for identifying the space group.
        angle_tolerance (float):
            Angle tolerance used for identifying the space group.
        cutoff (str):
            Radius of sphere in which atoms are considered as neighbors of a
            defect.
        defect_name (str):
            Defect name such as "Va_O1_2_inward", "Mg_i+Va_O1*2_2_coord1".
            Although directory name is usually used, this is needed for e.g.,
            unittest.

        An example of the yaml file.
            defect_structure: POSCAR
            perfect_structure: ../../defects/perfect/POSCAR
            displacement_distance (optional): 0.2
        """
        if yaml_filename is None:
            org_yaml = \
                Path(os.path.dirname(__file__)) / "default_defect_entry.yaml"
            yaml_filename = "defect_entry.yaml"
            shutil.copyfile(org_yaml, yaml_filename)

        abs_dir = Path(yaml_filename).parent

        with open(yaml_filename, "r") as f:
            yaml_data = yaml.safe_load(f)

        disp_dist = yaml_data.get("displacement_distance", disp_dist)
        # Perfect and perturbed defect structures.
        perfect_structure = \
            Structure.from_file(abs_dir / yaml_data["perfect_structure"])
        defect_structure = \
            Structure.from_file(abs_dir / yaml_data["defect_structure"])
        element_diff = element_diff_from_structures(
            defect_structure, perfect_structure)

        if not defect_name:
            defect_name = yaml_data.get("name", None)
            if not defect_name:
                _, defect_name = os.path.split(os.getcwd())
        name, charge, annotation = divide_dirname(defect_name)

        inserted_atom_indices = [i for i in range(defect_structure.num_sites)]
        removed_atoms = []
        for i, p_site in enumerate(perfect_structure):
            for j in inserted_atom_indices:
                d_site = defect_structure[j]
                distance = p_site.distance(d_site)
                # check displacement_distance and species for comparison
                if distance < disp_dist and p_site.specie == d_site.specie:
                    inserted_atom_indices.remove(j)
                    break
            else:
                removed_atoms.append({"element": str(p_site.specie),
                                      "index": i,
                                      "coords": list(p_site.frac_coords)})

        # check the consistency of the removed and inserted atoms
        if len(perfect_structure) + len(inserted_atom_indices) \
                - len(removed_atoms) != len(defect_structure):
            raise StructureError(
                "Atoms are not properly mapped within the displacement.")

        # pristine_defect_structure is a defect structure without displacement.
        pristine_defect_structure = deepcopy(perfect_structure)
        removed_atom_indices = [i["index"] for i in removed_atoms]
        removed_atom_coords = [i["coords"] for i in removed_atoms]
        for r in sorted(removed_atom_indices, reverse=True):
            pristine_defect_structure.pop(r)

        inserted_atoms = []
        for i in sorted(inserted_atom_indices):
            inserted_atom = defect_structure[i]
            inserted_atoms.append({"element": str(inserted_atom.specie),
                                   "index": i,
                                   "coords": list(inserted_atom.frac_coords)})
            pristine_defect_structure.insert(i, inserted_atom.specie,
                                             inserted_atom.frac_coords)

        inserted_atom_coords = [i["coords"] for i in inserted_atoms]
        defect_center_coords = cls.calc_defect_center(removed_atom_coords,
                                                      inserted_atom_coords,
                                                      defect_structure)
        neighboring_sites = []
        for i, site in enumerate(pristine_defect_structure):
            distance, _ = \
                site.distance_and_image_from_frac_coords(defect_center_coords)
            # Defect itself is not included to the neighboring sites
            if 1e-3 < distance < cutoff:
                neighboring_sites.append(i)

        sga = SpacegroupAnalyzer(
            pristine_defect_structure, symprec, angle_tolerance)
        initial_site_symmetry = sga.get_point_group_symbol()

        num_equiv_sites = None
        if calc_num_equiv_site:
            num_equiv_sites, _ = count_equivalent_clusters(perfect_structure,
                                                           inserted_atom_coords,
                                                           removed_atom_indices,
                                                           disp_dist)
        pristine_defect_structure.set_charge(charge)
        defect_structure.set_charge(charge)
        defect_type = determine_defect_type(inserted_atoms, removed_atoms)

        return cls(name=name,
                   defect_type=defect_type,
                   initial_structure=pristine_defect_structure,
                   perturbed_initial_structure=defect_structure,
                   removed_atoms=removed_atoms,
                   inserted_atoms=inserted_atoms,
                   changes_of_num_elements=element_diff,
                   charge=charge,
                   initial_site_symmetry=initial_site_symmetry,
                   cutoff=cutoff,
                   neighboring_sites=neighboring_sites,
                   annotation=annotation,
                   num_equiv_sites=num_equiv_sites)

    @classmethod
    def load_json(cls, filename="defect_entry.json"):
        return loadfn(filename)

    @property
    def atom_mapping_to_perfect(self):
        """ Returns a list of atom mapping from defect structure to perfect.

        Example of Mg32O32 supercell:
            When 33th atom, namely first O, is removed,
                return [0, 1, 2, .., 31, 33, 34, .., 62] (=mapping)
                len(mapping) = 63
        """
        total_nions_in_perfect = (len(self.initial_structure)
                                  - len(self.inserted_atoms)
                                  + len(self.removed_atoms))
        # initial atom mapping.
        mapping = list(range(total_nions_in_perfect))

        removed_atom_indices = [i["index"] for i in self.removed_atoms]
        for o in sorted(removed_atom_indices, reverse=True):
            mapping.pop(o)

        inserted_indices = [i["index"] for i in self.inserted_atoms]
        for i in sorted(inserted_indices, reverse=True):
            mapping.insert(i, None)

        return mapping

    def to_json_file(self, filename="defect_entry.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @staticmethod
    def calc_defect_center(removed_atom_coords: list,
                           inserted_atom_coords: list,
                           structure: Structure) -> list:
        """ Calculates arithmetic average to estimate center in frac coords."""
        defect_coords = removed_atom_coords + inserted_atom_coords
        return defect_center_from_coords(defect_coords, structure)

    @property
    def defect_center_coords(self) -> list:
        """ Return fractional coordinates of the defect center. """
        removed_atom_coords = [i["coords"] for i in self.removed_atoms]
        inserted_atom_coords = [i["coords"] for i in self.inserted_atoms]

        return self.calc_defect_center(
            removed_atom_coords, inserted_atom_coords, self.initial_structure)

    @property
    def anchor_atom_index(self) -> int:
        """ Returns an index of atom that is the farthest from the defect. """
        # distance_set = \
        #     self.initial_structure.lattice.get_all_distances(
        #         self.defect_center, self.initial_structure.frac_coords)[0]

        return anchor_atom_index(structure=self.initial_structure,
                                 center=self.defect_center_coords)


def anchor_atom_index(structure: Structure, center: np.array) -> int:
    """ Returns an index of atom that is the farthest from the defect.

    This atom is assumed not to displace in the defective supercell, and
    so used for analyzing local structure around a defect.
    Note that only the first occurrence is returned when using argmax.
    docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html
    """
    distance_set = \
        structure.lattice.get_all_distances(center, structure.frac_coords)[0]

    return int(np.argmax(distance_set))


def divide_dirname(dirname: str) -> Tuple[str, int, Optional[str]]:
    """Return the divided dirname to name, charge, and annotation

    There must be only one digit in the split dirname, which is a charge state,
    and front and back part correspond to name and annotation, respectively.

    "Va_Mg1_2" -> name = "Va_Mg1", charge = 2, annotation = None
    "Va_O1_2_inward" -> name = "Va_O1", charge = 2, annotation = "inward"
    "Mg_i+Va_O1*2_2_coord1"
            -> name = "Mg_i+Va_O1*2", charge = 2, annotation = "coord1"
    """
    split_dirname = dirname.split("_")

    def is_digit(n):
        try:
            int(n)
            return True
        except ValueError:
            return False

    digit_positions = [x for x, y in enumerate(split_dirname) if is_digit(y)]

    if len(digit_positions) != 1:
        raise ValueError(f"The dirname {dirname} is not valid")
    else:
        digit_pos = digit_positions[0]
        name = "_".join(split_dirname[:digit_pos])
        charge = int(split_dirname[digit_pos])
        annotation = "_".join(split_dirname[digit_pos + 1:])
        annotation = annotation if annotation else None

    return name, charge, annotation


def distances_from_defect_center(structure: Structure,
                                 defect_entry: DefectEntry) -> list:
    """ Returns a list of distances at atomic sites from defect center

    Note that in the case of an interstitial-type defect, zero is also set for
    the interstitial site itself.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        defect_entry (DefectEntry):
            DefectEntry class object considered
    """
    return distance_list(structure, defect_entry.defect_center_coords)