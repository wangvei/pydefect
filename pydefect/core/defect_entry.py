# -*- coding: utf-8 -*-
from copy import deepcopy
from enum import Enum, unique
import json
import numpy as np
import os
import ruamel.yaml as yaml
from typing import Union

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.defect_name import SimpleDefectName
from pydefect.core.error_classes import StructureError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import count_equivalent_clusters
from pydefect.vasp_util.util import element_diff_from_structures, \
    get_defect_charge_from_vasp
from pydefect.util.structure_tools import defect_center_from_coords
from pydefect.core.config \
    import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL, CUTOFF_RADIUS

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


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
    def from_string(cls, s):
        for m in cls:
            if m.value == s or m.name == s:
                return m
        raise AttributeError("Defect type " + str(s) + " is not proper.\n" +
                             "Supported info:\n" + cls.name_list())

    @classmethod
    def from_simple_defect_name(cls, name: Union[str, SimpleDefectName]):
        if isinstance(name, str):
            name = SimpleDefectName.from_str(name)

        if name.is_interstitial:
            return cls.interstitial
        elif name.is_vacancy:
            return cls.vacancy
        else:
            return cls.substituted

    @classmethod
    def name_list(cls):
        return ', '.join([e.name for e in cls])

    @property
    def is_defect_center_atom(self):
        return self in [DefectType.substituted, DefectType.interstitial]


def determine_defect_type(inserted_atoms: dict, removed_atoms: dict):
    d = DefectType.complex

    if len(removed_atoms) == 1:
        if len(inserted_atoms) == 1 and \
                list(removed_atoms.values())[0] == \
                list(inserted_atoms.values())[0]:
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
                 removed_atoms: dict,
                 inserted_atoms: dict,
                 changes_of_num_elements: dict,
                 charge: int,
                 initial_site_symmetry: str,
                 cutoff: float,
                 neighboring_sites: list,
                 annotation: str = None,
                 num_equiv_sites: int = None):
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
            removed_atoms (dict):
                Keys: Atom indices removed from the perfect supercell.
                      The index begins from 0.
                      For interstitials, set {}.
                Values: Fractional coordinates in the perfect positions.
            inserted_atoms (dict):
                Keys: Atom indices inserted to the supercell.
                      Indices are in  the defective supercell and begins from 0.
                      For vacancies, set [].
                Values: Fractional coordinates at the ideal positions.
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
        self.neighboring_sites = list(neighboring_sites)
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

        removed_atoms = {int(k): v for k, v in d["removed_atoms"].items()}
        inserted_atoms = {int(k): v for k, v in d["inserted_atoms"].items()}

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
            removed_atoms=removed_atoms,
            inserted_atoms=inserted_atoms,
            changes_of_num_elements=d["changes_of_num_elements"],
            charge=d["charge"],
            initial_site_symmetry=d["initial_site_symmetry"],
            cutoff=d["cutoff"],
            neighboring_sites=d["neighboring_sites"],
            annotation=d["annotation"],
            num_equiv_sites=d["num_equiv_sites"])

    def as_dict(self):
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
                  filename: str = None,
                  displacement_distance: float = 0.2,
                  symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                  angle_tolerance: float = ANGLE_TOL,
                  cutoff: float = CUTOFF_RADIUS,
                  calc_num_equiv_site: bool = True,
                  defect_name: str = None):
        """Construct the DefectEntry object from perfect and defective POSCARs.

        Note1: displacement_distance needs to be the same as the twice of max
               displacement distance for reducing the symmetry.
        Note2: Only unrelaxed but perturbed defective POSCAR structure is
               accepted.

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
        if filename is None:
            import shutil
            org_file = os.path.join(os.path.dirname(__file__),
                                    "default_defect_entry.yaml")
            shutil.copyfile(org_file, "defect_entry.yaml")
            filename = "defect_entry.yaml"

        abs_dir = os.path.split(os.path.abspath(filename))[0]

        with open(filename, "r") as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        displacement_distance = \
            yaml_data.get("displacement_distance", displacement_distance)

        # Perfect and perturbed defect structures.
        perfect_structure = Structure.from_file(
            os.path.join(abs_dir, yaml_data["perfect_structure"]))
        defect_structure = Structure.from_file(
            os.path.join(abs_dir, yaml_data["defect_structure"]))
        element_diff = element_diff_from_structures(
            defect_structure, perfect_structure)

        if not defect_name:
            _, defect_name = os.path.split(os.getcwd())
        name, charge, annotation = divide_dirname(defect_name)

        inserted_atom_indices = [i for i in range(defect_structure.num_sites)]
        removed_atoms = {}

        for i, p_site in enumerate(perfect_structure):
            for j in inserted_atom_indices:
                d_site = defect_structure[j]
                distance = p_site.distance(d_site)
                # check displacement_distance and species for comparison
                if distance < displacement_distance \
                        and p_site.specie == d_site.specie:
                    inserted_atom_indices.remove(j)
                    break
            else:
                removed_atoms[i] = list(p_site.frac_coords)

        # check the consistency of the removed and inserted atoms
        if len(perfect_structure) + len(inserted_atom_indices) \
                - len(removed_atoms) != len(defect_structure):
            raise StructureError(
                "Atoms are not properly mapped within the displacement.")

        # pristine_defect_structure is defect structure without displacement.
        pristine_defect_structure = deepcopy(perfect_structure)
        lattice = defect_structure.lattice
        for r in sorted(removed_atoms, reverse=True):
            pristine_defect_structure.pop(r)

        # If the inserted atom locates near an original atom within the
        # displacement_distance, it is assumed to be substituted.
        inserted_atoms = {}
        for i in sorted(inserted_atom_indices):
            inserted_atom = defect_structure[i]
            inserted_atoms[i] = list(inserted_atom.frac_coords)
            # substituted
            for removed_frac_coords in removed_atoms.values():
                if lattice.get_distance_and_image(
                        inserted_atom.frac_coords, removed_frac_coords)[0] \
                        < displacement_distance:
                    pristine_defect_structure.insert(i, inserted_atom.specie,
                                                     removed_frac_coords)
                    break
            # interstitial
            else:
                pristine_defect_structure.insert(i, inserted_atom.specie,
                                                 inserted_atom.frac_coords)

        defect_center_coords = cls.calc_defect_center(removed_atoms,
                                                      inserted_atoms,
                                                      defect_structure)
        neighboring_sites = []
        for i, site in enumerate(pristine_defect_structure):
            distance = \
                site.distance_and_image_from_frac_coords(defect_center_coords)[0]
            # Defect itself is not included to the neighboring sites
            if 1e-5 < distance < cutoff:
                neighboring_sites.append(i)

        sga = SpacegroupAnalyzer(structure=pristine_defect_structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance)

        initial_site_symmetry = sga.get_point_group_symbol()

        inserted_atom_coords = list(inserted_atoms.values())

        if calc_num_equiv_site:
            num_equiv_sites = \
                count_equivalent_clusters(perfect_structure,
                                          inserted_atom_coords,
                                          list(removed_atoms.keys()),
                                          displacement_distance)
        else:
            num_equiv_sites = None

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
        total_nions_in_perfect = \
            len(self.initial_structure) - len(self.inserted_atoms) \
            + len(self.removed_atoms)

        # initial atom mapping.
        mapping = list(range(total_nions_in_perfect))

        for o in sorted(self.removed_atoms.keys(), reverse=True):
            mapping.pop(o)

        for i in sorted(self.inserted_atoms.keys(), reverse=True):
            mapping.insert(i, None)

        return mapping

    def to_json_file(self, filename="defect_entry.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @staticmethod
    def calc_defect_center(removed_atoms: dict,
                           inserted_atoms: dict,
                           structure: Structure):
        """ Calculates arithmetic average to estimate center in frac coords."""
        defect_coords = \
            list(removed_atoms.values()) + list(inserted_atoms.values())
        return defect_center_from_coords(defect_coords, structure)

    @property
    def defect_center_coords(self):
        """ Return fractional coordinates of the defect center. """
        return self.calc_defect_center(self.removed_atoms, self.inserted_atoms,
                                       self.initial_structure)

    @property
    def anchor_atom_index(self):
        """ Returns an index of atom that is the farthest from the defect.

        This atom is assumed not to displace in the defective supercell, and
        so used for analyzing local structure around a defect.
        Note that only the first occurrence is returned when using argmax.
        docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html
        """
        # distance_set = \
        #     self.initial_structure.lattice.get_all_distances(
        #         self.defect_center, self.initial_structure.frac_coords)[0]

        return anchor_atom_index(self.initial_structure, self.defect_center_coords)


def anchor_atom_index(structure, center):
    """ Returns an index of atom that is the farthest from the defect.

    This atom is assumed not to displace in the defective supercell, and
    so used for analyzing local structure around a defect.
    Note that only the first occurrence is returned when using argmax.
    docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html
    """
    distance_set = \
        structure.lattice.get_all_distances(center, structure.frac_coords)[0]

    return np.argmax(distance_set)


def divide_dirname(dirname):
    """

    After the dirname is split by "_", there must be only one digit in the list,
    which is considered as a charge state, and Front (back) part corresponds to
    name (annotation).

    "Va_Mg1_2" -> name = "Va_Mg1", charge = 2, annotation = None
    "Va_O1_2_inward" -> name = "Va_O1", charge = 2, annotation = "inward"
    "Mg_i+Va_O1*2_2_coord1"
            -> name = "Mg_i+Va_O1*2", charge = 2, annotation = "coord1"
    """
    def is_digit(n):
        try:
            int(n)
            return True
        except ValueError:
            return False

    split_dirname = dirname.split("_")
    digit_positions = [x for x, y in enumerate(split_dirname) if is_digit(y)]

    if len(digit_positions) != 1:
        raise ValueError("The dirname {} is not valid".format(dirname))
    else:
        digit_pos = digit_positions[0]
        name = "_".join(split_dirname[:digit_pos])
        charge = int(split_dirname[digit_pos])
        annotation = "_".join(split_dirname[digit_pos + 1:])
        annotation = None if not annotation else annotation

    return name, charge, annotation


