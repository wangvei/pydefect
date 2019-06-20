# -*- coding: utf-8 -*-
from copy import deepcopy
import json
import numpy as np
import os
import ruamel.yaml as yaml

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.error_classes import StructureError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import count_equivalent_clusters
from pydefect.vasp_util.util import element_diff_from_poscar_files, \
    get_defect_charge_from_vasp
from pydefect.util.structure_tools import defect_center_from_coords
from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class DefectEntry(MSONable):
    """ Holds information related to the initial setting of a single defect. """

    def __init__(self,
                 name: str,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 removed_atoms: dict,
                 inserted_atoms: dict,
                 changes_of_num_elements: dict,
                 charge: int,
                 initial_site_symmetry: str,
                 neighboring_sites: list,
                 annotation: str = None,
                 num_equiv_sites: int = None):
        """
        Args:
            name (str):
                Name of a defect without charge.
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
            neighboring_sites (list):
                Indices of the perturbed site for reducing the symmetry
            annotation (str):

            num_equiv_sites (int):
                Number of equivalent sites in the given structure.
        """
        self.name = name
        self.initial_structure = initial_structure
        self.perturbed_initial_structure = perturbed_initial_structure
        self.removed_atoms = deepcopy(removed_atoms)
        self.inserted_atoms = deepcopy(inserted_atoms)
        self.changes_of_num_elements = deepcopy(changes_of_num_elements)
        self.charge = charge
        self.initial_site_symmetry = initial_site_symmetry
        self.neighboring_sites = list(neighboring_sites)
        self.num_equiv_sites = num_equiv_sites

    def __str__(self):
        outs = ["perturbed_initial_structure: \n" +
                str(self.perturbed_initial_structure),
                "name: " + str(self.name),
                "initial_site_symmetry: " + str(self.initial_site_symmetry),
                "removed_atoms: " + str(self.removed_atoms),
                "inserted_atoms: " + str(self.inserted_atoms),
                "changes_of_num_element: " + str(self.changes_of_num_elements),
                "charge: " + str(self.charge),
                "neighboring_sites: " + str(self.neighboring_sites),
                "num_equiv_sites: " + str(self.num_equiv_sites)]
        return "\n".join(outs)

    @classmethod
    def from_dict(cls, d: dict):
        # The keys need to be converted to integers.
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
            initial_structure=initial_structure,
            perturbed_initial_structure=perturbed_initial_structure,
            removed_atoms=removed_atoms,
            inserted_atoms=inserted_atoms,
            changes_of_num_elements=d["changes_of_num_elements"],
            charge=d["charge"],
            initial_site_symmetry=d["initial_site_symmetry"],
            neighboring_sites=d["neighboring_sites"],
            num_equiv_sites=d["num_equiv_sites"])

    @classmethod
    def from_yaml(cls,
                  filename: str,
                  displacement_distance: float = 0.2,
                  symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                  angle_tolerance: float = ANGLE_TOL,
                  perturbed_criterion: float = 0.001):
        """Construct the DefectEntry object from perfect and defective POSCARs.

        Note1: displacement_distance needs to be the same as the twice of max
               displacement distance for reducing the symmetry.
        Note2: Only unrelaxed but perturbed defective POSCAR structure is
               accepted.

        filename (str): yaml filename.
        displacement_distance (float): Tolerance to determine the same atom
        angle_tolerance (float):

        An example of the yaml file.
            name (optional): 2Va_O1 + Mg_i_2
            defect_structure: POSCAR
            perfect_structure: ../../defects/perfect/POSCAR
            charge (optional, if none, calculated from INCAR and POTCAR): 2
            displacement_distance (optional): 0.2
        """
        abs_dir = os.path.split(os.path.abspath(filename))[0]

        with open(filename, "r") as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        displacement_distance = \
            yaml_data.get("displacement_distance", displacement_distance)

        element_diff = element_diff_from_poscar_files(
            os.path.join(abs_dir, yaml_data["defect_structure"]),
            os.path.join(abs_dir, yaml_data["perfect_structure"]))

        # Perfect and perturbed defect structures.
        perfect_structure = Structure.from_file(
            os.path.join(abs_dir, yaml_data["perfect_structure"]))
        defect_structure = Structure.from_file(
            os.path.join(abs_dir, yaml_data["defect_structure"]))

        # set defect name
        if "name" in yaml_data.keys():
            name = yaml_data["name"]
        else:
            _, dirname = os.path.split(os.getcwd())
            name = "_".join(dirname.split("_")[:2])
            logger.warning("name: {} is set from the directory name.".
                           format(name))

        # set charge state
        if "charge" in yaml_data.keys():
            charge = yaml_data["charge"]
        else:
            charge = get_defect_charge_from_vasp(structure=defect_structure)
            logger.warning("charge {} is set from vasp input files.".
                           format(charge))

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
        neighboring_sites = []
        for i, site in enumerate(defect_structure):
            pristine_defect_site = pristine_defect_structure[i]
            distance = site.distance(pristine_defect_site)
            if perturbed_criterion < distance < displacement_distance:
                neighboring_sites.append(i)

        sga = SpacegroupAnalyzer(structure=pristine_defect_structure,
                                 symprec=symprec,
                                 angle_tolerance=angle_tolerance)

        initial_site_symmetry = sga.get_point_group_symbol()

        inserted_atom_coords = list(inserted_atoms.values())

        num_equiv_sites = \
            count_equivalent_clusters(perfect_structure,
                                      inserted_atom_coords,
                                      list(removed_atoms.keys()),
                                      displacement_distance)

        return cls(name=name,
                   initial_structure=pristine_defect_structure,
                   perturbed_initial_structure=defect_structure,
                   removed_atoms=removed_atoms,
                   inserted_atoms=inserted_atoms,
                   changes_of_num_elements=element_diff,
                   charge=charge,
                   initial_site_symmetry=initial_site_symmetry,
                   neighboring_sites=neighboring_sites,
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

    @property
    def defect_center(self):
        """ Return fractional coordinates of the defect center.

        Calculates the arithmetic average of the defect positions.
        """
        total_defect_coords = (list(self.removed_atoms.values()) +
                               list(self.inserted_atoms.values()))
        return defect_center_from_coords(total_defect_coords,
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

        return anchor_atom_index(self.initial_structure, self.defect_center)


def anchor_atom_index(structure, center):
    """ Returns an index of atom that is the farthest from the defect.

    This atom is assumed not to displace in the defective supercell, and
    so used for analyzing local structure around a defect.
    Note that only the first occurrence is returned when using argmax.
    docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html
    """
    distance_set = structure.lattice.get_all_distances(center,
                                                       structure.frac_coords)[0]

    return np.argmax(distance_set)

