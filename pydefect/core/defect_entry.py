# -*- coding: utf-8 -*-

import itertools
import json
import os
import ruamel.yaml as yaml
from pydefect.vasp_util.util import element_diff_from_poscar_files, \
    get_defect_charge_from_vasp

from pymatgen.core.structure import Structure

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def get_num_atoms_for_elements(structure):
    """
    Returns numbers of ions for elements in a structure.
    Example: Al1Mg31O32
        return: [1, 31, 32]
    """
    symbol_list = [site.specie.symbol for site in structure]

    return [len(tuple(a[1])) for a in itertools.groupby(symbol_list)]


class DefectEntry(MSONable):
    """

    This class object holds all the information related to initial setting of a
    single defect.
    """

    def __init__(self, name, initial_structure, perturbed_initial_structure,
                 removed_atoms, inserted_atoms, changes_of_num_elements, charge,
                 initial_symmetry, multiplicity, perturbed_sites):
        """
        Args:
            name (str):
                Name of a defect without charge. This is used when analyzing defect
                formation energy.
            initial_structure (Structure):
                Structure with a defect before the structure optimization.
            perturbed_initial_structure (Structure):
                Initial structure with perturbation of neighboring atoms.
            removed_atoms (dict):
                Keys: Atom indices removed from the perfect supercell.
                      The index begins from 0.
                      For interstitials, set {}.
                Values: DefectSupercell coordinates
            inserted_atoms (list):
                Atom indices inserted in the supercell after removing atoms.
                The index begins from 0.
                For vacancies, set [].
            changes_of_num_elements (dict):
                Keys: Element names
                Values: Change of the numbers of elements wrt perfect supercell.
            charge (int):
                Charge state of the defect. Charge is also added to the structure.
            initial_symmetry (str):
                Initial site symmetry such as D4h.
            multiplicity (int):
                Site multiplicity or number of equivalent sites in the given structure.
            perturbed_sites (list):
                Indices of the perturbed site for reducing the symmetry
        """
        self.name = name
        self.initial_structure = initial_structure
        self.perturbed_initial_structure = perturbed_initial_structure
        self.removed_atoms = removed_atoms
        self.inserted_atoms = inserted_atoms
        self.changes_of_num_elements = changes_of_num_elements
        self.charge = charge
        self.initial_symmetry = initial_symmetry
        self.multiplicity = multiplicity
        self.perturbed_sites = perturbed_sites

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        return self.as_dict() == other.as_dict()

    def __str__(self):
        outs = ["name: " + str(self.name),
                "initial_symmetry: " + str(self.initial_symmetry),
                "multiplicity: " + str(self.multiplicity),
                "removed_atoms: " + str(self.removed_atoms),
                "inserted_atoms: " + str(self.inserted_atoms),
                "changes_of_num_element: " + str(self.changes_of_num_elements),
                "charge: " + str(self.charge),
                "initial_structure: \n" + str(self.initial_structure),
                "perturbed_initial_structure: \n" +
                str(self.perturbed_initial_structure),
                "perturbed_sites: \n" + str(self.perturbed_sites)]
        return "\n".join(outs)

    @classmethod
    def from_dict(cls, d):
        """ Constructs a DefectEntry class object from a dictionary.
        """
        # The keys need to be converted to integers.
        removed_atoms = {int(k): v for k, v in d["removed_atoms"].items()}
        changes_of_num_elements = {k: int(v) for k, v in d["changes_of_num_elements"].items()}

        return cls(d["name"], d["initial_structure"],
                   d["perturbed_initial_structure"], removed_atoms,
                   d["inserted_atoms"], changes_of_num_elements, d["charge"],
                   d["initial_symmetry"], d["multiplicity"],
                   d["perturbed_sites"])

    # TODO: add perturbed_structure
    # @classmethod
    # def from_yaml(cls, filename, tolerance=0.1):
    #     """
    #     An example of the yaml file.
    #         name: 2Va_O1 + Mg_i_2
    #         initial_structure: POSCAR
    #         perturbed_initial_structure: POSCAR
    #         perfect_structure: ../../defects/perfect/POSCAR
    #         charge: 2 (optional, otherwise calc from INCAR and POTCAR)
    #         tolerance: 0.2 (optional)
    #     """

        # abs_dir = os.path.split(os.path.abspath(filename))[0]

        # with open(filename, "r") as yaml_file:
        #     yaml_data = yaml.safe_load(yaml_file)

        # if "tolerance" in yaml_data.keys():
        #     tolerance = yaml_data["tolerance"]

        # element_diff = element_diff_from_poscar_files(
        #     os.path.join(abs_dir, yaml_data["initial_structure"]),
        #     os.path.join(abs_dir, yaml_data["perfect_structure"]))

        # defect_structure = Structure.from_file(
        #     os.path.join(abs_dir, yaml_data["initial_structure"]))
        # perfect_structure = Structure.from_file(
        #     os.path.join(abs_dir, yaml_data["perfect_structure"]))

        # # set name
        # if "name" in yaml_data.keys():
        #     name = yaml_data["name"]
        # else:
        #     _, name = os.path.split(os.getcwd())
        #     print("name", name, "is set from the directory name.")

        # # set charge state
        # if "charge" in yaml_data.keys():
        #     charge = yaml_data["charge"]
        # else:
        #     nions = get_num_atoms_for_elements(defect_structure)
        #     charge = get_defect_charge_from_vasp(nions=nions)
        #     print("charge", charge, "is set from vasp input files.")

        # inserted_atoms = [i for i in range(defect_structure.num_sites)]
        # removed_atoms = {}

        # for i, p_site in enumerate(perfect_structure):
        #     for j in inserted_atoms:
        #         d_site = defect_structure[j]
        #         distance = p_site.distance(d_site)
        #         # check distance and species for comparison
        #         if distance < tolerance and p_site.specie == d_site.specie:
        #             inserted_atoms.remove(j)
        #             break
        #     # *else* block is active if *for j* loop is not broken.
        #     # else is not recommended in effective python as it's confusing.
        #     else:
        #         removed_atoms[i] = list(p_site.frac_coords)

        # # check the consistency of the removed and inserted atoms
        # if not (sum([i for i in element_diff.values() if i > 0])
        #         == len(inserted_atoms)
        #         and sum([-i for i in element_diff.values() if i < 0])
        #         == len(removed_atoms)):
        #     raise ImproperInputStructureError(
        #         "Atoms in two structures are not mapped in the tolerance.")

        # return cls(name, defect_structure, removed_atoms, inserted_atoms,
        #            element_diff, charge, None, None)

    @classmethod
    def load_json(cls, filename="defect_entry.json"):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def atom_mapping_to_perfect(self):
        """
        Returns a list of atom mapping from defect structure to perfect.
        Example of Mg32O32 supercell:
            When 33th atom, namely first O, is removed,
                mapping = [0, 1, 2, .., 31, 33, 34, .., 62]
                len(mapping) = 63

        """
        total_nions = (sum(get_num_atoms_for_elements(self.initial_structure))
                       - len(self.inserted_atoms)
                       + len(self.removed_atoms))

        # initial atom mapping.
        mapping = list(range(total_nions))

        for o in sorted(self.removed_atoms.keys(), reverse=True):
            mapping.pop(o)

        for i in sorted(self.inserted_atoms, reverse=True):
            mapping.insert(i, None)

        return mapping

    def as_dict(self):
        """
        Dict representation of DefectInput class object.
        """
        d = {"name": self.name,
             "initial_structure": self.initial_structure,
             "perturbed_initial_structure": self.perturbed_initial_structure,
             "removed_atoms": self.removed_atoms,
             "inserted_atoms": self.inserted_atoms,
             "changes_of_num_elements": self.changes_of_num_elements,
             "charge": self.charge,
             "initial_symmetry": self.initial_symmetry,
             "multiplicity": self.multiplicity,
             "perturbed_sites": self.perturbed_sites}
        return d

    def to_json_file(self, filename="defect_entry.json"):
        """
        Writes a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    # TODO: remove bugs below
    # def anchor_atom_index(self):
    #     """
    #     Returns an index of atom that is the farthest from the defect.
    #     This atom is assumed not to displace during the structure
    #     optimization, and so used for analyzing local defect structure.
    #     """
        # radius = max(self.initial_structure.lattice.abc) * 2
        # num_sites = len(self.initial_structure.sites)
        # shortest_distances = np.full(num_sites, radius, dtype=float)

        # distance_set = self.initial_structure.get_sites_in_sphere(
        #     self._defect_coords, radius, include_index=True)

        # for d in distance_set:
        #     atom_index = d[2]
        #     if d[1] < shortest_distances[atom_index]:
        #         shortest_distances[atom_index] = d[1]

        # farthest_atom_index = np.argmax(shortest_distances)
        # farthest_dist = shortest_distances[farthest_atom_index]

        # return farthest_atom_index, farthest_dist


class ImproperInputStructureError(Exception):
    pass


