#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import numpy as np

from monty.json import MontyEncoder
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def get_nions(structure):
    """
    Return numbers of ions for elements in a structure.
    Example: Al1Mg63O6
        return: [1, 63, 64]

    """
    return [int(i) for i in structure.to(fmt="poscar").split("\n")[6].split()]


class DefectEntry:
    """
    This class object holds all the information related to initial setting of a
    single defect.
    Args:
        initial_structure (Structure):
            Structure with a defect before the structure optimization.
        removed_atoms (dict):
            Keys: Atom indices removed from the perfect supercell.
                  The index begins from 0.
                  For interstitials, set {}.
            Values: DefectSupercell coordinates
        inserted_atoms (dict):
            Keys: Atom indices inserted in the supercell after removing atoms.
                  The index begins from 0.
                  For vacancies, set {}.
            Values: DefectSupercell coordinates
        changes_of_num_elements (dict):
            Keys: Element names
            Values: Change of the numbers of elements wrt perfect supercell.
        charge (int):
            Charge state of the defect
        in_name (str): Optional.
            Inserted element name.
            "Va" is inserted for vacancies.
            When in_name is not well defined, e.g, complex defects, set False.
        out_name (str): Optional.
            Removed atomic site name.
            "in", (n: integer), is inserted for interstitials. E.g., "i1".
            When out_name is not well defined, set False.
    """
    def __init__(self, initial_structure, removed_atoms, inserted_atoms,
                 changes_of_num_elements, charge, in_name=False,
                 out_name=False):
        self._initial_structure = initial_structure
        self._removed_atoms = removed_atoms
        self._inserted_atoms = inserted_atoms
        self._changes_of_num_elements = changes_of_num_elements
        self._charge = charge
        self._in_name = in_name
        self._out_name = out_name

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        return self.as_dict() == other.as_dict()

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectEntry class object from a dictionary.
        """
        # The keys need to be converted to integers.
        removed_atoms = {int(k): v for k, v in d["removed_atoms"].items()}
        inserted_atoms = {int(k): v for k, v in d["inserted_atoms"].items()}
        changes_of_num_elements = \
            {k: int(v) for k, v in d["changes_of_num_elements"].items()}

        return cls(d["initial_structure"], removed_atoms, inserted_atoms,
                   changes_of_num_elements, d["charge"], d["in_name"],
                   d["out_name"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def initial_structure(self):
        return self._initial_structure

    @property
    def removed_atoms(self):
        return self._removed_atoms

    @property
    def inserted_atoms(self):
        return self._inserted_atoms

    @property
    def changes_of_num_elements(self):
        return self._changes_of_num_elements

    @property
    def charge(self):
        return self._charge

    @property
    def in_name(self):
        return self._in_name

    @property
    def out_name(self):
        return self._out_name

    @property
    def defect_center(self):
        """
        Returns coordinates of defect center by calculating the averaged
        coordinates. If len(defect_coords) == 1, returns defect_coords[0].
        """
        defect_coords = (list(self._removed_atoms.values()) +
                         list(self._inserted_atoms.values()))
        return [np.mean(i) for i in np.array(defect_coords).transpose()]

    @property
    def atom_mapping_to_perfect(self):
        """
        Returns a list of atom mapping from defect structure to perfect.
        Example of Mg32O32 supercell:
            When 33th atom, namely first O, is removed,
                mapping = [0, 1, 2, .., 31, 33, 34, .., 62]
                len(mapping) = 63

        """
        total_nions = (sum(get_nions(self._initial_structure))
                       - len(self._inserted_atoms)
                       + len(self._removed_atoms))

        # initial atom mapping.
        mapping = list(range(total_nions))

        for o in sorted(self._removed_atoms.keys(), reverse=True):
            mapping.pop(o)

        for i in sorted(self._inserted_atoms.keys(), reverse=True):
            mapping.insert(i, None)

        return mapping

    def as_dict(self):
        """
        Dict representation of DefectInput class object.
        """
        d = {"initial_structure": self._initial_structure,
             "removed_atoms": self._removed_atoms,
             "inserted_atoms": self._inserted_atoms,
             "changes_of_num_elements": self._changes_of_num_elements,
             "charge": self._charge,
             "in_name": self._in_name,
             "out_name": self._out_name}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
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
        # radius = max(self._initial_structure.lattice.abc) * 2
        # num_sites = len(self._initial_structure.sites)
        # shortest_distances = np.full(num_sites, radius, dtype=float)

        # distance_set = self._initial_structure.get_sites_in_sphere(
        #     self._defect_coords, radius, include_index=True)

        # for d in distance_set:
        #     atom_index = d[2]
        #     if d[1] < shortest_distances[atom_index]:
        #         shortest_distances[atom_index] = d[1]

        # farthest_atom_index = np.argmax(shortest_distances)
        # farthest_dist = shortest_distances[farthest_atom_index]

        # return farthest_atom_index, farthest_dist
