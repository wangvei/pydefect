#!/usr/bin/env python

import json
from monty.json import MontyEncoder
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


class Defect:
    """
    This class object holds some properties related to a defect.
    Args:
        structure (Structure): pmg Structure/IStructure class object
        removed_atom_index (array of int): Atom index removed in the perfect
                                           supercell.
                                           For interstitial, set to None.
        inserted_atom_index (array of int): Atom index inserted in the supercell after
                                   removing an atom. For vacancy, set to None.
        defect_coords (Nx3 array): coordinates of defect position
        defect_center (Nx3 array): center of a defect.
                                   If len(defect_coords) == 1,
                                   same as defect_coords[0].
        in_name" (str): Inserted element name. "Va" is inserted for vacancies.
        out_name" (str): Removed site name. "in", where n is an integer,
                         is inserted for interstitials. E.g., "i1".
        charge (int): Charge state of the defect
    """
    def __init__(self, removed_atom_index, inserted_atom_index, defect_coords,
                 in_name, out_name, charge):
        self.removed_atom_index = removed_atom_index
        self.inserted_atom_index = inserted_atom_index
        self.defect_coords = defect_coords
        self.in_name = in_name
        self.out_name = out_name
        self.charge = charge

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a dictionary.
        """
        return cls(d["removed_atom_index"], d["inserted_atom_index"],
                   d["defect_coords"], d["in_name"], d["out_name"], d["charge"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a DefectSetting class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        """
        d = {"removed_atom_index": self.removed_atom_index,
             "inserted_atom_index": self.inserted_atom_index,
             "defect_coords": self.defect_coords,
             "in_name": self.in_name,
             "out_name": self.out_name,
             "charge": self.charge}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)


class IrreducibleSite:
    """
    This class object holds properties related to irreducible atom set.
    Note1: atomic indices need to be sorted. Thus, they can be written in one
           sequence. E.g., 17..32
    Note2: first_index atom is assumed to represent the irreducible atoms.

    Args:
        irreducible_name (str): element name with irreducible index (e.g., Mg1)
        element (str): element name (e.g., Mg)
        first_index (int): first index of irreducible_name.
        last_index (int): last index of irreducible_name.
        repr_coords (array): representative coordinates, namely the position
                            of first_index

    TODO1: Add the site symmetry information.
    """
    def __init__(self, irreducible_name, element, first_index, last_index,
                 repr_coords):
        self.irreducible_name = irreducible_name
        self.element = element
        self.first_index = first_index
        self.last_index = last_index
        self.repr_coords = repr_coords

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        return self.__dict__ == other.__dict__

    def as_dict(self):
        d = {"irreducible_name": self.irreducible_name,
             "element" : self.element,
             "first_index" : self.first_index,
             "last_index" : self.last_index,
             "repr_coords" : self.repr_coords}
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["irreducible_name"], d["element"], d["first_index"],
                   d["last_index"], d["repr_coords"])

    @property
    def natoms(self):
        """
        Returns number of atoms in a given (super)cell.
        """
        return self.last_index - self.first_index + 1
