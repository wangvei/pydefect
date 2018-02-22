#!/usr/bin/env python

import json

import numpy as np
from monty.json import MontyEncoder
from monty.serialization import loadfn
from copy import deepcopy

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def get_nions(defect_structure):
    """
    Return numbers of ions for elements in defect_structure.
    """
    nions = [int(i)
             for i in defect_structure.to(fmt="poscar").split("\n")[6].split()]
    return nions


class Supercell:
    """
    DFT result of supercell without any defect.
    """
    def __init__(self):
        self._structure = None
        self._energy = None
        self._eigenvalue = None
        self._electrostatic_potential = None
        self._ewald_param = None

    def set_calc_results(self, directory_path, poscar_name = "/CONTCAR",
            outcar_name = "/OUTCAR", vasprun_name = "/vasprun.xml"):
        """
        Args:
            directory_path (str): path of directory.
            poscar_name (str): Name of POSCAR file. Defaults to CONTCAR.
            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
            vasprun_name (str): Name of vasprun.xml file.
                                Defaults to vasprun.xml.
        """
        self._structure = contcar.structure
        path_contcar = directory_path + contcar_name
        path_outcar = directory_path + outcar_name
        path_vasprun = directory_path + vasprun_name
        contcar = Poscar.from_file(path_contcar)
        outcar = Outcar(path_outcar)
        vasprun = Vasprun(path_vasprun)
        self._structure = contcar.structure
        self._energy = outcar.final_energy
        self._eigenvalue = vasprun.eigenvalues
        self._electrostatic_potential = outcar.electrostatic_potential
        self._ewald_param = None

    @property
    def structure(self):
        return self._structure

    @property
    def energy(self):
        return self._energy

    @property
    def eigenvalue(self):
        return self._eigenvalue

    @property
    def electrostatic_potential(self):
        return self._electrostatic_potential

    @property
    def ewald_param(self):
        return self._ewald_param

    @ewald_param.setter
    def ewald_param_setter(self, ewald_param):
        self._ewald_param = ewald_param


class Perfect(Supercell):
    def __init__(self):
        super().__init__()


class Defect(Supercell):
    """
    This class object holds some properties related to a defect.
    Args:
        initial_structure (Structure): pmg Structure/IStructure class object.
            Defect structure before structure optimization
        removed_atom_index (array of int):
            Atom index removed in the perfect supercell, starting from 0.
            For interstitial, set to None.
        inserted_atom_index (array of int):
            Atom index inserted in the supercell after removing an atom.
            For vacancy, set to None.
        defect_coords (Nx3 array): coordinates of defect position
        in_name" (str): Inserted element name. "Va" is inserted for vacancies.
        out_name" (str): Removed site name. "in", where n is an integer,
                         is inserted for interstitials. E.g., "i1".
        charge (int): Charge state of the defect
    """
    def __init__(self, initial_structure, removed_atom_index,
                 inserted_atom_index, defect_coords, in_name, out_name, charge):
        super().__init__()
        self.initial_structure = initial_structure
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
        return cls(d["initial_structure"], d["removed_atom_index"],
                   d["inserted_atom_index"], d["defect_coords"], d["in_name"],
                   d["out_name"], d["charge"])

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
        d = {"initial_structure": self.initial_structure,
             "removed_atom_index": self.removed_atom_index,
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

    def defect_center(self):
        """
        Returns coords of defect center by calculating the average coords.
        If len(defect_coords) == 1, same as defect_coords[0].
        """
        return [np.mean(i) for i in np.array(self.defect_coords).transpose()]

    def atom_mapping_to_perfect(self):
        """
        Returns a list of atom mapping in defect structure to atoms in perfect.
        """
        total_nions = sum(get_nions(self.initial_structure))
        mapping = [i for i in range(total_nions)]
        in_index = self.inserted_atom_index
        out_index = self.removed_atom_index

        if in_index:
            mapping[in_index] = None
            for i in range(in_index + 1, total_nions):
                mapping[i] += -1

        if out_index is not None:
            for i in range(out_index, total_nions):
                mapping[i] += 1

        return mapping

    def anchor_atom_index(self):
        """
        Returns an index of atom that is the farthest from the defect.
        This atom is assumed not to displace during the structure
        optimization, and so used for analyzing defect structure.
        """

        radius = max(self.initial_structure.lattice.abc) * 2
        num_sites = len(self.initial_structure.sites)
        shortest_distances = np.full(num_sites, radius, dtype=float)

        distance_set = self.initial_structure.get_sites_in_sphere(
            self.defect_coords, radius, include_index=True)

        for d in distance_set:
            atom_index = d[2]
            if d[1] < shortest_distances[atom_index]:
                shortest_distances[atom_index] = d[1]

        farthest_atom_index = np.argmax(shortest_distances)
        farthest_dist = shortest_distances[farthest_atom_index]

        return farthest_atom_index, farthest_dist


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

