#!/usr/bin/env python

from abc import ABCMeta
import json

import numpy as np
from monty.json import MontyEncoder
from monty.serialization import loadfn
from pydefect.core.DFT_results import SupercellDftResults

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
    Example: Al1Mg63O6
        return: [1, 63, 64]

    """
    nions = [int(i)
             for i in defect_structure.to(fmt="poscar").split("\n")[6].split()]
    return nions


class Supercell(metaclass=ABCMeta):
    """
    """
    def __init__(self):
        self._dft_results = None

    @property
    def final_structure(self):
        return self._dft_results.final_structure

    @property
    def total_energy(self):
        return self._dft_results.total_energy

    @property
    def eigenvalues(self):
        return self._dft_results.eigenvalues

    @property
    def electrostatic_potential(self):
        return self._dft_results.electrostatic_potential

    @property
    def ewald_param(self):
        return self._dft_results.ewald_param

    def set_dft_results(self, directory_path, contcar_name="/CONTCAR",
                        outcar_name="/OUTCAR", vasprun_name="/vasprun.xml"):
        print("aaa")
        self._dft_results = \
            SupercellDftResults.from_vasp_files(
                directory_path, contcar_name, outcar_name, vasprun_name)


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
        self._initial_structure = initial_structure
        self._removed_atom_index = removed_atom_index
        self._inserted_atom_index = inserted_atom_index
        self._defect_coords = defect_coords
        self._in_name = in_name
        self._out_name = out_name
        self._charge = charge

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        # Expansion of _dft_results is necessary,
        _dft_results = []
        for i in d["irreducible_sites"]:
            irreducible_sites.append(IrreducibleSite.from_dict(i))

        return cls(d["initial_structure"], d["removed_atom_index"],
                   d["inserted_atom_index"], d["defect_coords"], d["in_name"],
                   d["out_name"], d["charge"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a Defect class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def initial_structure(self):
        return self._initial_structure

    @property
    def removed_atom_index(self):
        return self._removed_atom_index

    @property
    def inserted_atom_index(self):
        return self._inserted_atom_index

    @property
    def defect_coords(self):
        return self._defect_coords

    @property
    def in_name(self):
        return self._in_name

    @property
    def out_name(self):
        return self._out_name

    @property
    def charge(self):
        return self._charge

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        """
        d = {"initial_structure": self._initial_structure,
             "removed_atom_index": self._removed_atom_index,
             "inserted_atom_index": self._inserted_atom_index,
             "defect_coords": self._defect_coords,
             "in_name": self._in_name,
             "out_name": self._out_name,
             "charge": self._charge}
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
        return [np.mean(i) for i in np.array(self._defect_coords).transpose()]

    def atom_mapping_to_perfect(self):
        """
        Returns a list of atom mapping in defect structure to atoms in perfect.
        Example Mg32O32 supercell:
            Assuming that 33th atom (first O) is removed,
            mapping = [0, 1, 2, .., 31, 33, 34, .., 62]
            len(mapping) = 63
        """
        total_nions = sum(get_nions(self._initial_structure))
        mapping = [i for i in range(total_nions)]
        in_index = self._inserted_atom_index
        out_index = self._removed_atom_index

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

        radius = max(self._initial_structure.lattice.abc) * 2
        num_sites = len(self._initial_structure.sites)
        shortest_distances = np.full(num_sites, radius, dtype=float)

        distance_set = self._initial_structure.get_sites_in_sphere(
            self._defect_coords, radius, include_index=True)

        for d in distance_set:
            atom_index = d[2]
            if d[1] < shortest_distances[atom_index]:
                shortest_distances[atom_index] = d[1]

        farthest_atom_index = np.argmax(shortest_distances)
        farthest_dist = shortest_distances[farthest_atom_index]

        return farthest_atom_index, farthest_dist
