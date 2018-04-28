# -*- coding: utf-8 -*-

import itertools
import json
import os
import ruamel.yaml as yaml

from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.io.vasp import Incar
from pymatgen.io.vasp.inputs import Potcar

from monty.json import MontyEncoder
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


def element_diff_from_poscar_files(poscar1, poscar2):
    """
    Returns a dict of change of numbers of elements from poscar2 to poscar1
    For defect calculations, poscar2 should be "perfect".
    """
    c1 = Composition(
        Structure.from_file(poscar1).composition, allow_negative=True)
    c2 = Composition(
        Structure.from_file(poscar2).composition, allow_negative=True)
    c_diff = c1 - c2

    return {str(e): int(c_diff[e]) for e in c_diff}


def get_num_electrons_from_potcar(potcar, nions, charge=0):
    """
    Returns the number of electrons from POTCAR, number of ions, and charge
    state.
    """
    p = Potcar.from_file(potcar)
    # check only the number of ions written in potcar and nions.
    if not len(p) == len(nions):
        raise ValueError("Size of elements in POTCAR file is different")

    return sum([v.nelectrons * nions[i] for i, v in enumerate(p)]) - charge


def get_defect_charge_from_vasp(nions, potcar="POTCAR", incar="INCAR"):
    """
    Returns the defect charge by comparing nion, number of electrons in POTCAR,
    and NELECT in INCAR.
    """
    num_elect_neutral = get_num_electrons_from_potcar(potcar, nions)
    num_elect_incar = Incar.from_file(incar)["NELECT"]

    # charge is minus of difference of the electrons
    return int(num_elect_neutral - num_elect_incar)


class DefectEntry:
    """
    This class object holds all the information related to initial setting of a
    single defect.
    Args:
        name (str):
            Name of a defect without charge. This is used when analyzing defect
            formation energy.
        initial_structure (Structure):
            Structure with a defect before the structure optimization.
        removed_atoms (dict):
            Keys: Atom indices removed from the perfect supercell.
                  The index begins from 0.
                  For interstitials, set {}.
            Values: DefectSupercell coordinates
        inserted_atoms (list):
            Atom indices inserted in the supercell after removing atoms.
            The index begins from 0.
            For vacancies, set [].
        element_diff (dict):
            Keys: Element names
            Values: Change of the numbers of elements wrt perfect supercell.
        charge (int):
            Charge state of the defect
    """
    def __init__(self, name, initial_structure, removed_atoms, inserted_atoms,
                 element_diff, charge):
        self._name = name
        self._initial_structure = initial_structure
        self._removed_atoms = removed_atoms
        self._inserted_atoms = inserted_atoms
        self._element_diff = element_diff
        self._charge = charge

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        return self.as_dict() == other.as_dict()

    def __str__(self):
        outs = ["name:" + str(self._name),
                "removed_atoms:" + str(self._removed_atoms),
                "inserted_atoms:" + str(self._inserted_atoms),
                "element_diff:" + str(self._element_diff),
                "charge:" + str(self._charge),
                "structure: \n" + str(self._initial_structure)]
        return "\n".join(outs)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectEntry class object from a dictionary.
        """
        # The keys need to be converted to integers.
        removed_atoms = {int(k): v for k, v in d["removed_atoms"].items()}
        element_diff = {k: int(v) for k, v in d["element_diff"].items()}

        return cls(d["name"], d["initial_structure"], removed_atoms,
                   d["inserted_atoms"], element_diff, d["charge"])

    @classmethod
    def from_yaml(cls, filename, tolerance=0.1):
        """
        An example of the yaml file.
            name: 2Va_O1+Mg_i_2
            initial_structure: POSCAR
            perfect_structure: ../../defects/perfect/POSCAR
            charge: 2 (optional, otherwise calc from INCAR and POTCAR)
            tolerance: 0.2 (optional)
        """

        abs_dir = os.path.split(os.path.abspath(filename))[0]

        with open(filename, "r") as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        if "tolerance" in yaml_data.keys():
            tolerance = yaml_data["tolerance"]

        element_diff = element_diff_from_poscar_files(
            os.path.join(abs_dir, yaml_data["initial_structure"]),
            os.path.join(abs_dir, yaml_data["perfect_structure"]))

        defect_structure = Structure.from_file(
            os.path.join(abs_dir, yaml_data["initial_structure"]))
        perfect_structure = Structure.from_file(
            os.path.join(abs_dir, yaml_data["perfect_structure"]))

        # set name
        if "name" in yaml_data.keys():
            name = yaml_data["name"]
        else:
            _, name = os.path.split(os.getcwd())
            print("name", name, "is set from the directory name.")

        # set charge state
        if "charge" in yaml_data.keys():
            charge = yaml_data["charge"]
        else:
            nions = get_num_atoms_for_elements(defect_structure)
            charge = get_defect_charge_from_vasp(nions=nions)
            print("charge", charge, "is set from the initial files.")

        inserted_atoms = [i for i in range(defect_structure.num_sites)]
        removed_atoms = {}

        for i, p_site in enumerate(perfect_structure):
            for j in inserted_atoms:
                d_site = defect_structure[j]
                distance = p_site.distance(d_site)
                # check distance and species for comparison
                if distance < tolerance and p_site.specie == d_site.specie:
                    inserted_atoms.remove(j)
                    break
            # *else* block is active if *for j* loop is not broken.
            # else is not recommended in effective python as it's confusing.
            else:
                removed_atoms[i] = list(p_site.frac_coords)

        # check the consistency of the removed and inserted atoms
        if not (sum([i for i in element_diff.values() if i > 0])
                == len(inserted_atoms)
                and sum([-i for i in element_diff.values() if i < 0])
                == len(removed_atoms)):
            raise ImproperInputStructureError(
                "Atoms in two structures are not mapped in the tolerance.")

        return cls(name, defect_structure, removed_atoms, inserted_atoms,
                   element_diff, charge)

    @classmethod
    def load_json(cls, filename="defect_entry.json"):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def name(self):
        return self._name

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
    def element_diff(self):
        return self._element_diff

    @property
    def charge(self):
        return self._charge

    @property
    def atom_mapping_to_perfect(self):
        """
        Returns a list of atom mapping from defect structure to perfect.
        Example of Mg32O32 supercell:
            When 33th atom, namely first O, is removed,
                mapping = [0, 1, 2, .., 31, 33, 34, .., 62]
                len(mapping) = 63

        """
        total_nions = (sum(get_num_atoms_for_elements(self._initial_structure))
                       - len(self._inserted_atoms)
                       + len(self._removed_atoms))

        # initial atom mapping.
        mapping = list(range(total_nions))

        for o in sorted(self._removed_atoms.keys(), reverse=True):
            mapping.pop(o)

        for i in sorted(self._inserted_atoms, reverse=True):
            mapping.insert(i, None)

        return mapping

    def as_dict(self):
        """
        Dict representation of DefectInput class object.
        """
        d = {"name": self._name,
             "initial_structure": self._initial_structure,
             "removed_atoms": self._removed_atoms,
             "inserted_atoms": self._inserted_atoms,
             "element_diff": self._element_diff,
             "charge": self._charge}
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


class ImproperInputStructureError(Exception):
    pass


