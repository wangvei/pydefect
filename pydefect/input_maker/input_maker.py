# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from copy import deepcopy
import numpy as np
import re

from pymatgen.core.periodic_table import Element

from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def normalized_random_3d_vector():
    """
    Generates a random 3d unit vector with a uniform spherical distribution.
    stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """
    phi = np.random.uniform(0, np.pi * 2)
    costheta = np.random.uniform(-1, 1)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def random_vector(normed_vector, distance):
    """
    Returns a vector scaled by distance * x, where 0<x<1.

    Args:
        normed_vector (3x1 array): Normed 3d vector.
        distance (float): distance
    """
    return normed_vector * distance * np.random.random()


def perturb_neighbors(structure, center, cutoff, distance):
    """
    Randomly perturbs atoms around an input point in a structure.

    Args:
        structure (Structure): pmg Structure/IStructure class object
        center (3x1 array): Fractional coordinates of a central position.
        cutoff (float): Radius of a sphere in which atoms are perturbed [A].
        distance (float): Max distance for the perturbation [A].
    """
    if type(center) == list and len(center) == 3:
        cartesian_coords = structure.lattice.get_cartesian_coords(center)
        neighbors = structure.get_sites_in_sphere(
            cartesian_coords, cutoff, include_index=True)
    else:
        raise ValueError

    sites = []
    # Since translate_sites accepts only one vector, we need to iterate this.
    for i in neighbors:
        vector = random_vector(normalized_random_3d_vector(), distance)
        site = i[2]
        sites.append(site)
        structure.translate_sites(site, vector, frac_coords=False)

    return structure, sites


def get_int_from_string(x):
    """ 
    Returns joined integer number from a string.

    Args:
        x (str): a string
    """
    return int(''.join(i for i in x if i.isdigit()))


def parse_defect_name(defect_name):
    """ 
    Divides defect name by "_" to three part.
    E.g., "Va_Mg1_0" --> in_name="Va", out_name="Mg1", charge=0

    Args:
        defect_name (str): defect name defined in PyDefect, e.g., "Va_Mg1_2"
    """
    try:
        d = defect_name.split("_")
        in_name = d[0]
        out_name = d[1]
        charge = int(d[2])
    except ValueError:
        print("DefectSupercell {} is improper.".format(defect_name))

    if not re.match(r'^[a-xA-Z]+[1-9]+$', out_name):
        raise ValueError("DefectSupercell {} is improper.".format(defect_name))
    return in_name, out_name, charge


def print_already_exist(name):
    """
    Show the following message.
    Args:
        name (str): a string
    """
    print("{:>10} already exists, so nothing is done.".format(name))


def print_is_being_constructed(name):
    """
    Show the following message.
    Args:
        name (str): a string
    """
    print("{:>10} is being constructed.".format(name))


def filter_name_set(name_set, filtering_words):
    """
     Args:
        name_set (list): A set of defect names.
        filtering_words (list):

    When the following type names are given, constructs a set of defects.
        "Va"    --> A set of all the vacancies.
        "_i"     --> A set of all the interstitials.
        "Va_O"  --> A set of all the oxygen vacancies
        "Va_O1" --> A set of oxygen vacancies at O1 site
        "Mg_O"  --> A set of all the Mg-on-O antisite pairs.
        "Mg_O1" --> A set of Mg-on-O1 antisite pairs.

    When complete defect_name is given, constructs a particular defect.
        e.g., "Va_O1_2",  "Mg_O1_0"
    """
    filtered_name_set = []

    for p in filtering_words:
        pattern = r"" + re.escape(p)
        for d in name_set:
            if re.search(pattern, d):
                filtered_name_set.append(d)

    return list(set(filtered_name_set))


class DefectMaker:
    """
    Constructs a single DefectEntry class object from a given defect_name.

    Args:
        defect_name (str):
            defect name in PyDefect manner, e.g., "Va_Mg2_-2".
        structure (Structure):
            Structure class object for the *perfect* supercell.
        irreducible_sites ([IrreducibleSite])
        interstitial_coords (Nx3 list): coordinates of interstitial sites,
                                         e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ..]

    Parameters in use:
        in_name" (str): Inserted element name. "Va" is used for vacancies.
        out_name" (str): Removed site name. "in", where n is an integer,
                         is used for interstitials. E.g., "i1".
        charge (int): DefectSupercell charge state
        removed_atom_index (int): removed atom index from the perfect structure.
    """
    def __init__(self, defect_name, structure, irreducible_sites,
                 interstitial_coords):

        # deepcopy is required for modifying the original structure.
        defect_structure = deepcopy(structure)
        in_name, out_name, charge = parse_defect_name(defect_name)
        changes_of_num_elements = {}
        # -------------------- analyze out_name --------------------------------
        removed_atoms = {}
        # interstitial
        if re.match(r'^i[0-9]+$', out_name):
            interstitial_index = get_int_from_string(out_name)
            try:
                defect_coords = interstitial_coords[interstitial_index - 1]
            except ValueError:
                print("#{} interstitial not defined".format(interstitial_index))
        else:
            for i in irreducible_sites:
                if out_name == i.irreducible_name:
                    changes_of_num_elements[i.element] = -1
                    removed_index = i.first_index - 1
                    defect_coords = i.repr_coords
                    removed_atoms[removed_index] = i.repr_coords
                    break
            try:
                defect_structure.remove_sites([removed_index])
            except ValueError:
                print("{} in {} is improper.".format(out_name, defect_name))

        # -------------------- analyze in_name ---------------------------------
        inserted_atoms = []

        # This block needs to be run after finishing analyze_out_name because
        # defect coordinates is needed when inserting an in_name atom.
        if in_name == "Va":
            pass
        elif Element.is_valid_symbol(in_name):
            # There may be multiple irreducible sites for inserted element,
            # e.g., Mg1 and Mg2, element of in_name is inserted to just before
            # the same elements, otherwise to the 1st index.
            changes_of_num_elements[in_name] = 1
            if in_name in defect_structure.symbol_set:
                inserted_index = \
                    min(defect_structure.indices_from_symbol(in_name))
                inserted_atoms.append(inserted_index)
            else:
                inserted_index = 0
                inserted_atoms = [0]
            defect_structure.insert(inserted_index, in_name, defect_coords)
        else:
            raise ValueError("{} in {} is improper.".format(out_name,
                                                            defect_name))
        self.defect = DefectEntry(defect_structure, removed_atoms,
                                  inserted_atoms, changes_of_num_elements,
                                  charge)


class DefectInputSetMaker(metaclass=ABCMeta):
    """
    Abstract class that is subclassed by a particular first-principles code
    implementation. Constructs a set of DefectSupercell class objects based on
    given oxidation states.

    Args:
        defect_initial_setting (DefectInitialSetting):
            DefectInitialSetting class object.
        filtering_words (list):
            Specify a type of defects.
        particular_defects (list):
            It specifies (a) particular defect(s).

    Parameters in use:
        in_pattern (str): pattern for screening in_name
        out_pattern (str): pattern for screening out_name
    """

    def __init__(self, defect_initial_setting, filtering_words=None,
                 particular_defects=None):

        self._defect_initial_setting = defect_initial_setting
        self._filtering_words = filtering_words
        self._particular_defects = particular_defects

        if particular_defects:
            self._defect_name_set = particular_defects
        else:
            defect_all_name_set = defect_initial_setting.make_defect_name_set()
            if filtering_words:
                self._defect_name_set = \
                    filter_name_set(defect_all_name_set, filtering_words)
            else:
                defect_all_name_set.append("perfect")
                self._defect_name_set = defect_all_name_set

    def make_input(self):
        # Construct defect input files.
        for d in self._defect_name_set:
            if d == "perfect":
                self._make_perfect_input()
            else:
                self._make_defect_input(d)

    @abstractmethod
    def _make_perfect_input(self):
        pass

    @abstractmethod
    def _make_defect_input(self, defect_name):
        pass
