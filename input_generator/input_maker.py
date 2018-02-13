#!/usr/bin/env python
import os
import re
from abc import ABCMeta, abstractmethod
from copy import deepcopy

import numpy as np
from defect import Defect
from pymatgen.core.periodic_table import Element

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def normed_random_3d_vector():
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


def perturb_around_a_point(structure, center, cutoff, distance):
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
        vector = random_vector(normed_random_3d_vector(), distance)
        site = i[2]
        sites.append(site)
        structure.translate_sites(site, vector, frac_coords=False)

    return structure, sites


def _get_int_from_string(x):
    """ 
    Returns integer numbers from a string.

    Args:
        x (str): a string
    """
    return int(''.join(i for i in x if i.isdigit() or i == '.'))


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
    except:
        raise ValueError("Defect {} is improper.".format(defect_name))

    if not re.match(r'^[a-xA-Z]+[1-9]+$', out_name):
        raise ValueError("Defect {} is improper.".format(defect_name))
    return (in_name, out_name, charge)


def extended_range(i):
    """
    Extends range method used for positive/negative input value.
    The input value is included even for a positive number.
    E.g., extended_range(3) = [0, 1, 2, 3]
          extended_range(-3) = [-3, -2, -1, 0]

    Args:
        i (int): an integer
    """
    if not type(i) == int:
        raise AttributeError
    if i >= 0:
        return range(i + 1)
    else:
        return range(i, 1)


def print_already_exist(name):
    """
    Show the following message.
    Args:
        name (str): a string
    """
    print("{:>10} already exists, so nothing is done.".format(name))


def print_is_constructed(name):
    """
    Show the following message.
    Args:
        name (str): a string
    """
    print("{:>10} is being constructed.".format(name))


class DefectInputMaker:
    """
    Construct information for a defect. DefectSetting class object is required.

    Args:
        defect_name (str): defect name in PyDefect manner, e.g., "Va_Mg2_-2".
        defect_setting (DefectSetting) : DefectSetting class object.

    Parameters in use:
        in_name" (str): Inserted element name. "Va" is inserted for vacancies.
        out_name" (str): Removed site name. "in", where n is an integer,
                         is inserted for interstitials. E.g., "i1".
        charge (int): Charge state of the defect
    """

    def __init__(self, defect_name, defect_setting):

        self.is_directory = os.path.exists(defect_name)
        self.defect_name = defect_name
        self.defect_setting = defect_setting
        self.in_name, self.out_name, self.charge = \
            parse_defect_name(defect_name)
        # deepcopy is required for modifying the original structure.
        self.defect_structure = deepcopy(defect_setting.structure)

        # -------------------- analyze out_name --------------------------------
        # interstitial
        if re.match(r'^i[0-9]+$', self.out_name):
            interstitial_index = _get_int_from_string(self.out_name)
            try:
                self.defect_coords = \
                    self.defect_setting.interstitial_coords[
                        interstitial_index - 1]
                self.removed_atom_index = None
            except ValueError:
                print("{}th interstitial is not defined".format(
                    interstitial_index))
        else:
            # if out_name exists in irreducible_sites
            for i in self.defect_setting.irreducible_sites:
                if self.out_name == i.irreducible_name:
                    self.removed_atom_index = i.first_index - 1
                    self.defect_coords = i.repr_coords
            try:
                self.defect_structure.remove_sites([self.removed_atom_index])
            except ValueError:
                print("{} in defect {} is improper.".format(
                    self.out_name, self.defect_name))

        # -------------------- analyze in_name ---------------------------------
        # This method needs to be run after finishing analyze_out_name because
        # self.removed_atom_index and self.defect_coords are required for
        # determining self.defect_index and interstitial coordinates.
        if self.in_name == "Va":
            self.inserted_atom_index = None
        #            self.defect_index = self.removed_atom_index
        elif Element.is_valid_symbol(self.in_name):
            # There may be multiple irreducible sites for inserted element,
            # e.g., Mg1 and Mg2, element of in_name is inserted to just before
            # the same elements, otherwise to the 1st index.
            if self.in_name in self.defect_structure.symbol_set:
                self.inserted_atom_index = \
                    min(self.defect_structure.indices_from_symbol(self.in_name))
            else:
                inserted_atom_index = 0
            self.defect_structure.insert(
                self.inserted_atom_index, self.in_name, self.defect_coords)
        #    self.defect_coords = self.defect_structure.frac_coords[
        #                                       inserted_atom_index].tolist()
        else:
            raise ValueError("{} in {} is improper.". \
                             format(self.in_name, self.defect_name))

    def make_json(self, filename):
        """
        Needs to be modified by subclasses depending on the code.
        """
        d = Defect(self.removed_atom_index, self.inserted_atom_index,
                   self.defect_coords, self.in_name, self.out_name, self.charge)
        d.to_json_file(filename)


class DefectInputSetMaker(metaclass=ABCMeta):
    """
    Abstract class that must be subclassed by a particular first-principles
    code implementation.
    Construct a set of defect class object

    Args:
        defect_setting (DefectSetting): DefectSetting class object.
        specific_defects (str/list): It specifies a particular defect(s).
            Specify a type of defects.
                "Va" --> A set of all the vacancies.
                "i" --> A set of all the interstitials.
                "as" --> A set of all the antisites.
                "Va_O" --> A set of all the oxygen vacancies
                "Va_O1" --> A set of oxygen vacancies at O1 site
                "Va_O1_2" --> Only Va_O1_2 (
                "Mg_O" --> A set of all the Mg-on-O antisite pairs.
                "Mg_O1" --> A set of Mg-on-O1 antisite pairs.
                "Mg_O1_0" --> A set of Mg-on-O1 antisite pairs.
            List type is also allowed.
                ["Va", "i"] --> A set of all the vacancies and interstitials.
    Parameters in use:
    """

    def __init__(self, defect_setting, specific_defects=None):

        self.defect_setting = defect_setting

        if len(specific_defects.split("_")) == 3:
            self.defect_set = specific_defects
        else:
            self.elements = self.defect_setting.structure.symbol_set
            # Construct a set of defect class objects.
            self.defect_set = \
                self._vacancy_name_setter() + \
                self._interstitial_name_setter() + \
                self._substitutional_name_setter()

            for i in self.defect_setting.included:
                self.defect_set.append(i)
            for e in self.defect_setting.excluded:
                if e in self.defect_set:
                    self.defect_set.remove(e)
                else:
                    print("{} does not exist.".format(e))

    def _vacancy_name_setter(self, element=None):
        name_set = []
        for i in self.defect_setting.irreducible_sites:

            oxidation_state = self.defect_setting.oxidation_states[i.element]
            for o in extended_range(-oxidation_state):
                defect_name = "Va_" + i.irreducible_name + "_" + str(o)
                name_set.append(defect_name)
        return name_set

    def _interstitial_name_setter(self):
        name_set = []
        for e in self.elements:
            oxidation_state = self.defect_setting.oxidation_states[e]
            for j in range(len(self.defect_setting.interstitial_coords)):
                for o in extended_range(oxidation_state):
                    defect_name = e + "_i" + str(j + 1) + "_" + str(o)
                    name_set.append(defect_name)
        return name_set

    def _substitutional_name_setter(self):
        name_set = []
        for a in self.defect_setting.antisite_configs + \
                 self.defect_setting.dopant_configs:
            in_element, out_element = a
            oxidation_state_diff = \
                self.defect_setting.oxidation_states[in_element] \
                - self.defect_setting.oxidation_states[out_element]
            for i in self.defect_setting.irreducible_sites:
                if out_element == i.element:
                    for o in extended_range(oxidation_state_diff):
                        defect_name = in_element + "_" + \
                                      i.irreducible_name + "_" + str(o)
                        name_set.append(defect_name)
        return name_set

    @abstractmethod
    def make_perfect_input_files(self):
        pass
