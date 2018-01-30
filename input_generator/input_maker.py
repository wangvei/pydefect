#!/usr/bin/env python
import os
import shutil
import numpy as np
import warnings
import argparse
import json
import itertools as it
import sys
import re
import ruamel.yaml as yaml
from copy import deepcopy
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.core.periodic_table import Element
from pydefect.input_generator.defect_in import DefectSetting 

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def normed_random_3D_vector():
    """
    Generates a random 3D unit vector with a uniform spherical distribution.
    stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """
    phi = np.random.uniform(0, np.pi*2)
    costheta = np.random.uniform(-1, 1)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

def random_vector(normed_vector, distance):
    """
    Return a vector scaled by distance * x, where 0<x<1.

    Args:
        normed_vector (3x1 array): Normed 3D vector.
        distance (float): distance
    """
    return normed_vector * distance * np.random.random()

def perturb_around_a_point(structure, center, cutoff, distance):
    """
    Randomly perturb atoms around an input point in a structure.

    Args:
        center (3x1 array): Fractional coordinates of a central position.
        cutoff (float): Radius of a sphere in which atoms are perturbed [A]. 
        distance (float): Max distance for the perturbation [A].
    """
    if  type(center) == list and len(center) == 3:
        cartesian_coords = structure.lattice.get_cartesian_coords(center)
        neighbors = structure.get_sites_in_sphere(
                                  cartesian_coords, cutoff, include_index=True)
    else:
        raise ValueError
    sites = []
    # Since translate_sites accepts only one vector, we need to iterate this.
    for i in neighbors:
        vector = random_vector(normed_random_3D_vector(), distance)
        site = i[2]
        sites.append(site)
        structure.translate_sites(site, vector, frac_coords=False)

    return structure, sites
#    return {"structure": structure, "sites": sites}

def _get_int_from_string(x):
    """ 
    Return integer number from a string.
    """
    return int(''.join(i for i in x if i.isdigit() or i == '.'))

def _parse_defect_name(defect_name):
    """ 
    Divide defect name to three part.
    E.g., "Va_Mg1_0" --> in_name="Va", out_name="Mg1", charge=0
    """
    try:
        d = defect_name.split("_")
        in_name = d[0]
        out_name = d[1]
        charge = int(d[2])
    except:
        raise ValueError("Defect {} is improper.".format(defect_name))

    if not re.match(r'^[a-xA-Z]+[1-9]+$',out_name):
        raise ValueError("Defect {} is improper.".format(defect_name))
    return (in_name, out_name, charge)

def extended_range(i):
    """
    Extension of range method especially for negative input value.
    E.g., extended_range(3) = [0, 1, 2, 3]
          extended_range(-3) = [-3, -2, -1, 0]
    """
    if not type(i) == int:
        raise AttributeError
    if i >= 0: return range(i + 1)
    else : return range(i, 1)

def _print_already_exist(dirname):
    print("{:>10} alreadly exists, so nothing is done.".format(dirname))
    
def _print_is_constructed(dirname):
    print("{:>10} is constructed.".format(dirname))


class InputMaker():
    """
    Structure information is stored in defect_setting.

    Args:
        defect_name (str): defect name defined in PyDefect, e.g., "Va_Mg1_2"
        defect_setting: DefectSetting class object
    """
    def __init__(self, defect_name, defect_setting):

        if os.path.exists(defect_name):
            self.is_directory = True
        else:
            self.is_directory = False
        self.defect_name = defect_name
        self.defect_setting = defect_setting
        self.in_name, self.out_name, self.charge = \
                                            _parse_defect_name(defect_name)
        # deepcopy is needed for structure
        self.defect_structure = deepcopy(defect_setting.structure)
        # e.g., irreducible_names = ["Mg1", "O1"]
        self.irreducible_names = [i.irreducible_name 
                            for i in self.defect_setting.irreducible_sites]

    def analyze_name(self):
        if re.match(r'^i[0-9]+$', self.out_name):
            interstitial_index = _get_int_from_string(self.out_name)
            try:
                defect_coords = \
                self.defect_setting.interstitial_coords[interstitial_index - 1]
            except:
                raise ValueError(
                 "Interstitial # {} is not defined".format(interstitial_index))
        elif self.out_name in self.irreducible_names:
            # There may be multiple candidates for inserted element.
            for irreducible_site in self.defect_setting.irreducible_sites:
                if self.out_name == irreducible_site.irreducible_name:
                    removed_atomic_index = irreducible_site.first_index
                    defect_coords = irreducible_site.repr_coords
            self.defect_structure.remove_sites([removed_atomic_index - 1])
        else:
            raise ValueError("{} in {} is improper.".\
                                       format(self.out_name, self.defect_name))

        if self.in_name == "Va":
            self.defect_index = removed_atomic_index
            self.defect_coords = defect_coords
        elif Element.is_valid_symbol(self.in_name):
            # There may be multiple candidates for inserted element.
            # E.g., Mg exists in Mg1 and Mg2.
            # *in_name* element is inserted just before the same elements, 
            # othewise to the 1st index.
            candidate_atomic_indices = []
            # check all the irreducible_sites.
            if self.in_name in self.defect_structure.symbol_set:
                atomic_index = \
                   min(self.defect_structure.indices_from_symbol(self.in_name))
            else:
                atomic_index = 0
            self.defect_structure.insert(atomic_index, self.in_name, 
                                                                 defect_coords)
            self.defect_index = atomic_index + 1
            self.defect_coords = self.defect_structure.\
                                             frac_coords[atomic_index].tolist()
        else:
            raise ValueError("{} in {} is improper.".\
                                       format(self.in_name, self.defect_name))

    def make_directory_json(self):
        """
        Needs to be modified by subclasses depending on the code.
        """
#        if self.is_directory == False:
        os.makedirs(self.defect_name)
        # write a defect position to defect.json file.
        with open(self.defect_name + "/defect.json", 'w') as fw:
            json.dump({"defect_index": self.defect_index,
                       "defect_coords": self.defect_coords,
                       "in_name": self.in_name,
                       "out_name": self.out_name,
                       "charge": self.charge}, fw, indent=2)
#        else:
#            print("{} exists. is not defined".format(interstitial_index))

    def make_perturbed_defect_structure(self):
        # perturb neighboring atoms randomly.
       self.perturbed_defect_structure, self.perturbed_sites = \
                           perturb_around_a_point(self.defect_structure, 
                                                  self.defect_coords, 
                                                  self.defect_setting.cutoff, 
                                                  self.defect_setting.displace)
