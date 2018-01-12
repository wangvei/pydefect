#!/usr/bin/env python
import os
import shutil
import numpy as np
import warnings
import argparse
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
#import atom
from copy import deepcopy
import itertools as it
import sys
import re
import ruamel.yaml as yaml
#import parse_poscar as ppos
#import pydefect.input_generator.defect_input.DefectSetting as DefectSetting 

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pydefect.yaml")

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
        cutoff (float): Radius of sphere in which atoms are perturbed [A]. 
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
        random_vector = random_vector(normed_random_3D_vector(), distance)
        site = i[2]
        sites.append(site)
        structure.translate_sites(site, random_vector, frac_coords=False)

    return {"structure": structure, "sites": sites}


def potcar_dir():
    """    
    Return the name of POTCAR file directory.
    SETTINGS_FILE needs to be defined in the same module.
    """    
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except:
        raise IOError('.pydefect.yaml cannot be opened.')

    for k, v in d.items():
        if k == "DEFAULT_POTCAR":
            potcar_dir = v

    if not potcar_dir:
        raise ValueError('DEFAULT_POTCAR is not set in .pydefect.yaml')

    return potcar_dir

def make_POTCAR(dirname, elements, default_potcar_dir):
    """    
    Write POTCAR with a sequence of given elements names at *dirname*.
    So far, only default POTCAR files are supported.    
    """    
    with open(dirname + '/POTCAR', 'w') as potcar:
        for e in elements:
            potcar_file_name = default_potcar_dir + "/POTCAR_" + e
            with open(potcar_file_name) as pot:
                potcar.write(pot.read())

def _get_int_from_string(x):
    """ return int number only """
    return int(''.join(i for i in x if i.isdigit() or i == '.'))


class DefectName():
    """    
    Extract defect information from a defect directly name,
    e.g., "Va_Mg1_-2" --> in_name = "Va", out_name = "Mg1", charge = -2
    """    
    def __init__(self, defect_name):
        try:
            d = defect_name.split("_")
            self.in_name = d[0]
            self.out_name = d[1]
            self.charge = int(d[2])
        except:
            raise ValueError("Defect {} is inappropriate.", defect_name)

        if not re.match(r'^[a-xA-Z]+[1-9]+$', self.out_name):
            print(re.match(r'^[a-xA-Z]+[1-9]+$', self.out_name))
            raise ValueError("Defect {} is inappropriate.", defect_name)


class VaspInputMaker():
    """
    Construct a set of vasp input files.
    POTCAR files are fetched from ~/.pydefect.yaml.
    
    Args:
        defect_setting: DefectSetting class object
        defect_name (str): defect name defined in PyDefect, e.g., "Va_Mg1_2"
        poscar (str): DPOSCAR name
        incar (str): INCAR name
        kpoints (str): KPOINTS name
    
    """

    def __init__(self, defect_setting, defect_name, poscar="DPOSCAR", 
                 incar="INCAR", kpoints="KPONTS"):

        if os.path.exists(defect_name):
            print("{} alreadly exists, so nothing is done.".format(defect_name))     

        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints
        self.defect_name = DefectName(defect_name)
        # Need deepcopy for structure
        self.defect_structure = deepcopy(defect_setting.structure)

        irrep_element_names = [irrep_element.irrepname 
                            for irrep_element in defect_setting.irrep_elements]

        print(irrep_element_names)
        # Vacancy
        if self.defect_name.in_name == "Va" and \
                            self.defect_name.out_name in irrep_element_names:
            for irrep_element in defect_setting.irrep_elements:
                if self.defect_name.out_name == irrep_element.irrepname:
                    vacancy_atomic_index = irrep_element.first_index
                    self.defect_position = irrep_element.repr_coord

            if "vacancy_atomic_index" not in locals():
                raise ValueError("{} does not exist in defect setting.".format(
                                                    self.defect_name.out_name))
            self.defect_structure.remove_sites([vacancy_atomic_index - 1])
            self.defect_type = "vacancy"

        # Interstitial
        elif Element.is_valid_symbol(self.defect_name.in_name) and \
                             re.match(r'^i[0-9]+$', self.defect_name.out_name):

            interstitial_index = \
                            _get_int_from_string(self.defect_name.out_name) - 1
            try:
                coords = defect_setting.interstitial_coords[interstitial_index]
            except:
                raise ValueError(
             "Interstitial index {} does not exist".format(interstitial_index))
            for irrep_element in defect_setting.irrep_elements:
                candidate_atomic_indices = []
                if self.defect_name.in_name == irrep_element.element:
                    candidate_atomic_indices.append(irrep_element.first_index)
            # Intrinsic atom is inserted just before the same elements.
            if len(candidate_atomic_indices) > 0:
                  interstitial_atomic_index = min(candidate_atomic_indices) - 1
            # Extrinsic atom is appended.
            else:
                interstitial_atomic_index = 0
                
            self.defect_structure.insert(interstitial_atomic_index, 
                                             self.defect_name.in_name, coords)
            self.defect_type = "interstitial"

        # Substitutional + antisite
        elif Element.is_valid_symbol(self.defect_name.in_name) and \
                              self.defect_name.out_name in irrep_element_names:





            for irrep_element in defect_setting.irrep_elements:
                candidate_atomic_indices = []
                if self.defect_name.in_name == irrep_element.element:
                    candidate_atomic_indices.append(irrep_element.first_index)
            # Intrinsic atom is inserted just before the same elements.
            if len(candidate_atomic_indices) > 0:
                substituted_atomic_index = min(candidate_atomic_indices) - 1
            # Extrinsic atom is appended.
            else:
                interstitial_atomic_index = 0


        # randomly perturb neighboring atoms.
    

#    def construct_input_files(defect_setting, in_name, out_name, charge, incar, poscar, 
#                              potcar, kpoints, runshell):
#        """    
#        """    
#    
#    # TODO: this will be moved 
#    #    if os.path.exists(dirname):
#    #        print("{} alreadly exist, so nothing is done.", dirname)     
#    #    else:
#    
#        # POSCAR
#        new_structure = make_structure(defect_setting, in_name, out_name)
#        
#    
#        # POTCAR
#        elements = 
#        
#    
#        # INCAR
#        shutil.copyfile(incar, dirname + "/INCAR")
#        # NELECT
#    
#        # KPOINTS & runshell.sh
#        os.makedirs(dirname)
#        shutil.copyfile(kpoints, dirname + "/KPOINTS")
#        shutil.copyfile(runshell, dirname + "/" + runshell)
#
#
#    
class VaspInputSetMaker():
   pass
##
##    def __init__(self, defect_setting, incar="INCAR", structure="DPOSCAR",
##                 kpoints="KPONTS", run="run.sh", cutoff=4.0):
#
#
#        
#
##    make_POTCAR(dirname,        
#    
#
#
##    def _get_elements_from_names(self, host_atoms):
##        """
##        host_atoms[name(str)] = [int(rep), [charge(int)]]
##        Return element names with charges, which are one-size-fits-all.
##        E.g. host_atoms[Sn1] = [0, 1, 2],host_atoms[Sn2] = [2, 3, 4]
##             -> elements[Sn] = [0, 1, 2, 3, 4]
##        """
##        elements = {}
##
##        for name in names:
##            # Remove number. E.g., Mg1 -> Mg
##            element = name[0:-1]
##            if element not in elements:
##                elements[element] = set()
##
##            for charge in host_atoms[atom][1]:
##                elements[name].add(charge)
##
##        return elements
##
##
