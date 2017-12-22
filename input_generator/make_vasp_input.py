#!/usr/bin/env python
import os
import shutil
import numpy as np
import warnings
import argparse
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
#from pymatgen.core.periodic_table import Element, Specie, get_el_sp
#import atom
import itertools as it
import sys
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

#class MakeParticularInput():

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pydefect.yaml")

def random_three_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0, np.pi*2)
    costheta = np.random.uniform(-1, 1)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return (np.array([x, y, z]))

def perturb_around_defect(structure, defect_position, cutoff, distance):
    """
    Performs a random perturbation of the sites in a structure to break
    symmetries.

    Args:
        defect_position: Fractional coordinates (3x1 array)
        cutoff (float): Radius of sphere. 
        distance (float): Distance in angstroms by which to perturb each site.
    """

    if  type(defect_position) == list and len(defect_position) == 3:
        cartesian_defect_position = \
                        structure.lattice.get_cartesian_coords(defect_position)
        neighbors = structure.get_sites_in_sphere(
                         cartesian_defect_position, cutoff, include_index=True)
    else:
        raise ValueError

    def get_rand_vec():
        return random_three_vector() * distance * np.random.random()

    sites = []
    # Since translate_sites accept only one translation vector, we need to
    # iterate this.
    for i in neighbors:
        print(i[2])
        site = i[2]
        sites.append(site)
        structure.translate_sites(site, get_rand_vec(), frac_coords=False)

    return {"structure": structure, "sites": sites}

def _laad_potcar_dir():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except:
        assert IOError('.pydefect.yaml at the home directory cannot be opened.'):

    for k, v in d.items():
        if k == "DEFAULT_POTCAR":
            potcar_dir = v

    if not potcar_dir:
        assert ValueError('DEFAULT_POTCAR is not set in .pydefect.yaml'):

    return potcar_dir

def make_POTCAR(dirname, elements, default_potcar_dir):
    """    
    Construct POTCAR file from a given default_potcar_path and elements names.

    """    
    with open(dirname + '/POTCAR', 'w') as potcar:
        for s in elements:
            potcar_file_name = default_potcar_path + "/POTCAR_" + s
            shutil.copyfileobj(open(potcar_file_name), potcar)


class VaspInputMaker():
    """
    Constructs a set of vasp input files.
    
    Args:
        defect_setting: DefectSetting class object
        structure: pymatgen structure object
        incar (str): INCAR file name
        kpoints (str): KPOINTS file name
    
    + POTCAR files are fetched from ~/.pydefect.yaml
    """


    def __init__(self, defect_setting, structure="DPOSCAR", incar="INCAR", 
                 kpoints="KPONTS", cutoff=4.0):

#class VaspInputSetMaker():
#
#    def __init__(self, defect_setting, incar="INCAR", structure="DPOSCAR",
#                 kpoints="KPONTS", run="run.sh", cutoff=4.0):

def make_structure(defect_setting, in_name, out_name):

    structure = deepcopy(defect_setting.structure)        

        # Vacancy
        if in_name == "Va":
            site = 
            structure.remove_sites(
    
        # Interstitial
        if re.match("i[0-9]*" out_name):
            
        # Substitutional + antisite


    # randomly displace neighboring atoms.

    return {"structure": structure, "defect_position": defect_position}


def construct_input_files(defect_setting, in_name, out_name, charge, incar, poscar, 
                          potcar, kpoints, runshell):
    """    
    """    

    dirname = str(in_name) + "_" + str(out_name) + "_" + str(charge)

# TODO: this will be moved 
#    if os.path.exists(dirname):
#        print("{} alreadly exist, so nothing is done.", dirname)     
#    else:

    # POSCAR
    new_structure = make_structure(defect_setting, in_name, out_name)
    

    # POTCAR
    elements = 
    

    # INCAR
    shutil.copyfile(incar, dirname + "/INCAR")
    # NELECT

    # KPOINTS & runshell.sh
    os.makedirs(dirname)
    shutil.copyfile(kpoints, dirname + "/KPOINTS")
    shutil.copyfile(runshell, dirname + "/" + runshell)

    

        

#    make_POTCAR(dirname,        
    


#    def _get_elements_from_names(self, host_atoms):
#        """
#        host_atoms[name(str)] = [int(rep), [charge(int)]]
#        Return element names with charges, which are one-size-fits-all.
#        E.g. host_atoms[Sn1] = [0, 1, 2],host_atoms[Sn2] = [2, 3, 4]
#             -> elements[Sn] = [0, 1, 2, 3, 4]
#        """
#        elements = {}
#
#        for name in names:
#            # Remove number. E.g., Mg1 -> Mg
#            element = name[0:-1]
#            if element not in elements:
#                elements[element] = set()
#
#            for charge in host_atoms[atom][1]:
#                elements[name].add(charge)
#
#        return elements
#
#
    def _get_num(self, x):
        """ return number only """
        return float(''.join(i for i in x if i.isdigit() or i == '.'))
