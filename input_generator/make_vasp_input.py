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
        cartesian_defect_position = structure.lattice.get_cartesian_coords(defect_position)
        neighbors = structure.get_sites_in_sphere(cartesian_defect_position, cutoff, include_index=True)
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


class VaspInputMaker():

    def __init__(self, defect_setting, incar="INCAR", structure="DPOSCAR",
                 kpoints="KPONTS", run="run.sh", cutoff=4.0):
        """

        + POTCAR files are fetched from ~/.pydefect.yaml

        """
        # TODO1: check if the input files exist.

#        self.make_dir_three_input_files("perfect", incar, kpoints, run)
#        default_potcar_path = _laad_potcar_dir()
#
#        if type(structure) = "str":
#            structure = Structure.from_file(poscar)
#        elif type(structure) = "Structure":
#            pass
#        else:
#            raise TypeError 
#
#        random_disp = {"cutoff":self.cutoff, "" 
#
#        # Make perfect directory
#        if os.path.isdir("perfect"):  
#            print "     perfect alreadly exist, so is not made."
#        else:
#            os.mkdir("perfect")
#            self.make_dir_three_input_files("perfect", incar, kpoints, run)
#            shutil.copyfile(poscar, "perfect/POSCAR")
#            make_POTCAR("perfect", atomic_sites.keys(), default_potcar_path)
        
#       # construct vacancies
#        # Vacancies. Charges are sign reversed.
#        make_vacancies(atomic_sites, r_)
#        for name, key in atomic_sites:
#            for charge in range(int(key["charge"])):
##                element = ''.join([i for i in k if not i.isdigit()])
#                dirname="Va_"+ name + "_" + charge
#                os.mkdir(dirname)
#                make_dir_poscar("Va", name, -charge, "DPOSCAR", r_param,
#                                removed_atom=key["rep"])
#
#                
#        # Interstitials. 
#        for element in host_elements: 
#            for interstitial_index, coord in enumerate(interstitial_site):
#                for charge in host_elements[element]:
#        
#                    make_dir_poscar(element, "i" + str(interstitial_index + 1), charge,
#                                    opts.poscar, r_param, added_coord=coord,
#                                    added_atom_symbol=element)
#        # Antisites. 
#        for a in anti_site: 
#            # E.g., Mg1_O2 -> insert = Mg1, remove = O2
#            added_atom, removed_atom = [x for x in a.split("_")]
#            added_atom_charge = host_elements[added_atom]
#            removed_atom_charge = [-i for i in host_elements[removed_atom[0:-1]]]
#        
#            # Construct the possible charge from charge_b and charge_c
#            charges =list(set([sum(i) for i in 
#                             list(it.product(added_atom_charge, removed_atom_charge))]))
#        
#            for charge in charges:
#                make_dir_poscar(added_atom, removed_atom, charge, opts.poscar, r_param, 
#                                removed_atom=host_atoms[removed_atom][0], 
#                                added_atom_symbol=added_atom)
#        
#        # Dopants. Both substitutional and intersitials are considered.
#        # dopants[Mg] = [0, 1, 2]
#        for dopant in dopants:
#        #for dopant, dopant_charges in enumerate(dopants):
#            # substituted
#            for host_removed in host_atoms:
#                host_atom_charge = [-i for i in host_atoms[host_removed][1]]
#                charges = list(set([sum(i) for i in list(
#                                       it.product(dopants[dopant], host_atom_charge))]))
#                for charge in charges:
#        
#                    make_dir_poscar(dopant, host_removed, charge, opts.poscar, r_param,
#                                    removed_atom=host_atoms[host_removed][0], 
#                                    added_atom_symbol=dopant)
#            # Interstitials. 
#            for interstitial_index, coord in enumerate(interstitial_site):
#                for charge in dopants[dopant]:
#        
#                    make_dir_poscar(dopant, "i" + str(interstitial_index + 1), charge, 
#                                    opts.poscar, r_param, added_coord=coord, 
#
#    def make_vacancies():
#
#
#    def make_dir_poscar(inserted_element, removed_name, charge, poscar, r_param,
#                        removed_atom=False, added_coord=False, 
#                        added_atom_symbol=False):
#    """    
#    """    
#
#        name = inserted_element + "_" + removed_name + "_" + str(charge)
#    
#        if os.path.isdir(name): 
#            print "%12s alreadly exist, so is not made."% name
#        else: 
#            sp.call(["mkdir", name])
#            ppos.make_point_defects(poscar, removed_atom, added_coord,
#                                    added_atom_symbol, name + '/POSCAR',
#                                    r_param, header=name)
#
#

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


def make_POTCAR(dirname, elements, default_potcar_dir):
    """    
    Construct POTCAR file from a given default_potcar_path and elements names.

    """    
    with open(dirname + '/POTCAR', 'w') as potcar:
        for s in elements:
            potcar_file_name = default_potcar_path + "/POTCAR_" + s
            shutil.copyfileobj(open(potcar_file_name), potcar)

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



    # POTCAR



    # INCAR
    shutil.copyfile(incar, dirname + "/INCAR")

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
