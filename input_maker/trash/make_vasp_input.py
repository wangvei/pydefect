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
from pydefect.input_maker.defect_input import DefectSetting 

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

def make_potcar(dirname, elements, default_potcar_dir):
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
    """ 
    Return integer number from a string.
    """
    return int(''.join(i for i in x if i.isdigit() or i == '.'))


def _defect_name(defect_name):
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

def get_nions(defect_structure):
    """
    Return numbers of ions for elements in defect_structure. 
    """
    nions = [int(i) 
             for i in defect_structure.to(fmt="poscar").split("\n")[6].split()]
    return nions

def get_charge(potcar, nions, charge):
    """
    Return total charge from POTCAR file, number of ions, and charge state.
    """
    p = Potcar.from_file(potcar)
    # check only the length of potcar and nions.
    if not len(p) == len(nions):
        raise ValueError("Size of elements in POTCAR file is different")
    nelect = sum([v.nelectrons * nions[i] for i, v in enumerate(p)]) - charge
    return nelect        

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
        defect_setting: DefectInitialSetting class object
    """
    def __init__(self, defect_name, defect_setting):

        if os.path.exists(defect_name):
            self.is_directory = True
        else:
            self.is_directory = False
            self.defect_name = defect_name
            self.defect_setting = defect_setting
            self.in_name, self.out_name, self.charge = _defect_name(defect_name)
            # deepcopy is needed for structure
            self.defect_structure = deepcopy(defect_setting.structure)
            # e.g., irreducible_site_names = ["Mg1", "O1"]
            self.irreducible_site_names = [irreducible_site.irreducible_name 
                       for irreducible_site in self.defect_setting.irreducible_sites]

    def analyze_defect_name(self):
        # analyze out_name
        if re.match(r'^i[0-9]+$', self.out_name):
            interstitial_index = _get_int_from_string(self.out_name)
            try:
                defect_coords = \
                self.defect_setting.interstitial_coords[interstitial_index - 1]
            except:
                raise ValueError(
                 "Interstitial # {} is not defined".format(interstitial_index))
        elif self.out_name in self.irreducible_site_names:
            # There may be multiple candidates for inserted element.
            for irreducible_site in self.defect_setting.irreducible_sites:
                if self.out_name == irreducible_site.irreducible_name:
                    removed_atomic_index = irreducible_site.first_index
                    defect_coords = irreducible_site.repr_coords
            self.defect_structure.remove_sites([removed_atomic_index - 1])
        else:
            raise ValueError("{} is improper.".format(self.out_name))

        # analyze in_name
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
            raise ValueError("{} is improper.".format(self.in_name))

        def write_input_files(self):
            """
            Needs to be modified by subclasses depending on the code.
            """
            if self.is_directory == False:
                os.makedirs(self.defect_name)
                # write a defect position to defect.json file.
                with open(self.defect_name + "/defect.json", 'w') as fw:
                    json.dump({"defect_index": self.defect_index,
                               "defect_coords": self.defect_coords,
                               "in_name": self.in_name,
                               "out_name": self.out_name,
                               "charge": self.charge}, fw, indent=2)
            else:
                print("{} exists. is not defined".format(interstitial_index))
    
class VaspInputMaker(InputMaker):
    """
    Construct a set of vasp input files.
    POTCAR files are fetched from ~/.pydefect.yaml
    
    Args:
        defect_name (str): defect name defined in PyDefect, e.g., "Va_Mg1_2"
        defect_setting: DefectInitialSetting class object
        incar (str): INCAR name
        kpoints (str): KPOINTS name
    """

    def __init__(self, defect_name, defect_setting, incar="INCAR", 
                 kpoints="KPOINTS"):

        super().__init__(defect_name, defect_setting)
        self.incar = incar
        self.kpoints = kpoints

    def write_input_files(self):

        analyze_defect_name()            
        super().write_input_files()
        self.defect_structure.to(filename=self.defect_name + "/POSCAR-Initial")

        # vasp input
        # perturb neighboring atoms randomly.
        if self.defect_setting.displace is not None: 
            self.defect_structure, self.perturbed_sites = \
                   perturb_around_a_point(self.defect_structure, defect_coords, 
                      self.defect_setting.cutoff, self.defect_setting.displace)
            self.defect_structure.to(
                             filename=self.defect_name + "/POSCAR-DispInitial")
        
        self.defect_structure.to(filename=self.defect_name + "/POSCAR")

        elements = self.defect_structure.symbol_set
        nions = get_nions(self.defect_structure)
        # Construct POTCAR file
        make_potcar(self.defect_name, elements, potcar_dir()) 
        # Construct INCAR file
        shutil.copyfile(self.incar, self.defect_name + "/INCAR")
        nelect = get_charge(self.defect_name + "/POTCAR", nions, self.charge)
        with open(self.defect_name + '/INCAR', 'a') as i:
            i.write('NELECT = ' + str(nelect))
        # Construct KPOINTS file
        shutil.copyfile(self.kpoints, self.defect_name + "/KPOINTS")


class VaspInputSetMaker():

    def __init__(self, defect_setting, incar="INCAR", kpoints="KPOINTS"):

        #TODO: check INCAR KPOINTS
        
        if not os.path.exists(incar):
            raise IOError('{} does not exist.'.format(incar))
        if not os.path.exists(kpoints):
            raise IOError('{} does not exist.'.format(kpoints))

        self.defect_setting = defect_setting
        self.incar = incar
        self.kpoints = kpoints
        self.elements = self.defect_setting.structure.symbol_set

        self._perfect_constructor()
        self.defect_set = self._vacancy_setter() + self._interstitial_setter() + \
                          self._antisite_dopant_setter()
        for i in self.defect_setting.include:
            self.defect_set.append(i)

        for e in self.defect_setting.exclude:
            if e in self.defect_set:
                self.defect_set.remove(e)
            else:
                print("{} does not exist.".format(e))

        for d in self.defect_set:
            a = VaspInputMaker(d, self.defect_setting, self.incar, 
                               self.kpoints)
            if a.is_directory == True:
                _print_already_exist(d)
            else:
                _print_is_constructed(d)
                a.constructor()


    def _perfect_constructor(self):
        perfect = "perfect"
        if os.path.exists(perfect):
            _print_already_exist(perfect)
        else:
            _print_is_constructed(perfect)
            os.makedirs(perfect)
            self.defect_setting.structure.to(filename=perfect + "/POSCAR")
            shutil.copyfile(self.incar, perfect + "/INCAR")
            shutil.copyfile(self.kpoints, perfect + "/KPOINTS")
            make_potcar(perfect, self.elements, potcar_dir()) 

    def _vacancy_setter(self):
        name_set = []
        for i in self.defect_setting.irreducible_sites:
            oxidation_state = self.defect_setting.oxidation_states[i.element]
            for o in extended_range(oxidation_state):
                defect_name = "Va_" + i.irreducible_name + "_" + str(-o)
                name_set.append(defect_name)
        return name_set

    def _interstitial_setter(self):
        name_set = []
        for e in self.elements:
            oxidation_state = self.defect_setting.oxidation_states[e]
            for j in range(len(self.defect_setting.interstitial_coords)):
                for o in extended_range(oxidation_state):
                    defect_name = e + "_i"  + str(j + 1) + "_" + str(o)
                    name_set.append(defect_name)
        return name_set

    def _antisite_dopant_setter(self):
        name_set = []
        for a in self.defect_setting.antisite_configs + \
                                            self.defect_setting.dopant_configs:
            in_element, out_element = a.split("_")
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


def main():
    import argparse
    parser = argparse.ArgumentParser(
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--defectin", dest="defectin",
                        default="defect.in", type=str, help="defect.in name.")
    parser.add_argument("-p", "--dposcar", dest="dposcar", default="DPOSCAR",
                        type=str, help="DPOSCAR name.")
    parser.add_argument("--incar", dest="incar", default="INCAR", 
                        type=str, help="INCAR name.")
    parser.add_argument("--kpoints", dest="kpoints", default="KPOINTS", 
                        type=str, help="KPOINTS name.")
    parser.add_argument("--add", dest="add", type=str, nargs="+", 
                        help="Particular defect names added.")

    opts = parser.parse_args()
    defect_setting = DefectSetting.from_defect_in(poscar=opts.dposcar, 
                                                  defect_in_file=opts.defectin)
    if opts.add:
        for d in opts.add:
            a = VaspInputMaker(d, defect_setting, opts.incar, opts.kpoints)
            if a.is_directory == True:
                _print_already_exist(d)
            else:
                _print_is_constructed(d)
                a.constructor()
    else:
        VaspInputSetMaker(defect_setting, incar=opts.incar, 
                          kpoints=opts.kpoints)

if __name__ == "__main__": main()
