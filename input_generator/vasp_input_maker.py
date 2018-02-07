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
from defect_in import DefectSetting
from input_maker import InputMaker, extended_range, _print_already_exist, \
                        _print_is_constructed

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pydefect.yaml")

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

class VaspInputMaker(InputMaker):
    """
    Constructs a set of vasp input files.
    POSCARs are constructed automatically.
    POTCARs are fetched from ~/.pydefect.yaml
    
    Args:
        defect_name (str): defect name defined in PyDefect, e.g., "Va_Mg1_2"
        defect_setting: DefectSetting class object
        incar (str): INCAR name
        kpoints (str): KPOINTS name
    """

    def __init__(self, defect_name, defect_setting, incar="INCAR", 
                 kpoints="KPOINTS"):

        super().__init__(defect_name, defect_setting)
        self.incar = incar
        self.kpoints = kpoints

    def make_vasp_defect_input_files(self):
        self.analyze_defect_name()
        if self.is_directory == True:
            _print_already_exist(self.defect_name)
        else:
            _print_is_constructed(self.defect_name)

            self.make_directory_json()
            self.defect_structure.to(
                                 filename=self.defect_name + "/POSCAR-Initial")

            if not self.defect_setting.distance == 0.0:
                self.make_perturbed_defect_structure()
                self.perturbed_defect_structure.to(
                             filename=self.defect_name + "/POSCAR-DispInitial")
                self.perturbed_defect_structure.to(
                             filename=self.defect_name + "/POSCAR")
            else:
                self.defect_structure.to(filename=self.defect_name + "/POSCAR")
    
            elements = self.defect_structure.symbol_set
            nions = get_nions(self.defect_structure)
            # Construct POTCAR file
            make_potcar(self.defect_name, elements, potcar_dir()) 
            # Construct INCAR file
            shutil.copyfile(self.incar, self.defect_name + "/INCAR")
            nelect = get_charge(self.defect_name + "/POTCAR", 
                                                            nions, self.charge)
            with open(self.defect_name + '/INCAR', 'a') as i:
                i.write('NELECT = ' + str(nelect))
            # Construct KPOINTS file
            shutil.copyfile(self.kpoints, self.defect_name + "/KPOINTS")


class VaspInputSetMaker():

    def __init__(self, defect_setting, incar="INCAR", kpoints="KPOINTS"):

        if not os.path.exists(incar):
            raise IOError('{} does not exist.'.format(incar))
        if not os.path.exists(kpoints):
            raise IOError('{} does not exist.'.format(kpoints))

        self.defect_setting = defect_setting
        self.incar = incar
        self.kpoints = kpoints
        self.elements = self.defect_setting.structure.symbol_set

        self.make_vasp_perfect_input_files()

        self.defect_set = self._vacancy_setter() + \
                          self._interstitial_setter() + \
                          self._substitutional_setter()

        for i in self.defect_setting.included:
            self.defect_set.append(i)

        for e in self.defect_setting.excluded:
            if e in self.defect_set:
                self.defect_set.remove(e)
            else:
                print("{} does not exist.".format(e))

        print(self.defect_set)
        for d in self.defect_set:
            a = VaspInputMaker(d, self.defect_setting, incar, kpoints)
            a.make_vasp_defect_input_files()

    def make_vasp_perfect_input_files(self):
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

    def _substitutional_setter(self):
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
