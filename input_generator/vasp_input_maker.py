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
from input_maker import DefectInputMaker, DefectInputSetMaker, extended_range, \
              print_already_exist, print_is_constructed,  perturb_around_a_point

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
    return sum([v.nelectrons * nions[i] for i, v in enumerate(p)]) - charge


class VaspDefectInputMaker(DefectInputMaker):
    """
    Constructs a set of vasp input files.
    POSCARs are constructed automatically.
    POTCARs are fetched from ~/.pydefect.yaml
    
    Args:

        defect_setting: DefectSetting class object
        incar (str): INCAR name
        kpoints (str): KPOINTS name
    """

    def __init__(self, defect_name, defect_setting, incar="INCAR", 
                 kpoints="KPOINTS"):

        # Construct:
        # is_directory, defect_name, defect_setting, in_name, out_name, charge
        super().__init__(defect_name, defect_setting)

        self.incar = incar
        self.kpoints = kpoints

        # Construct: defect_structure, defect_coords, defect_index
        dir_name = self.defect_name + "/"

        if self.is_directory == True:
            print_already_exist(dir_name)
        else:
            print_is_constructed(dir_name)

            elements = self.defect_structure.symbol_set
            nions = get_nions(self.defect_structure)

            # Construct defect_name directory
            os.makedirs(dir_name)
            # Construct defect.json
            self.make_json(dir_name + "defect.json")
            # Construct POSCAR-type files
            # Three POSCAR-type files are created:
            # POSCAR-Initial: POSCAR with a defect
            # POSCAR-DispInitial: POSCAR with a defect and perturbation near
            #                      a defect
            # POSCAR: Same as POSCAR-DispInitial if neighboring atoms are
            #         perturbed, otherwize POSCAR-Initial
            self.defect_structure.to(filename=dir_name + "POSCAR-Initial")
            if not self.defect_setting.distance == 0.0:
                self.perturbed_defect_structure, self.perturbed_sites = \
                            perturb_around_a_point(self.defect_structure,
                                                   self.defect_coords,
                                                   self.defect_setting.cutoff,
                                                   self.defect_setting.distance)
                self.perturbed_defect_structure.to(
                                      filename= dir_name + "POSCAR-DispInitial")
                self.perturbed_defect_structure.to(filename=dir_name + "POSCAR")
            else:
                self.defect_structure.to(filename=dir_name + "POSCAR")
    
            # Construct POTCAR file
            make_potcar(dir_name, elements, potcar_dir())
            # Construct INCAR file
            shutil.copyfile(self.incar, dir_name + "INCAR")
            nelect = get_charge(dir_name + "POTCAR", nions, self.charge)
            with open(dir_name + 'INCAR', 'a') as fa:
                fa.write('NELECT = ' + str(nelect))
            # Construct KPOINTS file
            shutil.copyfile(self.kpoints, dir_name + "KPOINTS")


class VaspDefectInputSetMaker(DefectInputSetMaker):

    def __init__(self, defect_setting, incar="INCAR", kpoints="KPOINTS"):

        super().__init__(defect_setting)

        self.incar = incar
        self.kpoints = kpoints

        # Construct a set of defect names.
        self._make_perfect_input_files()
        for d in self.defect_set:
            VaspDefectInputMaker(d, defect_setting, incar, kpoints)

    def _make_perfect_input_files(self):
        dir_name = "perfect/"
        if os.path.exists(dir_name):
            print_already_exist(dir_name)
        else:
            print_is_constructed(dir_name)
            os.makedirs(dir_name)
            self.defect_setting.structure.to(filename=dir_name + "POSCAR")
            shutil.copyfile(self.incar, dir_name + "INCAR")
            shutil.copyfile(self.kpoints, dir_name + "KPOINTS")
            make_potcar(dir_name, self.elements, potcar_dir())


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
            a = VaspDefectInputMaker(d, defect_setting, opts.incar, opts.kpoints)
            if a.is_directory == True:
                print_already_exist(d)
            else:
                print_is_constructed(d)
                a.constructor()
    else:
        VaspDefectInputSetMaker(defect_setting, incar=opts.incar,
                          kpoints=opts.kpoints)

if __name__ == "__main__": main()
