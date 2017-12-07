#!/usr/bin/env python
from __future__ import print_function
import os
import shutil
import numpy as np
import warnings
import argparse 
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
#from pymatgen.core.periodic_table import Element, Specie, get_el_sp

import atom
import itertools as it
import sys  
import ruamel.yaml as yaml

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

#class MakeParticularInput():

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pydefect.yaml")

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


class MakeInput():

    def __init__(self, atomic_sites={}, dopant_charges={}, substitutions=[], 
                 interstitial_sites=[], incar="INCAR", poscar="DPOSCAR", 
                 kpoints="KPONTS", run="run.sh", cutoff=4.0, rand_distance=0.2):
        """
        
        + POTCAR files are fetched from ~/.pydefect.yaml

        Terminology:
            name: specie name + irreducible atom index
          specie: specie name
          charge: oxidation state, a single integer number
             Rep: representative position for *name* in structure object
            --------------------------------------------------------------
        Args:
            atomic_sites: a dictionary of atomic site index.
                          {name : [Rep, charge], ...}
          dopant_charges: {specie : charge, ...}
           substitutions: sum of antisites and dopant_sites ["Al_Mg1", ...]
      interstitial_sites: a list of lists with intersitial sites.
                          [[0, 0, 0], [0.1, 0.1, 0.1], ...]
                     ---->     i1            i2
        """
        # TODO1: check if the input files exist.

        default_potcar_path = _laad_potcar_dir()
        structure = Structure.from_file(poscar)                                 


       # construct perfect
        self.make_dir_three_input_files("perfect", incar, kpoints, run)
        shutil.copyfile(poscar, "perfect/POSCAR") 
        make_POTCAR("perfect", atomic_sites.keys(), default_potcar_path)

       # construct vacancies

        self.make_dir_three_input_files("perfect", incar, kpoints, run)


    def make_POTCAR(self, potcar_path, species, default_potcar_path):
        with open(potcar_path'/POTCAR', 'w') as potcar:
            for s in species:
                potcar_file_name = default_potcar_path + "/POTCAR_" + s
                shutil.copyfileobj(open(potcar_file_name), outfile)

        
    def make_dir_three_input_files(self,dirname,incar, kpoints, run):
        if not os.path.exists(dirname):
            os.makedirs(dirname) 
            shutil.copyfile(incar,dirname + "/INCAR")
            shutil.copyfile(kpoints,dirname + "/KPOINTS")
            shutil.copyfile(run, dirname + "/" + run)


    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in="defect.in"):
        """
        Construct four variables:
            name: specie name + irreducible atom index
          specie: specie name
          charge: oxidation state, a single integer number
             Rep: representative position for *name* in structure object
            --------------------------------------------------------------
            atomic_sites: {name : [Rep, charge], ...}
          dopant_charges: {specie : charge, ...}
           substitutions: sum of antisites and dopant_sites ["Al_Mg1", ...]
      interstitial_sites: [[0, 0, 0], [0.1, 0.1, 0.1], ...]
        """

        defect_in = open(defect_in)                                          

        atomic_sites = {}
        dopant_charges = {}
                                                                                    
        while True:                                                                 
            line = defects_in.readline().split()                                    
                                                                                    
            if line == []: continue                                                 
            if line[0] == "Name:":                                                  
                name = line[1]
                atomic_sites[name] = [None, None]
                while True: 
                    line = defects_in.readline().split()
                    if line[0] == "Rep:":
                        atomic_sites[name][0] = int(line[1])
                    elif line[0] == "Charge:":
                        atomic_sites[name][1] = int(line[1])
                    elif line == "": break
                if atomic_sites[name][0] == None or atomic_sites[name][1] == None:
                    raise ValueError("defect.in type file is not appropriate.")

            elif line[0] == "Dopant":
                name = line[1]
                while True: 
                    line = defects_in.readline().split()
                    if line[0] == "Charge:": 
                        dopant_charges[name] = line[1]
                    elif line == "": break
                dopant_charges[name] = int(line[1])

            elif line[0] == "Int_site:":                                              
                b = [self._get_num(line[i]) for i in range(1, len(line))]                
                interstitial_sites = [b[i:i + 3] for i in range(0, len(b), 3)]       
            elif line[0] == "Antisite:": antisites = line[1:]                         
            elif line[0] == "Dopant_site:": dopant_sites = line[1:]                         
            elif line[0] == "Sym_break:": sym_break = line[1]                         
            elif line[0] == "Irregular:": irredular_defects = line[1:]                                        

                break                                                               
        substitutions = antisites + dopant_sites
                                                                                    
        return cls(self, atomic_sites={}, dopant_charges={}, substitutions=[], 
                   interstitial_sites=[], incar="INCAR", poscar="DPOSCAR", 
                   kpoints="KPONTS", run="run.sh", cutoff=4.0, rand_distance=0.2)

                                                                                   
    def get_elements_from_host_atoms(host_atoms):                                   
        """                                                                         
        host_atoms[name(str)] = [int(rep), [charge(int)]]                           
        Return element names with charges, which are one-size-fits-all.             
        E.g. host_atoms[Sn1] = [0, 1, 2],host_atoms[Sn2] = [2, 3, 4]                
             -> elements[Sn] = [0, 1, 2, 3, 4]                                      
        """                                                                         
        elements = {}                                                               
                                                                                    
        for atom in host_atoms:                                                     
            # Remove number. E.g., Mg1 -> Mg                                        
            name = atom[0:-1]                                                       
            if name not in elements:                                                
                elements[name] = set()                                              
                                                                                    
            for charge in host_atoms[atom][1]:                                      
                elements[name].add(charge)                                          
                                                                                    
        return elements    


def _get_num(x):                                                                
    return float(''.join(i for i in x if i.isdigit() or i == '.'))  
