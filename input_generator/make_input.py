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

import atom
import itertools as it
import sys  

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

#class MakeParticularInput():

class MakeInput():

    def __init__(self, vacancy_sites=[], antisites=[], interstitial_sites=[], 
                 dopants=[], incar="INCAR", poscar="DPOSCAR", kpoints="KPONTS", run="run.sh",
                 cutoff=4.0, rand_distance=0.2):
        """
        Args:
            atomic_sites: a dictionary of atomic site index.
                          {"Mg1": 1, "O1":33}
            substitutions (antisites & substituted dopants): ["Al_Mg_0", "Mg_O1_2"]
            interstitial_sites: a list of lists with intersitial sites.
                               [[0,0,0], [0.5, 0.5, 0.5]]
                         ---->    i1           i2
            dopants (for intersititals): ["Al", "Ga"]
            oxidation_states: a dictionary of oxidation states for host and doped atoms.
                             {"Mg": 2, "O":-2, "Al": 3}

        + POTCAR files are fetched from ~/.pydefect.yaml

        """
        # TODO1: check if the input files exist.

        structure = Structure.from_file(poscar)                                 

       # construct perfect
        self.make_dir_three_input_files("perfect", incar, kpoints, run)
        shutil.copyfile(poscar, "perfect/POSCAR") 

       # construct vacancies
        
    def make_dir_three_input_files(self,dirname,incar, potcar, kpoints, run):
        if not os.path.exists(dirname):
            os.makedirs(dirname) 
            shutil.copyfile(incar,dirname + "/INCAR")
            shutil.copyfile(potcar,dirname + "/POTCAR")
            shutil.copyfile(kpoints,dirname + "/KPOINTS")
            shutil.copyfile(run, dirname + "/" + run)

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in="defect.in"):
        """
        Construct 
        """

        defect_in = open(defect_in)                                          
        atomic_sites = {} # host_atoms[name(str)] = [int(rep), [charge(int)]]         
        dopant_charges = {} # [name(str), ..]
                                                                                    
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
                dopants_charge[name] = int(line[1])

            elif line[0] == "Int_site:":                                              
                b = [_get_num(line[i]) for i in range(1, len(line))]                
                interstitial_sites = [b[i:i + 3] for i in range(0, len(b), 3)]       
            elif line[0] == "Antisite:": anti_sites = line[1:]                         
            elif line[0] == "Dopant_site:": dopant_sites = line[1:]                         
            elif line[0] == "Sym_break:": sym_break = line[1]                         
            elif line[0] == "Irregular:": irredular_defects = line[1:]                                        

                break                                                               
        substitutions = anti_sites + dopant_sites
                                                                                    
        return cls(structure, dopants=dopants, interstitials=interstitials, 
                   is_antisite=is_antisite, ElNeg_diff=ElNeg_diff, 
                   symbreak=symbreak, symprec=symprec)

                                                                                   
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

def get_num(x):                                                                
    return float(''.join(i for i in x if i.isdigit() or i == '.'))  
