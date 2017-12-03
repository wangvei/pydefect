#!/usr/bin/env python                                                            

import numpy as np
import argparse 
import parse_poscar as ppos
import parse_potcar as ppot
import atom
import itertools as it
import subprocess as sp
import sys  
import os

def _get_num(x):                                                                  
    return float(''.join(i for i in x if i.isdigit() or i == '.'))   

def parseDefectsIn(defects_in_name):
    defects_in = open(defects_in_name)
    host_atoms = {} # host_atoms[name(str)] = [int(rep), [charge(int)]]
    dopants = {} # dopants [name(str)] = [charge(int)]
    interstitial_site = [] # [[i1(float)], [i2(float)], .. ]
    anti_site = [] # e.g. [Mg1_O1, O1_Mg1, ...]

    while True:
        line = defects_in.readline().split()

        if line == []: continue 
        if line[0] == "Name:":
            if line[1][-1].isdigit():
                host_atoms[line[1]] = []
                # Representative atom index
                host_atoms[line[1]].append(int(defects_in.readline().split()[1]))
                # skip 3 lines
                for i in range(3): defects_in.readline()
                # Charge list
                host_atoms[line[1]].append(
                            [int(i) for i in defects_in.readline().split()[1:]])
            else:
                defects_in.readline()
                dopants[line[1]] = \
                           [int(i) for i in defects_in.readline().split()[1:]]

        if line[0] == "Int_site:": 
            b = [_get_num(line[i]) for i in range(1, len(line))] 
            interstitial_site = [b[i:i + 3] for i in range(0, len(b), 3)] 
        if line[0] == "Antisite:": anti_site = line[1:]
        if line[0] == "Sym_break:": sym_break = line[1]
        if line[0] == "Irregular:": 
            irredular_defects = line[1:]
            break

    return host_atoms, dopants, interstitial_site, anti_site, \
                                                   sym_break, irredular_defects

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

def make_dir_poscar(inserted, removed, charge, poscar, r_param,
                    removed_atom=False, added_coord=False, 
                    added_atom_symbol=False):

    name = inserted + "_" + removed + "_" + str(charge)

    if os.path.isdir(name): 
        print "%12s alreadly exist, so is not made."% name
    else: 
        sp.call(["mkdir", name])
        ppos.make_point_defects(poscar, removed_atom, added_coord,
                                added_atom_symbol, name + '/POSCAR',
                                r_param, header=name)

# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser("make input for point defect calculations")
parser.add_argument("-d", "--defects", dest="defects", nargs="+",
                    type=str, help="Input particular defect(s).")
parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                    type=str, help="POSCAR name.")
parser.add_argument("-pot", "--potcar", dest="potcar", default="POTCAR",
                    type=str, help="POTCAR name.")
parser.add_argument("--defect_in", dest="input", default="defects.in",
                    type=str, help="defects.in file.")
parser.add_argument("--cutoff", dest="cutoff", default=4.0,
                    type=float, help="Radial cutoff [A] for randomization.")
parser.add_argument("--distance", dest="distance", default=0.2,
                    type=float, help="Max disp. distance for randomization.")

opts = parser.parse_args()

# name has the element symbol with a number. E.g., Mg1
# host_atoms[name(str)] = [int(rep), [charge(int)]] 
host_atoms, dopants, interstitial_site, anti_site, sym_break, irredular_defects\
                                                   =  parseDefectsIn(opts.input)
if sym_break[0] == "T": 
    r_param = [opts.cutoff, opts.distance]
else:
    r_param = False

if opts.defects:

    for a in opts.defects:

        added_atom, removed_atom, charge = [x for x in a.split("_")]    
        charge = int(charge)
    
        if removed_atom[0] == "i":
            coord = interstitial_site[removed_atom[-1]]
        else:
            coord = False

        if added_atom == "Va":
            make_dir_poscar(added_atom, removed_atom, charge,
                            opts.poscar, r_param, 
                            removed_atom=host_atoms[removed_atom][0],
                            added_coord=coord)
        else:
            make_dir_poscar(added_atom, removed_atom, charge,
                            opts.poscar, r_param, 
                            removed_atom=host_atoms[removed_atom][0],
                            added_coord=coord,
                            added_atom_symbol=added_atom)
        exit(0)

# Make perfect directory
if os.path.isdir("perfect"):  
    print "     perfect alreadly exist, so is not made."
else:
    sp.call(["mkdir", "perfect"])
    sp.call(["cp", opts.poscar, "perfect/POSCAR"])

# Element symbols are obtained from POTCAR, also used for intersititals.
element_symbols = ppot.parsePOTCAR(opts.potcar)[-1]

# Begin constructon of directories and POSCAR files for defects.
# Vacancies. Charges are sign reversed.
for removed_atom in host_atoms:
    for charge in host_atoms[removed_atom][1]:

        make_dir_poscar("Va", removed_atom, -charge, opts.poscar, r_param,
                        removed_atom=host_atoms[removed_atom][0])
# Interstitials. 
# Need to obtain the element symbols from host_atoms 
# E.g., Mg1, Mg2, O1, O2 -> Mg, O
# host_elements[Mg] = {0, 1, 2} set
host_elements = get_elements_from_host_atoms(host_atoms) 

for element in host_elements: 
    for interstitial_index, coord in enumerate(interstitial_site):
        for charge in host_elements[element]:

            make_dir_poscar(element, "i" + str(interstitial_index + 1), charge,
                            opts.poscar, r_param, added_coord=coord,
                            added_atom_symbol=element)
# Antisites. 
for a in anti_site: 
    # E.g., Mg1_O2 -> insert = Mg1, remove = O2
    added_atom, removed_atom = [x for x in a.split("_")]
    added_atom_charge = host_elements[added_atom]
    removed_atom_charge = [-i for i in host_elements[removed_atom[0:-1]]]

    # Construct the possible charge from charge_b and charge_c
    charges =list(set([sum(i) for i in 
                     list(it.product(added_atom_charge, removed_atom_charge))]))

    for charge in charges:
        make_dir_poscar(added_atom, removed_atom, charge, opts.poscar, r_param, 
                        removed_atom=host_atoms[removed_atom][0], 
                        added_atom_symbol=added_atom)

# Dopants. Both substitutional and intersitials are considered.
# dopants[Mg] = [0, 1, 2]
for dopant in dopants:
#for dopant, dopant_charges in enumerate(dopants):
    # substituted
    for host_removed in host_atoms:
        host_atom_charge = [-i for i in host_atoms[host_removed][1]]
        charges = list(set([sum(i) for i in list(
                               it.product(dopants[dopant], host_atom_charge))]))
        for charge in charges:

            make_dir_poscar(dopant, host_removed, charge, opts.poscar, r_param,
                            removed_atom=host_atoms[host_removed][0], 
                            added_atom_symbol=dopant)
    # Interstitials. 
    for interstitial_index, coord in enumerate(interstitial_site):
        for charge in dopants[dopant]:

            make_dir_poscar(dopant, "i" + str(interstitial_index + 1), charge, 
                            opts.poscar, r_param, added_coord=coord, 
                            added_atom_symbol=dopant)
sys.exit(0)        
