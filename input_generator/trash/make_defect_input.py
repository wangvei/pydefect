#!/usr/bin/env python

import numpy as np
import argparse 
import parse_poscar as ppos
import parse_potcar as ppot
import atom
import itertools as it
import sys  

__author__ = "Yu Kumagai"
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

def _neatly_output_numbers(list_number):
    """
    Eg, [1, 2, 3, 5, 6, 8, 11, 12, 13] -> 1..3, 5, 6, 8, 11..13
    """
    length = len(list_number)
    out = str(list_number[0])

    is_dot = False

    for i in range(1,length - 1):
        a= list_number[i - 1]  
        b= list_number[i]  
        c= list_number[i + 1]  
        if b - a == 1 and c - b == 1:
            if is_dot == False:
                out += ".."
                is_dot = True
        else:
            if is_dot == False:
                out += " " + str(b)
            else:
                out += str(b)
                is_dot = False

    if is_dot == False: out += " " + str(list_number[length - 1])
    else:               out += str(list_number[length - 1])

    return out

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                    type=str, help="POSCAR name.")
parser.add_argument("--potcar", dest="potcar", default="POTCAR",
                    type=str, help="POTCAR name.")
parser.add_argument("-d","--dopants", dest="dopants", default=False, nargs="+", 
                    type=str, help="Dopant elements. Eg. Al Ga In.")
parser.add_argument("-i", dest="interstitials", default=False, nargs="+",
                    type=float, help="Inetrstitials. Eg. 0 0 0  0.5 0.5 0.5.")
parser.add_argument("-a","--antisite", dest="antisite", action="store_false",
                    help="Set if antisites are considered.")
parser.add_argument("-e","--endiff", dest="endiff", type=float, default=1.0,
                    help="Ectronegativity diff for antisites and impurities.")
parser.add_argument("-s","--symbreak", dest="symbreak", action="store_true",
                    help="Set if symmetry is not broken.")
parser.add_argument("--symprec",dest="symprec", type=float,
                    help="Set the symprec [A].")

opts = parser.parse_args()

symbols, num_atoms, atomic_pos = ppos.parsePOSCAR(opts.poscar)[1:4]

if symbols[0] == "X":
    symbols = ppot.parsePOTCAR(opts.potcar)[8]

symbols_all = []
for i in range(len(num_atoms)):
    for j in range(num_atoms[i]):
        symbols_all.append(symbols[i])

# Get the dictionary of the atomic mapping.
atom_map = ppos.get_atom_mapping(opts.poscar,symprec=opts.symprec)

# *name* includes the element name with irrep number.
name = {}
for i in sorted(atom_map):
    symbol = symbols_all[i - 1]
    num = 1    

    while True:
        if symbol + str(num) in name.values(): 
            num += 1
        else:
            break

    if num > 9:
        print "The number of irreducible sites must be less than 10."
        sys.exit(1)

    name[i] = symbol + str(num) 

# Obtain electron negativity and formal charges.
electron_negativity = {}
charge = {}

for i in name:
    try: 
        electron_negativity[name[i][0:-1]] = \
                                  float(atom.electron_negativity[name[i][0:-1]])
    except:
        print "The electron negativity of", name[i][0:-1], "cannot be obtained."
        electron_negativity[name[i][0:-1]] = "N.A."

    try: 
        charge[name[i][0:-1]] = int(atom.formal_charges[name[i][0:-1]])
    except:
        print "The charge of", name[i][0:-1], "cannot be obtained."
        charge[name[i][0:-1]] = "N.A."

# Note that "lattice_vectors" is a list so that the index should be reduced by 1.
for i in sorted(atom_map):
    print " Name:", name[i]
    print "  Rep:", i
    print " Equi:", _neatly_output_numbers(atom_map[i])
    print "Coord: %7.5f %7.5f %7.5f" % tuple(atomic_pos[i-1])
    print "ElNeg:", electron_negativity[name[i][0:-1]]
    print "Charg:",
    if charge[name[i][0:-1]] > 0:
        print ' '.join(['%-2s' % (str(i),) 
                                  for i in range(0, charge[name[i][0:-1]] + 1)])
    else:
        print ' '.join(['%-2s' % (str(i),) 
                                  for i in range(charge[name[i][0:-1]], 1)])
    print ""

if opts.interstitials:
    if not len(opts.interstitials) % 3 == 0:
        print "The interstitial coordinates are not proper. Bye."
        sys.exit(1)
    print " Int_site:", tuple([opts.interstitials[i:i + 3]
                                 for i in range(0, len(opts.interstitials), 3)])

if opts.antisite:
    print " Antisite:",

    antisite = []
    for i in list(it.combinations(name.keys(), 2)):
        if name[i[0]][0:-1] == name[i[1]][0:-1]: continue
    
        ENdiff = abs(electron_negativity[name[i[0]][0:-1]] - \
                     electron_negativity[name[i[1]][0:-1]])
        if ENdiff < opts.antisite:
            name1 = name[i[0]][0:-1] + "_" + name[i[1]]
            name2 = name[i[1]][0:-1] + "_" + name[i[0]]

            if not name1 in antisite: antisite.append(name1)
            if not name2 in antisite: antisite.append(name2)

    print " ".join(antisite)

if opts.dopants:
    for i, elem in enumerate(opts.dopants):
        try: 
            electron_negativity[elem] = atom.electron_negativity[elem]
        except:
            print "The electron negativity of", elem, "cannot be obtained."
            electron_negativity[elem] = "N.A."
    
        try: 
            charge[elem] = atom.formal_charges[elem]
        except:
            print "The electron negativity of", elem, "cannot be obtained."
            charge[elem] = "N.A."
    
    for i in range(len(opts.dopants)):
        print " Name:", opts.dopants[i]
        print "ElNeg:", electron_negativity[opts.dopants[i]]
        print "Charg:",
        if charge[opts.dopants[i]] > 0:
            print ' '.join(['%-2s' % (str(i),) 
                                for i in range(0, charge[opts.dopants[i]] + 1)])
        else:
            print ' '.join(['%-2s' % (str(i),) 
                                for i in range(charge[opts.dopants[i]],1)])
        print ""

    print "Dopant_site:",
    for i in opts.dopants:
        for j, k in name.iteritems():
            ENdiff = abs(electron_negativity[i] - electron_negativity[k[0:-1]])
            if ENdiff < opts.endiff:
                print i + "_" + k,
    print ""

print "Sym_break:", opts.symbreak
print "Irregular:"
