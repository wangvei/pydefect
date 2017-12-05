#!/usr/bin/env python

import numpy as np
import warnings
import argparse 
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
#from pymatgen.core.periodic_table import Element, Specie, get_el_sp

#import parse_poscar as ppos
#import parse_potcar as ppot
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

class DefectIn():

    def __init__(self, structure, charges=False, dopants=False, interstitials=False, 
                 antisites=False, ElNeg_diff=1.0, irreguar=None, symbreak=True, symprec=1e-5):

        # only vasp5 POSCAR format is supported.
        self.dopants = dopants
        self.interstitials = interstitials
        self.antisites = antisites
        self.ElNeg_diff = ElNeg_diff
        self.symbreak = symbreak
        self.symprec = symprec

        symm_finder = SpacegroupAnalyzer(structure)
        self.symm_structure = symm_finder.get_symmetrized_structure()
        self.equiv_site_seq = self.symm_structure.equivalent_sites
        
        elements = []
        frac_coords = []
        # repr_frac_coords e.g.,[["Mg1", "Mg", 2, [0, 0, 0]],
        #                        ["Mg2", "Mg", 4, [0.5, 0.5, 0.5]]]
        repr_frac_coords = []
        are_elements = {}
        electron_negativity = {}
        oxidation_states = {}

        for s in structure.symbol_set: 
            print(s)
            are_elements[s] = 0
        
            try: 
                electron_negativity[s] = float(atom.electron_negativity[s])
            except:
                warnings.warn("The electron negativity of", s, "is unavailable.")
                electron_negativity[s] = "N.A."

            try: 
                oxidation_states[s] = Element(s).common_oxidation_states[-1]
            except:
                warnings.warn("The oxidation state of", s, "is unavailable.")
                oxidation_states[s] = "N.A."

        self.electron_negativity = electron_negativity
        self.oxidation_states = oxidation_states
                
        for inequiv_site in self.equiv_site_seq:
            element = inequiv_site[0].species_string
            are_elements[element] += 1 # increment nr of inequiv site
            nr_sites = len(inequiv_site)
            coords = inequiv_site[0].frac_coords
            repr_frac_coords.append([element+str(are_elements[element]), 
                                     element, nr_sites, coords])
            for equiv_site in inequiv_site:
                elements.append(equiv_site.species_string)
                frac_coords.append(equiv_site.frac_coords)

        # Atoms are sorted by symmetry.
        self.structure = Structure(structure.lattice, elements, frac_coords)
        self.repr_frac_coords = repr_frac_coords

        if antisites = True:
            x = []
            for i in list(it.combinations(name.keys(), 2)):
                if name[i[0]][0:-1] == name[i[1]][0:-1]: continue
            
                ENdiff = abs(electron_negativity[name[i[0]][0:-1]] - \
                             electron_negativity[name[i[1]][0:-1]])
                if ENdiff < opts.antisite:
                    name1 = name[i[0]][0:-1] + "_" + name[i[1]]
                    name2 = name[i[1]][0:-1] + "_" + name[i[0]]
        
                    if not name1 in antisite: antisite.append(name1)
                    if not name2 in antisite: antisite.append(name2)

        # Obtain electron negativity and formal charges.
        
        
#        # Note that "lattice_vectors" is a list so that the index should be reduced by 1.
#        for i in sorted(atom_map):
#            print " Name:", name[i]
#            print "  Rep:", i
#            print " Equi:", _pretty_printing_numbers(atom_map[i])
#            print "Coord: %7.5f %7.5f %7.5f" % tuple(atomic_pos[i-1])
#            print "ElNeg:", electron_negativity[name[i][0:-1]]
#            print "Charg:",
#            if charge[name[i][0:-1]] > 0:
#                print ' '.join(['%-2s' % (str(i),) 
#                                          for i in range(0, charge[name[i][0:-1]] + 1)])
#            else:
#                print ' '.join(['%-2s' % (str(i),) 
#                                          for i in range(charge[name[i][0:-1]], 1)])
#            print ""
#        
#        if opts.interstitials:
#            if not len(opts.interstitials) % 3 == 0:
#                print "The interstitial coordinates are not proper. Bye."
#                sys.exit(1)
#            print " Int_site:", tuple([opts.interstitials[i:i + 3]
#                                         for i in range(0, len(opts.interstitials), 3)])
#        
#        if opts.antisite:
#            print " Antisite:",
#        
#            antisite = []
#            for i in list(it.combinations(name.keys(), 2)):
#                if name[i[0]][0:-1] == name[i[1]][0:-1]: continue
#            
#                ENdiff = abs(electron_negativity[name[i[0]][0:-1]] - \
#                             electron_negativity[name[i[1]][0:-1]])
#                if ENdiff < opts.antisite:
#                    name1 = name[i[0]][0:-1] + "_" + name[i[1]]
#                    name2 = name[i[1]][0:-1] + "_" + name[i[0]]
#        
#                    if not name1 in antisite: antisite.append(name1)
#                    if not name2 in antisite: antisite.append(name2)
#        
#            print " ".join(antisite)
#        
#        if opts.dopants:
#            for i, elem in enumerate(opts.dopants):
#                try: 
#                    electron_negativity[elem] = atom.electron_negativity[elem]
#                except:
#                    print "The electron negativity of", elem, "cannot be obtained."
#                    electron_negativity[elem] = "N.A."
#            
#                try: 
#                    charge[elem] = atom.formal_charges[elem]
#                except:
#                    print "The electron negativity of", elem, "cannot be obtained."
#                    charge[elem] = "N.A."
#            
#            for i in range(len(opts.dopants)):
#                print " Name:", opts.dopants[i]
#                print "ElNeg:", electron_negativity[opts.dopants[i]]
#                print "Charg:",
#                if charge[opts.dopants[i]] > 0:
#                    print ' '.join(['%-2s' % (str(i),) 
#                                        for i in range(0, charge[opts.dopants[i]] + 1)])
#                else:
#                    print ' '.join(['%-2s' % (str(i),) 
#                                        for i in range(charge[opts.dopants[i]],1)])
#                print ""
#        
#            print "Dopant_site:",
#            for i in opts.dopants:
#                for j, k in name.iteritems():
#                    ENdiff = abs(electron_negativity[i] - electron_negativity[k[0:-1]])
#                    if ENdiff < opts.endiff:
#                        print i + "_" + k,
#            print ""
#        
#        print "Sym_break:", opts.symbreak
#        print "Irregular:"


    @classmethod
    def from_file(cls, poscar="POSCAR",dopants=False, interstitials=False, 
                  antisites=False, ElNeg_diff=1.0, symbreak=True, symprec=1e-5):
        """
        Construct DefectIn class object from a POSCAR file.
        """
        structure = Structure.from_file(poscar)
        return cls(structure, dopants=False, interstitials=False, 
                   antisites=False, ElNeg_diff=1.0, symbreak=True, symprec=1e-5)


    @classmethod
    def from_defect_in(cls, defectin="defect.in"):
        """
        Construct DefectIn class object from a defect.in file.
        """
        pass

    def to(self, filename1="defect.in",filename2="DPOSCAR"):
        """
        Print readable defect.in file.
        """
        _pretty_printing_defect_in(filename=filename1) 
        self.structure.to(fmt="poscar",filename=filename2)
                
    def _pretty_printing_defect_in(self,filename="defect.in"):
        pass



    def _pretty_printing_numbers(list_number):
        """
        Return pretty format of a series of numbers.
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

