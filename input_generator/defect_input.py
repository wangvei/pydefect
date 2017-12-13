#!/usr/bin/env python

#from __future__ import print_function
import itertools as it
import numpy as np
import sys
import warnings
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import atom

# from pymatgen.core.periodic_table import Element, Specie, get_el_sp

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


class DefectIn():

    def __init__(self, structure, dopants=[], interstitials=False,
                 is_antisite=False, ElNeg_diff=1.0, irregular=None, symbreak=True, symprec=1e-5):

        # only vasp5 POSCAR format is supported.
        self.dopants = dopants
        self.interstitials = interstitials
        self.is_antisite = is_antisite
        self.ElNeg_diff = ElNeg_diff
        self.irregular = irregular
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

        # intrinsic elements + dopants
        for s in structure.symbol_set + tuple(dopants):
            are_elements[s] = 0

            # Obtain electron negativity and oxidation states.
            try:
                electron_negativity[s] = float(atom.electron_negativity[s])
            except:
                warnings.warn("The electron negativity of " + s + " is unavailable.")
                electron_negativity[s] = "N.A."

            try:
                oxidation_states[s] = Element(s).common_oxidation_states[-1]
            except:
                warnings.warn("The oxidation state of " + s + " is unavailable.")
                oxidation_states[s] = "N.A."

        self.electron_negativity = electron_negativity
        self.oxidation_states = oxidation_states

        for inequiv_site in self.equiv_site_seq:
            element = inequiv_site[0].species_string
            are_elements[element] += 1  # increment nr of inequiv site
            nr_sites = len(inequiv_site)
            coords = inequiv_site[0].frac_coords
            repr_frac_coords.append([element + str(are_elements[element]),
                                     element, nr_sites, coords])
            for equiv_site in inequiv_site:
                elements.append(equiv_site.species_string)
                frac_coords.append(equiv_site.frac_coords)

        # Atoms are sorted by symmetry.
        self.structure = Structure(structure.lattice, elements, frac_coords)
        # repr_frac_coords = ["Mg1", "Mg", nr_sites, repr_coords]
        self.repr_frac_coords = repr_frac_coords

        if is_antisite:
            # ex.             Mg_O1
            # antisites = [[ "Mg", "O1"], ...]
            antisites = []
            for r in self.repr_frac_coords:
                for s in structure.symbol_set:
                    if r[1] == s: continue
                    try:
                        ENdiff = electron_negativity[r[1]] - electron_negativity[s]
                    except:
                        # It's a fictitious number
                        ENdiff = 100
                    if abs(ENdiff) < ElNeg_diff:
                        antisites.append([s, r[0]])
            self.antisites = antisites

        if dopants:
            dopant_sites = []
            for r in self.repr_frac_coords:
                for d in self.dopants:
                    ENdiff = electron_negativity[r[1]] - electron_negativity[d]
                    if abs(ENdiff) < ElNeg_diff:
                        dopant_sites.append([d, r[0]])
            self.dopant_sites = dopant_sites

    #    def _return_defect_type_by_ElNeg_diff(self, ElNeg_diff, a, b):
    #        kkkkkkkkkkk

    @classmethod
    def from_str_file(cls, poscar="POSCAR", dopants=[], interstitials=False,
                      is_antisite=False, ElNeg_diff=1.0, irregular=None, symbreak=True, symprec=1e-5):
        """
        Construct DefectIn class object from a POSCAR file.
        """

        structure = Structure.from_file(poscar)

        return cls(structure, dopants=dopants, interstitials=interstitials,
                   is_antisite=is_antisite, ElNeg_diff=ElNeg_diff, irregular=irregular,
                   symbreak=symbreak, symprec=symprec)

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defectin="defect.in"):
        """
        Construct DefectIn class object from a defect.in file.
        """
 
        structure = Structure.from_file(poscar)
 
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
 
        return cls(structure, dopants=dopants, interstitials=interstitials,
                   is_antisite=is_antisite, ElNeg_diff=ElNeg_diff,
                   symbreak=symbreak, symprec=symprec)

    #    def get_elements_from_host_atoms(host_atoms):
    #        """
    #        host_atoms[name(str)] = [int(rep), [charge(int)]]
    #        Return element names with charges, which are one-size-fits-all.
    #        E.g. host_atoms[Sn1] = [0, 1, 2],host_atoms[Sn2] = [2, 3, 4]
    #             -> elements[Sn] = [0, 1, 2, 3, 4]
    #        """
    #        elements = {}
    #
    #        for atom in host_atoms:
    #            # Remove number. E.g., Mg1 -> Mg
    #            name = atom[0:-1]
    #            if name not in elements:
    #                elements[name] = set()
    #
    #            for charge in host_atoms[atom][1]:
    #                elements[name].add(charge)
    #
    #        return elements

    def to(self, filename1="defect.in", filename2="DPOSCAR"):
        """
        Print readable defect.in file.
        """
        self._print_defect_in(filename=filename1)
        self.structure.to(fmt="poscar", filename=filename2)

    def print_defect_in(self, filename="defect.in"):
        i = 1
        #        print(self.repr_frac_coords)
        for r in self.repr_frac_coords:
            print("  Name: {}".format(r[0]))
            print("   Rep: {}".format(str(i)))
            print(" Equiv: {}".format(str(i) + ".." + str(i + r[2] - 1)))
            i += r[2]
            print(" Coord: %9.7f %9.7f %9.7f" % tuple(r[3]))
            self._print_ElNeg(self.electron_negativity[r[1]])
            self._print_charges(self.oxidation_states[r[1]])
            print("")

        if self.interstitials:
            coords = [float(i) for i in self.interstitials.split()]
            if not len(coords) % 3 == 0:
                raise ValueError("The interstitial coordinates are not proper.")
            print("Int_site: ", end="")
            print([coords[i:i + 3] for i in range(0, len(coords), 3)])

        if self.is_antisite:
            print("Antisite: ", end="")
            print(' '.join(i[0] + "_" + i[1] for i in self.antisites), end='\n\n')

        if self.dopants:
            for d in self.dopants:
                print("Dopant: {}".format(d))
                self._print_ElNeg(self.electron_negativity[d])
                self._print_charges(self.oxidation_states[d])
            print("")
            print("Dopant_site: ", end="")
            print(' '.join(i[0] + "_" + i[1] for i in self.dopant_sites), end='\n\n')

        print("Sym_break: {}".format(self.symbreak))

        if self.irregular:
            print("Irregular: {}".format(self.irregular))
        else:
            print("Irregular: ")

    def _print_ElNeg(self, ElNeg):
        print("EleNeg: {}".format(ElNeg))

    def _print_charges(self, os):
        print("Charge: {}".format(os))

class NoElementNameError(Exception):
    pass 

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name.")
    parser.add_argument("--potcar", dest="potcar", default="POTCAR",
                        type=str, help="POTCAR name.")
    parser.add_argument("-d","--dopants", dest="dopants", default="", nargs="+", 
                        type=str, help="Dopant elements. Eg. Al Ga In.")
    parser.add_argument("-i", dest="interstitials", default=False, nargs="+",
                        type=float, help="Inetrstitials. Eg. 0 0 0  0.5 0.5 0.5.")
    parser.add_argument("-a","--antisite", dest="is_antisite", action="store_false",
                        help="Set if antisites are considered.")
    parser.add_argument("-e","--ElNeg_diff", dest="ElNeg_diff", type=float, default=1.0,
                        help="Criterion of the electronegativity difference for constructing antisites and impurities.")
    parser.add_argument("-s","--symbreak", dest="symbreak", action="store_true",
                        help="Set if symmetry is not broken.")
    parser.add_argument("--symprec",dest="symprec", type=float, default=0.01,
                        help="Set the symprec [A].")

    opts = parser.parse_args()

    defect_in = DefectIn.from_str_file(
                poscar=opts.poscar, dopants=opts.dopants, interstitials=opts.interstitials, 
                is_antisite=opts.is_antisite, ElNeg_diff=opts.ElNeg_diff,irregular=None, 
                symbreak=opts.symbreak, symprec=opts.symprec)

    defect_in.print_defect_in() 

if __name__ == "__main__": main()
