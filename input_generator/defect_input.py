#!/usr/bin/env python

#import itertools as it
#import numpy as np
import sys
import warnings
import json
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import pydefect.input_generator.atom as atom
from monty.json import MontyEncoder

# from pymatgen.core.periodic_table import Element, Specie, get_el_sp

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

class IrrepElement():
    """
    This class object has various properties related to atoms.
    """
    def __init__(self, irrepname, element, multiplicity, repr_coord):
        self.irrepname = irrepname
        self.element = element
        self.multiplicity = multiplicity
        self.repr_coord = repr_coord

    def as_dict(self):
        d = {"irrepname": self.irrepname,
             "element" : self.element,
             "multiplicity" : self.multiplicity,
             "repr_coord" : self.repr_coord}
        return d


class DefectSetting():
    def __init__(self, structure, irrep_elements, dopant_sites, interstitials,
                 antisites, include, exclude, symbreak, symprec, 
                 oxidation_states, electron_negativity):

        self.structure = structure
        self.irrep_elements = irrep_elements
        self.dopant_sites = dopant_sites
        self.interstitials = interstitials
        self.antisites = antisites
        self.include = include
        self.exclude = exclude
        self.symbreak = symbreak
        self.symprec = symprec
        self.oxidation_states = oxidation_states
        self.electron_negativity = electron_negativity

    def as_dict(self):
        """
        Dict representation of DefectIn class object.
        """
        d = {"structure": self.structure,
             "irrep_elements": self.irrep_elements,
             "dopant_sites": self.dopant_sites,
             "interstitials": self.interstitials,
             "antisites": self.antisites,
             "include": self.include,
             "exclude": self.exclude,
             "symbreak": self.symbreak,
             "symprec": self.symprec,
             "oxidation_states": self.oxidation_states,
             "electron_negativity": self.electron_negativity}
        return d

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    @classmethod
    def from_dict(cls, d):
        """ NEED TO BE MODIFIED """
        return cls(d["structure"], d["irrep_elements"], d["dopant_sites"],
                   d["interstitials"], d["antisites"], d["include"], 
                   d["exclude"], d["symbreak"], d["symprec"], 
                   d["oxidation_states"], d["electron_negativity"])

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """
        Construct DefectSetting class object from a defect.in file.
        """
 
        structure = Structure.from_file(poscar)
        defect_in = open(defect_in_file)
        irrep_elements = []
        electron_negativity = {}
        oxidation_states = {}
        dopant_sites = []
        interstitials = [] # [[i1(float)], [i2(float)], .. ]

        while True:
            line = defect_in.readline().split()
 
            if line == []: continue
            elif line[0] == "Name:":
                irrepname = line[1]
                # remove number from irrepname
                element = ''.join([i for i in irrepname if not i.isdigit()])
                # Representative atom index
                first, last = defect_in.readline().split()[1].split("..")
                multiplicity = int(last) - int(first)
                repr_coord = \
                           [float(i) for i in defect_in.readline().split()[1:]]
                irrep_elements.append(IrrepElement(irrepname, element, multiplicity, repr_coord).as_dict())
                electron_negativity[element] = \
                                    float(defect_in.readline().split()[1])
                oxidation_states[element] = \
                                    float(defect_in.readline().split()[1])

            elif line[0] == "Dopant:":
                dopant_element = dopants.append(line[1])
                electron_negativity[dopant_element] = \
                                        float(defect_in.readline().split()[1])
                oxidation_states[dopant_element] = \
                                        float(defect_in.readline().split()[1])
 
            elif line[0] == "Int_site:":
                b = [_get_num(line[i]) for i in range(1, len(line))]
                interstitials = [b[i:i + 3] for i in range(0, len(b), 3)]

            elif line[0] == "Antisite:": 
                antisites = line[1:]
                if anti_site is []: anti_site = False

            elif line[0] == "Symbreak:": 
                symbreak = line[1]
                if symbreak is not False:
                    symprec = float(symbreak)
                    symbreak = True

            elif line[0] == "Include:":
                include = line[1:]

            elif line[0] == "Exclude:":
                exclude = line[1:]

            else:
                raise NotSupportedFlagError

        return cls(structure, irrep_elements, dopant_sites, interstitials,
                   antisites, include, exclude, symbreak, symprec, 
                   oxidation_states, electron_negativity)

class NotSupportedFlagError(Exception):
    pass 

class DefectIn():
    """
    It generates and parses defect.in file, which controles the input setting 
    of defect calculations.
    """

#    def __init__(self, structure, dopants=[], interstitials=False,
#                 is_antisite=False, ElNeg_diff=1.0, include=None, 
#                 exclude=None, symbreak=True, symprec=1e-5):
    def __init__(self, structure, dopants, interstitials, is_antisite, 
                 ElNeg_diff, include, exclude, symbreak, symprec):

        self.dopants = dopants
        self.interstitials = interstitials
        self.is_antisite = is_antisite
        self.ElNeg_diff = ElNeg_diff
        self.include = include
        self.exclude = exclude
        self.symbreak = symbreak
        self.symprec = symprec
        # only vasp5 POSCAR format is supported.
        self.equiv_site_seq = SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites

        # e.g.,  irrep_elements_index["Mg"] = 2 means Mg element has two inequivalent sites
        self._irrep_element_index = {}
        self.electron_negativity = {}
        self.oxidation_states = {}

        # Obtain electron negativity and oxidation states.
        # intrinsic elements + dopants
        for s in structure.symbol_set + tuple(self.dopants):
            self._irrep_element_index[s] = 0
            try:
                self.electron_negativity[s] = atom.electron_negativity[s]
            except:
                warnings.warn("Electron negativity of " + s + " is unavailable.")
                self.electron_negativity[s] = "N.A."
            try:
                self.oxidation_states[s] = atom.charge[s]
            # if one wants to use pmg oxidation states.
#                oxidation_states[s] = Element(s).common_oxidation_states[-1]
            except:
                warnings.warn("Oxidation state of " + s + " is unavailable.")
                self.oxidation_states[s] = "N.A."

#        self.electron_negativity = electron_negativity
#        self.oxidation_states = oxidation_states

        # e.g.,  irrep_elements[0] = "IrrepSpecie object"
        self.irrep_elements = []
        # list of elements, eg [Mg, Mg, O, O]
        self._elements = []
        # list of frac_coords, eg [[0,0,0],[0.5,0,0],...]
        self._frac_coords = []

        for inequiv_site in self.equiv_site_seq:
            self._element = inequiv_site[0].species_string
            self._irrep_element_index[self._element] += 1  # increment number of inequiv site
            self._multiplicity = len(inequiv_site)
            self._repr_coord = inequiv_site[0].frac_coords
            self.irrep_elements.append(IrrepElement(self._element + str(self._irrep_element_index[self._element]), 
                                     self._element, self._multiplicity, self._repr_coord).as_dict())
            for equiv_site in inequiv_site:
                self._elements.append(equiv_site.species_string)
                self._frac_coords.append(equiv_site.frac_coords)

        # Atoms are sorted by symmetry, which will be written in DPOSCAR
        self.structure = Structure(structure.lattice, self._elements, self._frac_coords)

        self.antisites = []
        if is_antisite is True:
            # eg. antisites = [[ "Mg", "O1"], ...]
            for i in self.irrep_elements:
                for s in structure.symbol_set:
                    if i["element"] == s: 
                        continue
                    try:
                        ENdiff = electron_negativity[i["element"]] - electron_negativity[s]
                    except:
                        # It's a fictitious number
                        ENdiff = 100
                    if abs(ENdiff) < ElNeg_diff:
                        self.antisites.append([s, i["element"]])

        self.dopant_sites = []
        if dopants:
            for i in self.irrep_elements:
                for d in dopants:
                    try:
                        ENdiff = electron_negativity[i["element"]] - electron_negativity[d]
                    except:
                        # It's a fictitious number
                        ENdiff = 100
                    if abs(ENdiff) < ElNeg_diff:
                        self.dopant_sites.append([d, i["element"]])

        self.setting = DefectSetting(self.structure, self.irrep_elements, self.dopant_sites, self.interstitials,              
                         self.antisites, self.include, self.exclude, self.symbreak, self.symprec, self.oxidation_states, self.electron_negativity)

    @classmethod
    def from_str_file(cls, poscar, dopants=[], interstitials=False,
                 is_antisite=False, ElNeg_diff=1.0, include=None, 
                 exclude=None, symbreak=True, symprec=1e-5):
        """
        Construct DefectIn class object from a POSCAR file.
        Some parameters are set by default.
        """

        structure = Structure.from_file(poscar)

        return cls(structure, dopants, interstitials, is_antisite, ElNeg_diff, 
                   include, exclude, symbreak, symprec)

    def to(self, filename1="defect.in", filename2="DPOSCAR"):
        """
        Print readable defect.in file.
        """
        self._print_defect_in(filename=filename1)
        self.structure.to(fmt="poscar", filename=filename2)

    def _print_defect_in(self, filename="defect.in"):
        file = open(filename, 'w')                                             

        i = 1
        #        print(self.repr_frac_coords)
        for e in self.irrep_elements:
            file.write("  Name: {}\n".format(e.irrepname))
            file.write("   Rep: {}\n".format(e.repr_coord))
            file.write(" Equiv: {}\n".format(str(i) + ".." + str(i + e.multiplicity - 1)))
            i += e.multiplicity
            file.write(" Coord: %9.7f %9.7f %9.7f\n" % tuple(e.repr_coord))
            file.write("EleNeg: {}\n".format(self.electron_negativity[e.element]))
            file.write("Charge: {}\n\n".format(self.oxidation_states[e.element]))

        if self.interstitials:
            coords = [float(i) for i in self.interstitials.split()]
            if not len(coords) % 3 == 0:
                raise ValueError("The interstitial coordinates are not proper.")
            file.write("Int_site: ")
            file.write(str([coords[i:i + 3] for i in range(0, len(coords), 3)]) +"\n")

        if self.antisites is not []:
            file.write("Antisite: ")
            file.write(' '.join(i[0] + "_" + i[1] for i in self.antisites)+ "\n\n")

        if self.dopants:
            for d in self.dopants:
                file.write("Dopant: {}\n".format(d))
                file.write("EleNeg: {}\n".format(self.electron_negativity[d]))
                file.write("Charge: {}\n".format(self.oxidation_states[d]))
            file.write("\n")
            file.write("Dopant_site: ")
            file.write(' '.join(i[0] + "_" + i[1] for i in self.dopant_sites) + "\n\n")

        if self.symbreak == True:
            file.write("Symbreak: {}\n".format(self.symprec))
        else:
            file.write("Symbreak: {}\n".format(self.symbreak))
            

        if self.include:
            file.write("Include: {}\n".format(self.include))
        else:
            file.write("Include: \n")

        if self.exclude:
            file.write("Exclude: {}\n".format(self.exclude))
        else:
            file.write("Exclude: \n")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name.")
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
                is_antisite=opts.is_antisite, ElNeg_diff=opts.ElNeg_diff, include=None, exclude=None,
                symbreak=opts.symbreak, symprec=opts.symprec)

    defect_in.print_defect_in() 

if __name__ == "__main__": main()
