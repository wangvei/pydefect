#!/usr/bin/env python

import sys
import warnings
import json
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import pydefect.input_generator.atom as atom
from monty.json import MontyEncoder

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

class IrrepElement():
    """
    >> IrreducibleSite()
    This class object holds properties related to irreducible (irrep) atom set.
    Note1: atomic indices need to be sorted. Thus, they can be written in one 
           sequence.
    Note2: first_index atom is assumed to represent the irreducible atoms.

    Args:
        irrepname (str): element name with irreducible index (e.g., Mg1)
        element (str): element name (e.g., Mg)
        first_index (int): first index of irrepname. 
        last_index (int): last index of irrepname. 
        repr_coords (array): representative coordinates, namely the position
                            of first_index

    TODO1: Add the site symmetry information.
    """
    def __init__(self, irrepname, element, first_index, last_index, 
                 repr_coords):
        self.irrepname = irrepname
        self.element = element
        self.first_index = first_index
        self.last_index = last_index
        self.repr_coords = repr_coords

    def __eq__(self, other):
        if other is None or type(self) != type(other): 
            raise TypeError
        return self.__dict__ == other.__dict__

    def as_dict(self):
        d = {"irrepname": self.irrepname,
             "element" : self.element,
             "first_index" : self.first_index,
             "last_index" : self.last_index,
             "repr_coords" : self.repr_coords}
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["irrepname"], d["element"], d["first_index"], 
                   d["last_index"], d["repr_coords"]) 

    @property
    def natoms(self):
        """
        Return number of atoms in a given (super)cell.
        """
        return self.last_index - self.first_index + 1


class DefectSetting():
    """
    This class object holds full information on the setting of the point 
    defect calculations.

    Args:
        structure: pmg Structure/IStructure class object
        irrep_elements: IrrepElement class object
        dopant_configs (array): dopant configurations,
            e.g., ["Al_Mg", "N_O"].
        antisite_configs (array): antisite configuraions, 
            e.g., ["Mg_O", "O_Mg"].
        interstitial_coords (3x1 array): coordinations of interstitial sites,
            e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ...].
        include (array): exceptionally added defect type with a charge state.
            e.g., ["Va_O1_-1", "Va_O1_-2"]                     
        exclude (array): exceptionally removed defect type with a charge state.
                         In case some of them don't exist, they'll be ignored.
            e.g., ["Va_O1_1", "Va_O1_2"]                     
        is_sym_broken (bool): Whether to perturb defect's neighboring atoms for 
#        symbreak (bool): Whether to perturb defect's neighboring atoms for 
                         symmetry breaking.
        displace (float): Maximum displacement distance in angstrom.
        cutoff (float): Cutoff radius in which atoms are displaced.
        symprec (float): Precision used for symmetry analysis.
        oxidation_states (dict): Oxidation states for relevant elements.
        electronegativity (dict): Electronegativity for relevant elements.
    """
    def __init__(self, structure, irrep_elements, dopant_configs, 
                 antisite_configs, interstitial_coords, include, exclude, 
                 symbreak, displace, cutoff, symprec, oxidation_states, 
                 electronegativity):

        self.structure = structure
        self.irrep_elements = irrep_elements
        self.dopant_configs = dopant_configs
        self.antisite_configs = antisite_configs
        self.interstitial_coords = interstitial_coords
        self.include = include
        self.exclude = exclude
        self.symbreak = symbreak
        self.displace = displace
        self.cutoff = cutoff
        self.symprec = symprec
        self.oxidation_states = oxidation_states
        self.electronegativity = electronegativity

    def __eq__(self, other):
        if other is None or type(self) != type(other): 
            raise TypeError
        return self.__dict__ == other.__dict__
#        return self.as_dict() == other.as_dict()

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        """
        d = {"structure": self.structure,
             "irrep_elements": [i.as_dict() for i in self.irrep_elements],
             "dopant_configs": self.dopant_configs,
             "antisite_configs": self.antisite_configs,
             "interstitial_coords": self.interstitial_coords,
             "include": self.include,
             "exclude": self.exclude,
             "symbreak": self.symbreak,
             "displace": self.displace,
             "cutoff": self.cutoff,
             "symprec": self.symprec,
             "oxidation_states": self.oxidation_states,
             "electronegativity": self.electronegativity}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Construct a DefectSetting class object from a dictionary.
        """
        irrep_elements = []
        for i in d["irrep_elements"]:
            irrep_elements.append(IrrepElement.from_dict(i))

        return cls(d["structure"], irrep_elements, d["dopant_configs"],
                   d["interstitial_coords"], d["antisite_configs"], 
                   d["include"], d["exclude"], d["symbreak"], d["displace"], 
                   d["cutoff"], d["symprec"], d["oxidation_states"], 
                   d["electronegativity"])

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        fw = open(filename, 'w')
        json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    #def from_json():

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """
        Construct DefectSetting class object from a defect.in file.
        Currently, format of the defect.in file is not flexible, 
        so be careful for manipulating it.   
        """
        structure = Structure.from_file(poscar)
        irrep_elements = []
        electronegativity = {}
        oxidation_states = {}
        dopants = []
        dopant_configs = []
        displace = None
        cutoff = None
        symprec = None

        with open(defect_in_file) as defect_in:
            for l in defect_in:
                line = l.split()
     
                if line == []: 
                    continue
                elif line[0] == "Name:":
                    irrepname = line[1]
                    # remove irreducible index number from irrepname
                    element = ''.join(
                                     [i for i in irrepname if not i.isdigit()])
                    # Skip one line which shows the representative index.
                    # First_index atom is assumed to represent irrepname atoms.
                    defect_in.readline()
                    # Representative atom index
                    first_index, last_index = [int(i) 
                          for i in defect_in.readline().split()[1].split("..")]
                    repr_coords = \
                           [float(i) for i in defect_in.readline().split()[1:]]
                    irrep_elements.append(IrrepElement(irrepname, element, 
                                    first_index, last_index, repr_coords))
    
                    electronegativity[element] = \
                                         float(defect_in.readline().split()[1])
                    oxidation_states[element] = \
                                         int(defect_in.readline().split()[1])
                elif line[0] == "Dopant:":
                    for d in line[1:]:
                        electronegativity[d] = \
                                         float(defect_in.readline().split()[1])
                        oxidation_states[d] = \
                                         int(defect_in.readline().split()[1])
                elif line[0] == "Dopant_site:":
                    dopant_configs = line[1:]
                elif line[0] == "Int_site:":
                    b = [float(''.join(i for i in line[i] if i.isdigit() or 
                                       i == '.')) for i in range(1, len(line))]
                    interstitial_coords = [b[i:i + 3] 
                                                  for i in range(0, len(b), 3)]
                elif line[0] == "Antisite:": 
                    antisite_configs = line[1:]
                    if antisite_configs is []: antisite_configs = False
                elif line[0] == "Symbreak:": 
                    if line[1] == "False":
                        symbreak = False
                    else: 
                        try:
                            displace = float(line[1])
                            symbreak = True
                        except:
                            raise NotSupportedFlagError(
                                           "Symbreak flag is not appropriate.")
                elif line[0] == "Cutoff:":
                    cutoff = float(line[1])
                elif line[0] == "Symprec:":
                    symprec = float(line[1])
                elif line[0] == "Include:":
                    include = line[1:]
                elif line[0] == "Exclude:":
                    exclude = line[1:]
                else:
                    raise NotSupportedFlagError(
                                           line[0] + " flag is not supported!")

        return cls(structure, irrep_elements, dopant_configs, antisite_configs, 
                   interstitial_coords, include, exclude, symbreak, displace, 
                   cutoff, symprec, oxidation_states, electronegativity)


class DefectInMaker():
    """
    This class generates defect.in file, which shows the setting of point 
    defect calculations.

    Args:
        dopants (array): dopant names, e.g., ["Al", "N"]
        interstitial_coords (3x1 array): coordinations of interstitial sites,
            e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ...].
        is_antisite (bool): Whether to consider antisite defects.
        ElNeg_diff (float): Electronegativity difference for determining sets 
                            of antisites and dopant sites. 
        include (array): exceptionally added defect type with a charge state.
            e.g., ["Va_O1_-1", "Va_O1_-2"]                     
        exclude (array): exceptionally removed defect type with a charge state.
            e.g., ["Va_O1_1", "Va_O1_2"]                     
        symbreak (bool): Whether to perturb defect's neighboring atoms for 
            symmetry breaking.
        displace (float): Maximum displacement distance in angstrom.
        cutoff (float): Cutoff radius for detemining atoms displaced.
        symprec (float): Precision used for symmetry analysis.
    """

    def __init__(self, structure, dopants, interstitial_coords, is_antisite, 
                 ElNeg_diff, include="", exclude="", symbreak=False, 
                 displace=0.2, cutoff=3.0, symprec=0.01):

        self.dopants = dopants
        if not len(interstitial_coords) % 3 == 0:
                warnings.warn("Be careful. Interstital site is not proper.")
        self.interstitial_coords = interstitial_coords
        self.is_antisite = is_antisite
        self.ElNeg_diff = ElNeg_diff
        self.include = include
        self.exclude = exclude
        self.symbreak = symbreak
        self.displace = displace
        self.cutoff = cutoff
        self.symprec = symprec
        self.electronegativity = {}
        self.oxidation_states = {}

        # Get electronegativity and oxidation states.
        # intrinsic elements + dopants
        for s in structure.symbol_set + tuple(self.dopants):
            try:
                self.electronegativity[s] = atom.electronegativity[s]
            except:
                warnings.warn("Electronegativity of " + s + " is unavailable.")
                self.electronegativity[s] = "N.A."
            try:
                self.oxidation_states[s] = atom.charge[s]
            # To use pmg oxidation states, uncomment the following
#               self.oxidation_states[s] = \
#                                        Element(s).common_oxidation_states[-1]
            except:
                warnings.warn("Oxidation state of " + s + " is unavailable.")
                self.oxidation_states[s] = "N.A."

        self.symmetrized_structure = \
                      SpacegroupAnalyzer(structure).get_symmetrized_structure()
        # num_irrep_elements["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irrep_elements = {}
        # irrep_elements (aray): a set of IrrepElement class objects
        self.irrep_elements = []
        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        last = 0
        equiv_sites = self.symmetrized_structure.equivalent_sites
        for i, e in enumerate(equiv_sites):
            element = e[0].species_string
            if element not in num_irrep_elements.keys():
                num_irrep_elements[element] = 1
            else:
                # increment number of inequiv site for element
                num_irrep_elements[element] += 1
            first = last + 1
            last = last + len(e)
            repr_coords = e[0].frac_coords
            irrepname = element + str(num_irrep_elements[element]) 
            self.irrep_elements.append(
                     IrrepElement(irrepname, element, first, last, repr_coords))

        EN_keys = self.electronegativity.keys()

        # E.g., antisite_configs = [["Mg, "O"], ...]
        self.antisite_configs = []
        if is_antisite is True:
            for s1 in structure.symbol_set:
                for s2 in structure.symbol_set:
                    if s1 == s2: 
                        continue
                    if s1 in EN_keys and s2 in EN_keys:
                        if abs(self.electronegativity[s1] - 
                               self.electronegativity[s2]) < ElNeg_diff:
                            self.antisite_configs.append([s1, s2])
                    else:
                        self.electronegativity_not_defined(s1, s2)

        # E.g., dopant_configs = [["Al", "Mg"], ...]
        self.dopant_configs = []
        if dopants:
            for d in dopants:
                if d in structure.symbol_set:
                    warnings.warn("Dopant " + d + " exists in host.")
                    continue
                for s1 in structure.symbol_set:
                    if s1 in EN_keys and d in EN_keys:
                        if abs(self.electronegativity[s1] - 
                               self.electronegativity[d]) < ElNeg_diff:
                            self.dopant_configs.append([d, s1])
                    else:
                        self.electronegativity_not_defined(d, s1)

        self.setting = DefectSetting(
                        self.symmetrized_structure, self.irrep_elements, 
                        self.dopant_configs, self.antisite_configs, 
                        self.interstitial_coords, self.include, self.exclude, 
                        self.symbreak, self.displace, self.cutoff, 
                        self.symprec, self.oxidation_states, 
                        self.electronegativity)

    def electronegativity_not_defined(self, element1, element2):
        print("Electronegativity of {} and/or {} is not defined".\
                                            format(element1, element2))

    @classmethod
    def from_str_file(cls, poscar, dopants=[], interstitial_coords=False,
         is_antisite=False, ElNeg_diff=1.0, include=None,
         exclude=None, symbreak=True, displace=0.2, cutoff=3.0, symprec=0.01):
        """
        Construct DefectInMaker class object from a POSCAR file.
        VERY IMPORTANT: Some parameters are set by default.
        """
        structure = Structure.from_file(poscar)

        return cls(structure, dopants, interstitial_coords, is_antisite, 
                   ElNeg_diff, include, exclude, symbreak, displace, cutoff, 
                   symprec)

    def to(self, defectin_file="defect.in", poscar_file="DPOSCAR"):
        """
        Print readable defect.in file.
        """
        self._write_defect_in(defectin_file)
        # HACK:  pmg has a bug, Symmetrized structure object cannot be poscar
        Structure.from_str(self.symmetrized_structure.to(fmt="cif"),
                                fmt="cif").to(fmt="poscar", filename="DPOSCAR")

    def _write_defect_in(self, defectin_file="defect.in"):
        with open(defectin_file, 'w') as fw:

            for e in self.irrep_elements:
                fw.write("  Name: {}\n".format(e.irrepname))
                fw.write("   Rep: {}\n".format(e.first_index))
                fw.write(" Equiv: {}\n".format(
                              str(e.first_index) + ".." + str(e.last_index)))
                fw.write(" Coord: %9.7f %9.7f %9.7f\n" % tuple(e.repr_coords))
                fw.write("EleNeg: {}\n".format(
                                           self.electronegativity[e.element]))
                fw.write("Charge: {}\n\n".format(
                                             self.oxidation_states[e.element]))
            fw.write("Int_site: ")

            if self.interstitial_coords:
                if type(self.interstitial_coords) == str:
                    coords = [float(i) 
                                     for i in self.interstitial_coords.split()]
                    if not len(coords) % 3 == 0:
                        raise ValueError(
                                "The interstitial coordinates are not proper.")
                else:
                    coords = self.interstitial_coords
                fw.write(str([coords[i:i + 3] 
                                     for i in range(0, len(coords), 3)]) +"\n")
            else: fw.write("\n")

            if self.antisite_configs is not []:
                fw.write("Antisite: ")
                fw.write(' '.join(i[0] + "_" + i[1] 
                                       for i in self.antisite_configs)+ "\n\n")
            if self.dopants:
                for d in self.dopants:
                    fw.write("Dopant: {}\n".format(d))
                    fw.write("EleNeg: {}\n".format(self.electronegativity[d]))
                    fw.write("Charge: {}\n".format(self.oxidation_states[d]))
                    fw.write("\n")
    
                fw.write("Dopant_site: ")
                fw.write(' '.join(i[0] + "_" + i[1] 
                                        for i in self.dopant_configs) + "\n\n")
            if self.symbreak == True:
                fw.write("Symbreak: {}\n".format(self.displace))
            else:
                fw.write("Symbreak: {}\n".format(self.symbreak))
            fw.write("Include: {}\n".format(self.include))
            fw.write("Exclude: {}\n".format(self.exclude))
            fw.write("Cutoff: {}\n".format(self.cutoff))
            fw.write("Symprec: {}\n".format(self.symprec))

    @staticmethod
    def print_dopant_info(dopant):
        """
        This is used for adding dopant information a posteriori.
        """
        try:
            electronegativity = atom.electronegativity[dopant]
        except:
            warnings.warn("Electronegativity of " + s + " is unavailable.")
            electronegativity[dopant] = "N.A."
        try:
            oxidation_states = atom.charge[dopant]
        except:
            warnings.warn("Oxidation state of " + s + " is unavailable.")
            oxidation_states = "N.A."

        print("Dopant: {}".format(dopant))
        print("EleNeg: {}".format(electronegativity))
        print("Charge: {}".format(oxidation_states))


class NotSupportedFlagError(Exception):
    pass 


def main():
    import argparse
    parser = argparse.ArgumentParser(
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name.")
    parser.add_argument("-d","--dopants", dest="dopants", default="", nargs="+", 
                        type=str, help="Dopant elements. Eg. Al Ga In.")
    parser.add_argument("-i", dest="interstitial_coords", default=False, 
                        nargs="+", type=float, 
                        help="Inetrstitials. Eg. 0 0 0  0.5 0.5 0.5.")
    parser.add_argument("-a","--antisite", dest="is_antisite", 
                        action="store_false",
                        help="Set if antisites are considered.")
    parser.add_argument("-e","--ElNeg_diff", dest="ElNeg_diff", type=float, 
                        default=1.0, help="Criterion of the electronegativity \
                        difference for constructing antisite and impurities.")
    parser.add_argument("--include", dest="include", type=str, default="",
                        help="Exceptionally included defects. E.g. Va_O2_-1.")
    parser.add_argument("--exclude", dest="exclude", type=str, default="",
                        help="Exceptionally excluded defects. E.g. Va_O2_0.")
    parser.add_argument("-s","--symbreak", dest="symbreak", action="store_true",
                        help="Set if symmetry is not broken.")
    parser.add_argument("--displace", dest="displace", type=float, default=0.2,
                        help="Displacement distance.")
    parser.add_argument("--cutoff",dest="cutoff", type=float, default=3.0,
                        help="Set the cutoff [A].")
    parser.add_argument("--symprec",dest="symprec", type=float, default=0.01,
                        help="Set the symprec [A].")
    parser.add_argument("--print_dopant", dest="print_dopant", default=None,
                        type=str, help="Print Dopant information.")

    opts = parser.parse_args()

    if opts.print_dopant:
        DefectInMaker.print_dopant_info(opts.print_dopant)
    else:
        defect_in = DefectInMaker.from_str_file(opts.poscar, opts.dopants, 
                    opts.interstitial_coords, opts.is_antisite, 
                    opts.ElNeg_diff, opts.include, opts.exclude, opts.symbreak,
                    opts.displace, opts.cutoff, opts.symprec)
        defect_in.to()

if __name__ == "__main__": main()
