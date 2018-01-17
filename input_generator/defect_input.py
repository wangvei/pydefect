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
    This class object has some properties related to irreducible atoms.
    In the structure, the atomic indices need to be sorted.

    **CAUTION**: first_index atomic is assumed to represent irreducible atoms.

    Args:
        irrepname (str): element name with irreducible index (e.g., Mg1)
        element (str): element name (e.g., Mg)
        first_index (int): first index of *irrepname. 
        last_index (int): last index of *irrepname. 
        repr_coord (array): representative coordination

    TODO1: Add the site symmetry information.
    """
    def __init__(self, irrepname, element, first_index, last_index, 
                 repr_coord):
        self.irrepname = irrepname
        self.element = element
        self.first_index = first_index
        self.last_index = last_index
        self.repr_coord = repr_coord

    def __eq__(self, other):
        if other is None or type(self) != type(other): 
            raise TypeError
        return self.__dict__ == other.__dict__

    def as_dict(self):
        d = {"irrepname": self.irrepname,
             "element" : self.element,
             "first_index" : self.first_index,
             "last_index" : self.last_index,
             "repr_coord" : self.repr_coord}
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["irrepname"], d["element"], d["first_index"], 
                   d["last_index"], d["repr_coord"]) 

    @property
    def multiplicity(self):
        return self.last_index - self.first_index + 1


class DefectSetting():
    """
    This class object holds all the information on the setting of the point 
    defect calculations.

    Args:
        structure: pmg Structure/IStructure class object
        irrep_elements: IrrepElement class object
        dopant_configs (array): dopant configurations,
            e.g., ["Al_Mg", "N_O"].
        antisite_configs (2x1 array): antisite configuraions, 
            e.g., [["Mg", "O"], ["O", "Mg"]].
        interstitial_coords (3x1 array): coordinations of interstitial sites,
            e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ...].
        include (array): exceptionally added defect type with a charge state.
            e.g., ["Va_O_-1", "Va_O_-2"]                     
        exclude (array): exceptionally removed defect type with a charge state.
            e.g., ["Va_O_1", "Va_O_2"]                     
        symbreak (bool): Whether to perturb defect's neighboring atoms for 
            symmetry breaking.
        displace (float): Maxmum displacement distance in angstrom.
        cutoff (float): Cutoff radius for detemining atoms displaced.
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
#        print(self.__dict__, other.__dict__)
        return self.__dict__ == other.__dict__

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

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def to_json_file(self, filename):
        """
        Returns a json string representation of the MSONable object.
        """
        fw = open(filename, 'w')
        json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def from_dict(cls, d):
        """ NEED TO BE MODIFIED """
        irrep_elements = []
        for i in d["irrep_elements"]:
            irrep_elements.append(IrrepElement.from_dict(i))

        return cls(d["structure"], irrep_elements, d["dopant_configs"],
                   d["interstitial_coords"], d["antisite_configs"], 
                   d["include"], d["exclude"], d["symbreak"], d["displace"], 
                   d["cutoff"], d["symprec"], d["oxidation_states"], 
                   d["electronegativity"])

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """
        Construct DefectSetting class object by parsing a defect.in file.
        Currently, format of the defect.in file is inflexible, so be careful
        for manipulating it.   
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
                    repr_coord = \
                           [float(i) for i in defect_in.readline().split()[1:]]
                    irrep_elements.append(IrrepElement(irrepname, element, 
                                    first_index, last_index, repr_coord))
    
                    electronegativity[element] = \
                                         float(defect_in.readline().split()[1])
                    oxidation_states[element] = \
                                         float(defect_in.readline().split()[1])
                elif line[0] == "Dopant:":
                    for d in line[1:]:
                        electronegativity[d] = \
                                         float(defect_in.readline().split()[1])
                        oxidation_states[d] = \
                                         float(defect_in.readline().split()[1])
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
                            displace = float(symbreak)
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


class NotSupportedFlagError(Exception):
    pass 


class DefectInMaker():
    """
    This class generates defect.in file, which controles the setting of 
    point defect calculations.

    Args:
        dopants (array): dopant names, e.g., ["Al", "N"]
        interstitial_coords (3x1 array): coordinations of interstitial sites,
            e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ...].
        is_antisite (bool): Whether to consider antisite defects.
        ElNeg_diff (float): Electronegativity difference used for determining
            sets of antisites and dopant sites. 
        include (array): exceptionally added defect type with a charge state.
            e.g., ["Va_O_-1", "Va_O_-2"]                     
        exclude (array): exceptionally removed defect type with a charge state.
            e.g., ["Va_O_1", "Va_O_2"]                     
        symbreak (bool): Whether to perturb defect's neighboring atoms for 
            symmetry breaking.
        displace (float): Maxmum displacement distance in angstrom.
        cutoff (float): Cutoff radius for detemining atoms displaced.
        symprec (float): Precision used for symmetry analysis.
    """

    def __init__(self, structure, dopants, interstitial_coords, is_antisite, 
                 ElNeg_diff, include="", exclude="", symbreak=False, 
                 displace=0.2, cutoff=3.0, symprec=0.01):

        self.dopants = dopants
        self.interstitial_coords = interstitial_coords
        self.is_antisite = is_antisite
        self.ElNeg_diff = ElNeg_diff
        self.include = include
        self.exclude = exclude
        self.symbreak = symbreak
        self.displace = displace
        self.cutoff = cutoff
        self.symprec = symprec
        # only vasp5 POSCAR format is supported for getting the elements names.
        # equivalent_positions: Equivalent positions from SpacegroupAnalyzer.
        self.equiv_site_seq = \
     SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites

        # irrep_elements_index["Mg"] = 2 means that Mg element has two 
        # inequivalent sites
        self._irrep_element_index = {}
        self.electronegativity = {}
        self.oxidation_states = {}

        # Get electron negativity and oxidation states.
        # intrinsic elements + dopants
        for s in structure.symbol_set + tuple(self.dopants):
            try:
                self.electronegativity[s] = atom.electronegativity[s]
            except:
                warnings.warn("Electronegativity of " + s + " is unavailable.")
                self.electronegativity[s] = "N.A."
            try:
                self.oxidation_states[s] = atom.charge[s]
            # if one wants to use pmg oxidation states. Then, needs to import 
            # Element class
#                oxidation_states[s] = Element(s).common_oxidation_states[-1]
            except:
                warnings.warn("Oxidation state of " + s + " is unavailable.")
                self.oxidation_states[s] = "N.A."

        # irrep_elements (aray): a set of IrrepElement class objects
        self.irrep_elements = []
        # list of element names for all atoms, e.g., ["Mg", "Mg", "O", "O"]
        self._elements = []
        # an array of frac_coords for all atoms, e.g., [[0, 0, 0], ...]
        self._frac_coords = []

        for inequiv_site in self.equiv_site_seq:
            self._element = inequiv_site[0].species_string
            if self._element not in self._irrep_element_index.keys():
                self._irrep_element_index[self._element] = 1
            else:
                # increment number of inequiv site
                self._irrep_element_index[self._element] += 1
            self._first = len(self._elements) + 1
            self._last = len(self._elements) + len(inequiv_site)
            self._repr_coord = inequiv_site[0].frac_coords
            # E.g., IrrepElement("Mg1", "Mg", 1, 32, [0, 0, 0])
            self.irrep_elements.append(IrrepElement(self._element 
                               + str(self._irrep_element_index[self._element]), 
                                 self._element, self._first, self._last, 
                                 self._repr_coord))
#                                 self._repr_coord).as_dict())
            for equiv_site in inequiv_site:
                self._elements.append(equiv_site.species_string)
                self._frac_coords.append(equiv_site.frac_coords)

        # Atoms are sorted by symmetry, which will be written in DPOSCAR
        self.structure = Structure(structure.lattice, 
                                             self._elements, self._frac_coords)
        # E.g., antisite_configs = [["Mg", "O1"], ...]
        self.antisite_configs = []
        if is_antisite is True:
            for i in self.irrep_elements:
                for s in structure.symbol_set:
                    if i.element == s: 
                        continue
                    try:
                        ENdiff = electronegativity[i.element] \
                               - electronegativity[s]
                    except:
                        # A fictitious number is inserted for atoms w/o 
                        # electronegativity.
                        ENdiff = 100
                    if abs(ENdiff) < ElNeg_diff:
                        self.antisite_configs.append([s, i.element])

        # E.g., dopant_configs = [["Al", "Mg"], ...]
        self.dopant_configs = []
        if dopants:
            for i in self.irrep_elements:
                for d in dopants:
                    try:
                        ENdiff = electronegativity[i.element] \
                               - electronegativity[d]
                    except:
                        # A fictitious number is inserted for atoms w/o 
                        # electronegativity.
                        ENdiff = 100
                    if abs(ENdiff) < ElNeg_diff:
                        self.dopant_configs.append([d, i.element])

        self.setting = DefectSetting(self.structure, self.irrep_elements, 
                        self.dopant_configs, self.antisite_configs, 
                        self.interstitial_coords, self.include, self.exclude, 
                        self.symbreak, self.displace, self.cutoff, 
                        self.symprec, self.oxidation_states, 
                        self.electronegativity)

    @classmethod
    def from_str_file(cls, poscar, dopants=[], interstitial_coords=False,
         is_antisite=False, ElNeg_diff=1.0, include=None,
         exclude=None, symbreak=True, displace=0.2, cutoff=3.0, symprec=0.01):
        """
        Constructs DefectInMaker class object from a POSCAR file.
        VERY IMPORTANT: Some parameters are set by default.
        """
        structure = Structure.from_file(poscar)

        return cls(structure, dopants, interstitial_coords, is_antisite, 
                   ElNeg_diff, include, exclude, symbreak, displace, cutoff, 
                   symprec)

    def to(self, defectin_file="defect.in", poscar_file="DPOSCAR"):
        """
        Prints readable defect.in file.
        """
        self._write_defect_in(defectin_file)
        self.structure.to(fmt="poscar", filename=poscar_file)

    def _write_defect_in(self, defectin_file="defect.in"):
        with open(defectin_file, 'w') as f:

            for e in self.irrep_elements:
                f.write("  Name: {}\n".format(e.irrepname))
                f.write("   Rep: {}\n".format(e.first_index))
                f.write(" Equiv: {}\n".format(
                              str(e.first_index) + ".." + str(e.last_index)))
                f.write(" Coord: %9.7f %9.7f %9.7f\n" % tuple(e.repr_coord))
                f.write("EleNeg: {}\n".format(
                                           self.electronegativity[e.element]))
                f.write("Charge: {}\n\n".format(
                                             self.oxidation_states[e.element]))
            f.write("Int_site: ")
            if self.interstitial_coords:
                if type(self.interstitial_coords) == str:
                    coords = [float(i) 
                                     for i in self.interstitial_coords.split()]
                    if not len(coords) % 3 == 0:
                        raise ValueError(
                                "The interstitial coordinates are not proper.")
                else:
                    coords = self.interstitial_coords
                f.write(str([coords[i:i + 3] 
                                     for i in range(0, len(coords), 3)]) +"\n")
            if self.antisite_configs is not []:
                f.write("Antisite: ")
                f.write(' '.join(i[0] + "_" + i[1] 
                                       for i in self.antisite_configs)+ "\n\n")
            if self.dopants:
                for d in self.dopants:
                    f.write("Dopant: {}\n".format(d))
                    f.write("EleNeg: {}\n".format(self.electronegativity[d]))
                    f.write("Charge: {}\n".format(self.oxidation_states[d]))
                    f.write("\n")
    
                f.write("Dopant_site: ")
                f.write(' '.join(i[0] + "_" + i[1] 
                                        for i in self.dopant_configs) + "\n\n")
            if self.symbreak == True:
                f.write("Symbreak: {}\n".format(self.displace))
            else:
                f.write("Symbreak: {}\n".format(self.symbreak))
            f.write("Include: {}\n".format(self.include))
            f.write("Exclude: {}\n".format(self.exclude))
            f.write("Symprec: {}\n".format(self.symprec))

    @staticmethod
    def print_dopant_info(dopant):
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

def main():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name.")
    parser.add_argument("-d","--dopants", dest="dopants", default="", nargs="+", 
                        type=str, help="Dopant elements. Eg. Al Ga In.")
    parser.add_argument("-i", dest="interstitial_coords", default=False, nargs="+",
                        type=float, help="Inetrstitials. Eg. 0 0 0  0.5 0.5 0.5.")
    parser.add_argument("-a","--antisite", dest="is_antisite", action="store_false",
                        help="Set if antisites are considered.")
    parser.add_argument("-e","--ElNeg_diff", dest="ElNeg_diff", type=float, default=1.0,
                        help="Criterion of the electronegativity difference for constructing antisite and impurities.")
    parser.add_argument("--include", dest="include", type=str, default="",
                        help="Exceptionally included defect type. Eg Va_O2_-1.")
    parser.add_argument("--exclude", dest="exclude", type=str, default="",
                        help="Exceptionally excluded defect type. Eg Va_O2_0.")
    parser.add_argument("-s","--symbreak", dest="symbreak", action="store_true",
                        help="Set if symmetry is not broken.")
    parser.add_argument("--displace", dest="displace", type=float, default=0.2,
                        help="Displacement distance.")
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
                    opts.displace, opts.symprec)
        defect_in.to()

if __name__ == "__main__": main()
