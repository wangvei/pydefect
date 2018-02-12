#!/usr/bin/env python

import warnings
import json

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pydefect.input_generator import atom
from pydefect.input_generator.defect import IrreducibleSite
from monty.json import MontyEncoder
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

# Following defaults determine the condition of automatic defect calculations.
_EN_DIFF = 1.0
_DISTANCE = 0.2
_CUTOFF = 3.0
_SYMPREC = 0.01

class DefectSetting():
    """
    This class object holds full information on the setting of the point 
    defect calculations.

    Args:
        structure (Structure): pmg Structure/IStructure class object
        irreducible_sites: IrreducibleSite class objects
        dopant_configs (2x1 array): dopant configurations,
                                    e.g., [["Al", Mg"], ["N", "O"]
        antisite_configs (array): antisite configurations,
                                  e.g., [["Mg","O"], ["O", "Mg"]]
        interstitial_coords (3x1 array): coordinates of interstitial sites,
                                         e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ..]
        included (array): exceptionally added defect type with a charge state.
                          e.g., ["Va_O1_-1", "Va_O1_-2"]
        excluded (array): exceptionally removed defect type with a charge state.
                          In case some of them don't exist, they'll be ignored.
                          e.g., ["Va_O1_1", "Va_O1_2"]
        distance (float): Maximum displacement distance in angstrom. 0 means
                          random displacement is not considered.
        cutoff (float): Cutoff radius in which atoms are displaced.
        symprec (float): Precision used for symmetry analysis.
        oxidation_states (dict): Oxidation states for relevant elements.
        electronegativity (dict): Electronegativity for relevant elements.
    """

    def __init__(self, structure, irreducible_sites, dopant_configs,
                 antisite_configs, interstitial_coords, included, excluded,
                 distance, cutoff, symprec, oxidation_states,
                 electronegativity):

        self.structure = structure
        self.irreducible_sites = irreducible_sites
        self.dopant_configs = dopant_configs
        # dopant element names are derived from in_name of dopant_configs.
        self.dopants = set([d[0] for d in dopant_configs])
        self.antisite_configs = antisite_configs
        self.interstitial_coords = interstitial_coords
        self.included = included
        self.excluded = excluded
        self.distance = distance
        self.cutoff = cutoff
        self.symprec = symprec
        self.oxidation_states = oxidation_states
        self.electronegativity = electronegativity

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        # This cannot be "return self.__dict__ == other.__dict__",
        # because irreducible_sites shows just pointers.
        # We need to make dictionary for irreducible_sites.
        return self.as_dict() == other.as_dict()

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectSetting class object from a dictionary.
        """
        # Expansion of irreducible_sites is necessary,
        # because irreducible_sites shows just pointers.
        irreducible_sites = []
        for i in d["irreducible_sites"]:
            irreducible_sites.append(IrreducibleSite.from_dict(i))

        return cls(d["structure"], irreducible_sites, d["dopant_configs"],
                   d["antisite_configs"], d["interstitial_coords"],
                   d["included"], d["excluded"], d["distance"],
                   d["cutoff"], d["symprec"], d["oxidation_states"],
                   d["electronegativity"])

    @classmethod
    def json_load(cls, filename):
        """
        Construct a DefectSetting class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """
        Construct DefectSetting class object from a defect.in file.
        Currently, format of the defect.in file is not flexible,
        so be careful when modifying it.
        """
        structure = Structure.from_file(poscar)
        irreducible_sites = []
        electronegativity = {}
        oxidation_states = {}
        dopant_configs = []
        distance = None
        cutoff = None
        symprec = None

        with open(defect_in_file) as di:
            for l in di:
                line = l.split()

                if not line:
                    continue

                elif line[0] == "Irreducible":
                    irreducible_name = line[2]
                    # remove index from irreducible_name, e.g., "Mg1" --> "Mg"
                    element = ''.join(
                        [i for i in irreducible_name if not i.isdigit()])
                    first_index, last_index = [int(i)
                                               for i in
                                               di.readline().split()[2].split(
                                                   "..")]
                    repr_coords = [float(i) for i in di.readline().split()[2:]]
                    irreducible_sites.append(IrreducibleSite(irreducible_name,
                                                             element,
                                                             first_index,
                                                             last_index,
                                                             repr_coords))
                    electronegativity[element] = float(di.readline().split()[1])
                    oxidation_states[element] = int(di.readline().split()[2])

                elif line[0] == "Interstitial":
                    # Interstitial coordinates: [[0, 0, 0], [0.2, 0.2, 0.2]]
                    # coords_native = ["[0,", "0," "0],", "0.2", "0.2", "0.2]]"]
                    coords_native = [i for i in range(2, len(line))]
                    # coords_str = ["0", "0", "0", "0.2", "0.2", "0.2"]
                    coords_str = [(''.join(i for i in line[i]
                                           if i.isdigit() or i == '.'))
                                  for i in coords_native]
                    # coords_float = [0, 0, 0, 0.2, 0.2, 0.2]
                    coords_float = [float(i) for i in coords_str]
                    try:
                        # interstitial_coords = [[0, 0, 0], [0.2, 0.2, 0.2]]
                        interstitial_coords = [b[i:i + 3]
                                     for i in range(int(len(coords_float) / 3))]
                    except ValueError:
                        print("Interstitial coordinates isn't a multiple of 3.")

                elif line[0] == "Antisite":
                    antisite_configs = line[2:]

                elif line[0] == "Dopant":
                    d = line[2]
                    electronegativity[d] = float(di.readline().split()[1])
                    oxidation_states[d] = int(di.readline().split()[2])

                elif line[0] == "Substituted":
                    dopant_configs = [i.split("_") for i in  line[2:]]

                elif line[0] == "Maximum":
                    distance = float(line[2])

                elif line[0] == "Cutoff":
                    cutoff = float(line[5])

                elif line[0] == "Symprec:":
                    symprec = float(line[1])

                elif line[0] == "Exceptionally":
                    if line[1] == "included:":
                        included = line[2:]
                    elif line[1] == "excluded:":
                        excluded = line[2:]

                else:
                    raise NotSupportedFlagError(line[0] + " is not supported!")

        return cls(structure, irreducible_sites, dopant_configs,
                   antisite_configs, interstitial_coords, included, excluded,
                   distance, cutoff, symprec, oxidation_states,
                   electronegativity)

    @classmethod
    def from_basic_settings(cls, poscar, dopants=[], interstitial_coords=[],
                            is_antisite=True, EN_diff=_EN_DIFF, included="",
                            excluded="", distance=_DISTANCE, cutoff=_CUTOFF,
                            symprec=_SYMPREC):
        """
        Generates DefectSetting object with some default settings.

        Args:
            dopants (array): dopant element names, e.g., ["Al", "N"]
            is_antisite (bool): Whether to consider antisite defects.
            EN_diff (float): Electronegativity (EN) difference for determining
                             sets of antisites and dopant sites.
        """


class DefectInMaker:
    """
    This class generates DefectSetting object with some default settings.

    Args:
        dopants (array): dopant element names, e.g., ["Al", "N"]
        interstitial_coords (3x1 array): coordinations for interstitial sites,
                                         e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ..]
        is_antisite (bool): Whether to consider antisite defects.
        en_diff (float): Electronegativity (EN) difference for determining sets
                         of antisites and dopant sites.
        included (array): exceptionally added defect type with a charge state.
                          e.g., ["Va_O1_-1", "Va_O1_-2"]
        excluded (array): exceptionally removed defect type with a charge state.
                          e.g., ["Va_O1_1", "Va_O1_2"]
        distance (float): Maximum displacement distance in angstrom.
        cutoff (float): Cutoff radius for determining atoms displaced.
        symprec (float): Precision used for symmetry analysis.
    """

    def __init__(self, structure, dopants, interstitial_coords, is_antisite,
                 en_diff=_EN_DIFF, included="", excluded="", distance=_DISTANCE,
                 cutoff=_CUTOFF, symprec=_SYMPREC):

        self.dopants = dopants
        if interstitial_coords:
            if not len(interstitial_coords) % 3 == 0:
                raise ValueError("Interstitial coords is not a multiple of 3.")
        self.interstitial_coords = interstitial_coords
        self.is_antisite = is_antisite
        self.en_diff = en_diff
        self.included = included
        self.excluded = excluded
        self.distance = distance
        self.cutoff = cutoff
        self.symprec = symprec
        self.electronegativity = {}
        self.oxidation_states = {}

        self.antisite_configs = []
        self.dopant_configs = []

        # Set electronegativity and oxidation states for both intrinsic elements
        # and dopants
        elements_involved = structure.symbol_set + tuple(self.dopants)
        for e in elements_involved:
            try:
                self.electronegativity[e] = atom.electronegativity[e]
            except KeyError:
                warnings.warn("Electronegativity of " + e + " is unavailable.")
                self.electronegativity[e] = "N.A."
            try:
                self.oxidation_states[e] = atom.charge[e]
            except KeyError:
                warnings.warn("Oxidation state of " + e + " is unavailable.")
                self.oxidation_states[e] = "N.A."

        self.symmetrized_structure = \
            SpacegroupAnalyzer(structure).get_symmetrized_structure()
        # num_irreducible_sites (dict): number of irreducible sites for elements
        # num_irreducible_sites["Mg"] = 2: Mg has 2 inequivalent sites
        num_irreducible_sites = {}
        # irreducible_sites (array): a set of IrreducibleSite class objects
        self.irreducible_sites = []
        # equivalent_sites (list): List of list of pymatgen PeriodicSite class
        #                          object from SpacegroupAnalyzer.
        equiv_sites = self.symmetrized_structure.equivalent_sites
        # construct atomic indices of equivalent sites. E.g., 1..32

        for element in structure.symbol_set:
            num_irreducible_sites[element] = 0

        last_index = 0
        for i, e in enumerate(equiv_sites):
            element = e[0].species_string
            # increment the number of irreducible site for element
            num_irreducible_sites[element] += 1
            first_index = last_index + 1
            last_index = last_index + len(e)
            repr_coords = e[0].frac_coords
            irreducible_name = element + str(num_irreducible_sites[element])
            self.irreducible_sites.append(
                IrreducibleSite(irreducible_name, element, first_index,
                                last_index, repr_coords))

        en_keys = self.electronegativity.keys()

        # E.g., antisite_configs = [["Mg, "O"], ...]
        if is_antisite is True:
            for s1 in structure.symbol_set:
                for s2 in structure.symbol_set:
                    if s1 == s2:
                        continue
                    if s1 in en_keys and s2 in en_keys:
                        abs_diff = abs(self.electronegativity[s1] -
                                       self.electronegativity[s2])
                        if abs_diff < en_diff:
                            self.antisite_configs.append([s1, s2])
                    else:
                        cls.electronegativity_not_defined(s1, s2)

        # E.g., dopant_configs = [["Al", "Mg"], ...]
        if dopants:
            for d in dopants:
                if d in structure.symbol_set:
                    warnings.warn("Dopant " + d + " exists in host.")
                    continue
                for s1 in structure.symbol_set:
                    if d in en_keys and s1 in en_keys:
                        abs_diff = abs(self.electronegativity[d] -
                                       self.electronegativity[s1])
                        if abs_diff < en_diff:
                            self.dopant_configs.append([d, s1])
                    else:
                        cls.electronegativity_not_defined(d, s1)

        self.setting = \
            DefectSetting(self.symmetrized_structure, self.irreducible_sites,
                          self.dopant_configs, self.antisite_configs,
                          self.interstitial_coords, self.included,
                          self.excluded, self.distance, self.cutoff,
                          self.symprec, self.oxidation_states,
                          self.electronegativity)

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        """
        d = {"structure": self.structure,
             "irreducible_sites": [i.as_dict() for i in self.irreducible_sites],
             "dopant_configs": self.dopant_configs,
             "antisite_configs": self.antisite_configs,
             "interstitial_coords": self.interstitial_coords,
             "included": self.included,
             "excluded": self.excluded,
             "distance": self.distance,
             "cutoff": self.cutoff,
             "symprec": self.symprec,
             "oxidation_states": self.oxidation_states,
             "electronegativity": self.electronegativity}
        return d

    @classmethod
    def from_structure_file(cls, poscar, dopants=[], interstitial_coords=False,
                            is_antisite=False, en_diff=_EN_DIFF, included=None,
                            excluded=None, distance=_EN_DIFF, cutoff=_CUTOFF,
                            symprec=_SYMPREC):
        """
        Returns a json file.
        """
        structure = Structure.from_file(poscar)
        return cls(structure, dopants, interstitial_coords, is_antisite,
                   en_diff, included, excluded, distance, cutoff,
                   symprec)

    def to(self, defect_in_file="defect.in", poscar_file="DPOSCAR"):
        """
        Prints readable defect.in file.
        """
        self._write_defect_in(defect_in_file)
        # HACK: pmg has a bug, Symmetrized structure object cannot be converted
        #       to poscar
        Structure.from_str(self.symmetrized_structure.to(fmt="cif"),
                           fmt="cif").to(fmt="poscar", filename="DPOSCAR")

    def _write_defect_in(self, defect_in_file="defect.in"):
        with open(defect_in_file, 'w') as fw:

            for e in self.irreducible_sites:
                fw.write(
                    "  Irreducible element: {}\n".format(e.irreducible_name))
                fw.write("     Equivalent atoms: {}\n".format(
                    str(e.first_index) + ".." + str(e.last_index)))
                fw.write("Factional coordinates: %9.7f %9.7f %9.7f\n" \
                         % tuple(e.repr_coords))
                fw.write("    Electronegativity: {}\n".format(
                    self.electronegativity[e.element]))
                fw.write("      Oxidation state: {}\n\n".format(
                    self.oxidation_states[e.element]))
            fw.write("Interstitial coordinates: ")
            if self.interstitial_coords:
                fw.write(str([self.interstitial_coords[i:i + 3] for i in
                            range(0, len(self.interstitial_coords), 3)]) + "\n")
            else:
                fw.write("\n")

            if self.antisite_configs is not []:
                fw.write("Antisite defects: ")
                fw.write(' '.join(i[0] + "_" + i[1]
                                  for i in self.antisite_configs) + "\n\n")

            if self.dopants:
                for d in self.dopants:
                    if d not in self.symmetrized_structure.symbol_set:
                        fw.write("   Dopant element: {}\n".format(d))
                        fw.write("Electronegativity: {}\n".format(
                            self.electronegativity[d]))
                        fw.write("  Oxidation state: {}\n\n".format(
                            self.oxidation_states[d]))

                fw.write("Substituted defects: ")
                fw.write(' '.join(i[0] + "_" + i[1]
                                  for i in self.dopant_configs) + "\n\n")
            fw.write("Maximum Displacement: {}\n\n".format(self.distance))
            fw.write("Exceptionally included: {}\n".format(self.included))
            fw.write("Exceptionally excluded: {}\n".format(self.excluded))
            fw.write(
                "Cutoff region of atoms perturbed: {}\n".format(self.cutoff))
            fw.write("Symprec: {}\n".format(self.symprec))

    @staticmethod
    def print_dopant_info(dopant):
        """
        This is used for adding dopant information a posteriori.
        """
        try:
            electronegativity = atom.electronegativity[dopant]
        except KeyError:
            warnings.warn("Electronegativity of " + dopant + " is unavailable.")
            electronegativity = "N.A."

        try:
            oxidation_states = atom.charge[dopant]
        except KeyError:
            warnings.warn("Oxidation state of " + dopant + " is unavailable.")
            oxidation_states = "N.A."

        print("   Dopant element: {}".format(dopant))
        print("Electronegativity: {}".format(electronegativity))
        print("  Oxidation state: {}".format(oxidation_states))


class NotSupportedFlagError(Exception):
    pass


def get_electronegativity(s):
    try:
        return atom.electronegativity[s]
    except:
        warnings.warn("Electronegativity of " + s + " is unavailable.")
        return None


def get_oxidation_states(s):
    try:
        return atom.charge[s]
    except:
        warnings.warn("Oxidation state of " + s + " is unavailable.")
        return None

def main():
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name for a supercell used.")
    parser.add_argument("-d", "--dopants", dest="dopants", default="",
                        nargs="+", type=str,
                        help="Dopant elements. E.g., Ga In"
    parser.add_argument("-i", dest="interstitial_coords", nargs="+",
                        default=None, type=float,
                        help="Interstitial coordinates. E.g., 0.5 0.5 0.5.")
    parser.add_argument("-a", "--antisite", dest="is_antisite",
                        action="store_false",
                        help="Set when antisites are not considered.")
    parser.add_argument("-e", dest="en_diff", type=float, default=_EN_DIFF,
                        help="Criterion of the electronegativity difference \
                              determining antisites and/or impurities.")
    parser.add_argument("--included", dest="included", type=str, default="",
                        nargs="+",
                        help="Exceptionally included defects. E.g., Va_O2_-1.")
    parser.add_argument("--excluded", dest="excluded", type=str, default="",
                        nargs="+",
                        help="Exceptionally excluded defects. E.g., Va_O2_0.")
    parser.add_argument("--distance", dest="distance", type=float,
                        default=_DISTANCE,
                        help="Displacement distance. 0 means that \
                             random displacement is not considered.")
    parser.add_argument("--cutoff", dest="cutoff", type=float, default=_CUTOFF,
                        help="Set the cutoff radius [A] in which atoms are \
                              displaced.")
    parser.add_argument("--symprec", dest="symprec", type=float,
                        default=_SYMPREC,
                        help="Set precision used for symmetry analysis [A].")
    parser.add_argument("--print_dopant", dest="print_dopant", default=None,
                        type=str, help="Print dopant information that can be \
                                        added a posteriori.")
    opts = parser.parse_args()

    if opts.print_dopant:
        print_dopant_info(opts.print_dopant)
    else:
        defect_in = DefectInMaker.from_structure_file(
            opts.poscar, opts.dopants, opts.interstitial_coords,
            opts.is_antisite, opts.en_diff, opts.included,
            opts.excluded, opts.distance, opts.cutoff, opts.symprec)

        defect_in.to()

if __name__ == "__main__": main()
