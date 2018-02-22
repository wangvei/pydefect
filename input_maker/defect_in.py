#!/usr/bin/env python

import warnings
import json

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pydefect.input_maker.atom as atom
from pydefect.input_maker.defect import IrreducibleSite
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


def extended_range(i):
    """
    Extends range method used for positive/negative input value.
    The input value is included even for a positive number.
    E.g., extended_range(3) = [0, 1, 2, 3]
          extended_range(-3) = [-3, -2, -1, 0]

    Args:
        i (int): an integer
    """
    if not type(i) == int:
        raise AttributeError
    if i >= 0:
        return range(i + 1)
    else:
        return range(i, 1)


class DefectSetting:
    """
    This class object holds full information on the setting of the point 
    defect calculations.

    Args:
        structure (Structure): pmg Structure/IStructure class object
        irreducible_sites: IrreducibleSite class objects
        dopant_configs (Nx2 array): dopant configurations,
                                    e.g., [["Al", Mg"], ["N", "O"]
        antisite_configs (Nx2 array): antisite configurations,
                                  e.g., [["Mg","O"], ["O", "Mg"]]
        interstitial_coords (Nx3 array): coordinates of interstitial sites,
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
        # This must not be "return self.__dict__ == other.__dict__",
        # because irreducible_sites shows just pointers.
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
                    element = \
                        ''.join(
                            [i for i in irreducible_name if not i.isdigit()])
                    first_index, last_index = \
                        [int(i) for i in di.readline().split()[2].split("..")]
                    repr_coords = [float(i) for i in di.readline().split()[2:]]
                    irreducible_sites.append(IrreducibleSite(irreducible_name,
                                                             element,
                                                             first_index,
                                                             last_index,
                                                             repr_coords))
                    electronegativity[element] = float(di.readline().split()[1])
                    oxidation_states[element] = int(di.readline().split()[2])

                elif line[0] == "Interstitial":
                    # "Interstitial coordinates: 0 0 0 0.25 0.25 0.25"
                    #  --> [0, 0, 0, 0.25, 0.25, 0.25]
                    b = [float(''.join(i for i in line[i] if i.isdigit()
                                       or i == '.')) for i in
                         range(2, len(line))]
                    # If the number for interstitial coords is not divided by 3,
                    # return Error.
                    #  [0, 0, 0, 0.25, 0.25, 0.25]
                    #            --> [[0, 0, 0], [0.25, 0.25, 0.25]]
                    try:
                        interstitial_coords = [b[i:i + 3]
                                               for i in range(int(len(b) / 3))]
                    except ValueError:
                        print("Interstitial coordinates isn't a multiple of 3.")

                elif line[0] == "Antisite":
                    antisite_configs = [i.split("_") for i in line[2:]]

                elif line[0] == "Dopant":
                    d = line[2]
                    electronegativity[d] = float(di.readline().split()[1])
                    oxidation_states[d] = int(di.readline().split()[2])

                elif line[0] == "Substituted":
                    dopant_configs = [i.split("_") for i in line[2:]]

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
                            is_antisite=True, en_diff=_EN_DIFF, included="",
                            excluded="", distance=_DISTANCE, cutoff=_CUTOFF,
                            symprec=_SYMPREC):
        """
        Generates DefectSetting object with some default settings.

        Args:
            dopants (array): dopant element names, e.g., ["Al", "N"]
            is_antisite (bool): Whether to consider antisite defects.
            en_diff (float): Electronegativity (EN) difference for determining
                             sets of antisites and dopant sites.
        """

        structure = Structure.from_file(poscar)
        # Set electronegativity and oxidation states of elements and dopants
        electronegativity = {}
        oxidation_states = {}
        for s in structure.symbol_set + tuple(dopants):
            electronegativity[s] = get_electronegativity(s)
            oxidation_states[s] = get_oxidation_states(s)

        symmetrized_structure = \
            SpacegroupAnalyzer(structure,
                               symprec=symprec).get_symmetrized_structure()
        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irreducible_sites = {}
        # irreducible_sites (array): a set of IrreducibleSite class objects
        irreducible_sites = []
        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        last = 0
        equiv_sites = symmetrized_structure.equivalent_sites
        for i, e in enumerate(equiv_sites):
            element = e[0].species_string
            if element not in num_irreducible_sites.keys():
                num_irreducible_sites[element] = 1
            else:
                # increment number of inequivalent sites for element
                num_irreducible_sites[element] += 1
            first = last + 1
            last = last + len(e)
            repr_coords = e[0].frac_coords
            irreducible_name = element + str(num_irreducible_sites[element])
            irreducible_sites.append(IrreducibleSite(
                irreducible_name, element, first, last, repr_coords))

        en_keys = electronegativity.keys()

        # E.g., antisite_configs = [["Mg, "O"], ...]
        antisite_configs = []
        if is_antisite is True:
            for s1 in structure.symbol_set:
                for s2 in structure.symbol_set:
                    if s1 == s2:
                        continue
                    if s1 in en_keys and s2 in en_keys:
                        if abs(electronegativity[s1] - electronegativity[s2]) \
                                < en_diff:
                            antisite_configs.append([s1, s2])
                    else:
                        cls.electronegativity_not_defined(s1, s2)

        # E.g., dopant_configs = [["Al", "Mg"], ...]
        dopant_configs = []
        if dopants:
            for d in dopants:
                if d in structure.symbol_set:
                    warnings.warn("Dopant " + d + " exists in host.")
                    continue
                for s1 in structure.symbol_set:
                    if s1 in en_keys and d in en_keys:
                        if abs(electronegativity[s1] - electronegativity[d]) \
                                < en_diff:
                            dopant_configs.append([d, s1])
                    else:
                        cls.electronegativity_not_defined(d, s1)

        return cls(symmetrized_structure, irreducible_sites, dopant_configs,
                   antisite_configs, interstitial_coords, included, excluded,
                   distance, cutoff, symprec, oxidation_states,
                   electronegativity)

    def as_dict(self):
        """
        Dict representation of DefectSetting class object.
        """
        d = {"structure":           self.structure,
             "irreducible_sites":   [i.as_dict() for i in
                                     self.irreducible_sites],
             "dopant_configs":      self.dopant_configs,
             "antisite_configs":    self.antisite_configs,
             "interstitial_coords": self.interstitial_coords,
             "included":            self.included,
             "excluded":            self.excluded,
             "distance":            self.distance,
             "cutoff":              self.cutoff,
             "symprec":             self.symprec,
             "oxidation_states":    self.oxidation_states,
             "electronegativity":   self.electronegativity}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def to(self, defectin_file="defect.in", poscar_file="DPOSCAR"):
        """
        Prints readable defect.in file.
        """
        self._write_defect_in(defectin_file)
        # HACK: pmg has a bug, Symmetrized structure object cannot be poscar
        Structure.from_str(self.structure.to(fmt="cif"), fmt="cif").to(
            fmt="poscar", filename=poscar_file)

    def make_defect_name_set(self):
        """
        Returns a set of defect names by default.
        """
        name_set = []
        # Vacancies
        for i in self.irreducible_sites:
            os = self.oxidation_states[i.element]
            for o in extended_range(-os):
                name_set.append("Va_" + i.irreducible_name + "_" + str(o))

        # Interstitials
        inserted_elements = \
            tuple(self.structure.symbol_set) + tuple(self.dopants)
        for e in inserted_elements:
            os = self.oxidation_states[e]
            for j in range(len(self.interstitial_coords)):
                for o in extended_range(os):
                    name_set.append(e + "_i" + str(j + 1) + "_" + str(o))

        # antisites + substituted dopants
        for in_elem, out_elem in self.antisite_configs + self.dopant_configs:
            for i in self.irreducible_sites:
                if out_elem == i.element:
                    os_diff = self.oxidation_states[in_elem] - \
                              self.oxidation_states[out_elem]
                    for o in extended_range(os_diff):
                        name_set.append(
                            in_elem + "_" + i.irreducible_name + "_" + str(o))

        for i in self.included:
            if i not in name_set:
                name_set.append(i)
            else:
                print("{} is set to be included, but exists.".format(e))

        for e in self.excluded:
            if e in name_set:
                name_set.remove(e)
            else:
                print("{} is set to be excluded, but does not exist.".format(e))

        return name_set

    def _write_defect_in(self, defectin_file="defect.in"):
        with open(defectin_file, 'w') as fw:

            for e in self.irreducible_sites:
                fw.write("   Irreducible element: {}\n".format(
                    e.irreducible_name))
                fw.write("      Equivalent atoms: {}\n".format(
                    str(e.first_index) + ".." + str(e.last_index)))
                fw.write("Fractional coordinates: %9.7f %9.7f %9.7f\n" %
                         tuple(e.repr_coords))
                fw.write("     Electronegativity: {}\n".format(
                    self.electronegativity[e.element]))
                fw.write("       Oxidation state: {}\n\n".format(
                    self.oxidation_states[e.element]))

            fw.write("Interstitial coordinates: ")
            if self.interstitial_coords:
                fw.write(str(self.interstitial_coords) + "\n")
            else:
                fw.write("\n")

            if self.antisite_configs is not []:
                fw.write("Antisite defects: ")
                fw.write(' '.join(i[0] + "_" + i[1]
                                  for i in self.antisite_configs) + "\n\n")

            if self.dopants:
                for d in self.dopants:
                    if not d in self.structure.symbol_set:
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
            fw.write("Cutoff region of atoms perturbed: {}\n".format(
                self.cutoff))
            fw.write("Symprec: {}\n".format(self.symprec))

    @staticmethod
    def electronegativity_not_defined(e1, e2):
        print("Electronegativity of {} and/or {} is not defined".format(e1, e2))


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


def print_dopant_info(dopant):
    """
    This is used for adding dopant information a posteriori.
    """
    # TODO: check if the dopant is element.
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

    print("   Dopant element: {}".format(dopant))
    print("Electronegativity: {}".format(electronegativity))
    print("  Oxidation state: {}".format(oxidation_states))


def main():
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name.")
    parser.add_argument("-d", "--dopants", dest="dopants", default="",
                        nargs="+", type=str, help="Dopant elements. Eg. Ga In.")
    parser.add_argument("-i", dest="interstitial_coords", nargs="+",
                        default=None, type=float,
                        help="Interstitial coordinates. Eg., 0.5 0.5 0.5.")
    parser.add_argument("-a", "--antisite", dest="is_antisite",
                        action="store_false",
                        help="Set if antisites are not considered.")
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
        defect_setting = DefectSetting.from_basic_settings(
            opts.poscar, opts.dopants, opts.interstitial_coords,
            opts.is_antisite, opts.en_diff, opts.included,
            opts.excluded, opts.distance, opts.cutoff,
            opts.symprec)
        defect_setting.to()


if __name__ == "__main__": main()
