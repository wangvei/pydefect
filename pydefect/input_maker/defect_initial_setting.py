# -*- coding: utf-8 -*-
import json
import warnings

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from monty.json import MontyEncoder
from monty.serialization import loadfn

import pydefect.core.atom as atom
from pydefect.core.irreducible_site import IrreducibleSite

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

# Following defaults determine the condition of automatic defect calculations.
# electronegativity difference for antisites and substitutional impurities
_EN_DIFF = 1.0
# Maximum displacement distance
_DISTANCE = 0.2
# Cutoff radius in which atoms are perturbed.
_CUTOFF = 3.0
_SYMPREC = 0.01


def extended_range(i):
    """
    Extended range method used for positive/negative input value.
    The input value is included for both positive and negative numbers.
    E.g., extended_range(3) = [0, 1, 2, 3]
          extended_range(-3) = [-3, -2, -1, 0]

    Args:
        i (int): an integer
    """
    if type(i) is not int:
        raise AttributeError
    if i >= 0:
        return range(i + 1)
    else:
        return range(i, 1)


def get_electronegativity(s):
    try:
        return atom.electronegativity[s]
    except KeyError:
        warnings.warn("Electronegativity of " + s + " is unavailable.")
        return None


def get_oxidation_state(s):
    try:
        return atom.charge[s]
    except KeyError:
        warnings.warn("Oxidation state of " + s + " is unavailable.")
        return None


def print_dopant_info(dopant):
    """
    This method is used to add dopant information a posteriori.
    """
    if Element.is_valid_symbol(dopant):
        electronegativity = get_electronegativity(dopant)
        oxidation_state = get_oxidation_state()

        print("   Dopant element: {}".format(dopant))
        print("Electronegativity: {}".format(electronegativity))
        print("  Oxidation state: {}".format(oxidation_state))
    else:
        warnings.warn(dopant + " is not a proper element name.")


class DefectInitialSetting:
    """
    This class object holds full information used for setting of a series of
    point defect calculations for a particular material.

    Args:
        structure (Structure): pmg Structure/IStructure class object for
                               perfect supercell
        irreducible_sites (list): list of IrreducibleSite class objects
        dopant_configs (Nx2 list): dopant configurations,
                                    e.g., [["Al", Mg"], ["N", "O"]
        antisite_configs (Nx2 list): antisite configurations,
                                  e.g., [["Mg","O"], ["O", "Mg"]]
        interstitial_coords (Nx3 list): coordinates for interstitial sites,
                                         e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ..]
        included (list): exceptionally added defects with charges
                          e.g., ["Va_O1_-1", "Va_O1_-2"]
        excluded (list): exceptionally removed defects with charges.
                          If they don't exist, this flag does nothing.
                          e.g., ["Va_O1_1", "Va_O1_2"]
        distance (float): Maximum displacement in angstrom.
                          0 means random displacement is switched off.
        cutoff (float): Cutoff radius in which atoms are displaced.
        symprec (float): Precision used for symmetry analysis.
        oxidation_states (dict): Oxidation states for relevant elements.
                                 Used to determine the default defect charges.
        electronegativity (dict): Electronegativity for relevant elements.
                                 Used to determine the substitutional defects.
    """

    def __init__(self, structure, irreducible_sites, dopant_configs,
                 antisite_configs, interstitial_coords, included, excluded,
                 distance, cutoff, symprec, oxidation_states,
                 electronegativity):

        self._structure = structure
        self._irreducible_sites = irreducible_sites
        self._dopant_configs = dopant_configs
        # dopant element names are derived from in_name of dopant_configs.
        self._dopants = set([d[0] for d in dopant_configs])
        self._antisite_configs = antisite_configs
        self._interstitial_coords = interstitial_coords
        self._included = included
        self._excluded = excluded
        self._distance = distance
        self._cutoff = cutoff
        self._symprec = symprec
        self._oxidation_states = oxidation_states
        self._electronegativity = electronegativity

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        # "return self.__dict__ == other.__dict__" is inapplicable,
        # because irreducible_sites returns pointers.
        print(self.as_dict())
        print(other.as_dict())
        return self.as_dict() == other.as_dict()

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectInitialSetting class object from a dictionary.
        """
        # Expansion of irreducible_sites is necessary.
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
        Constructs a DefectInitialSetting class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """
        Constructs DefectInitialSetting class object from a defect.in file.
        Currently, the file format of defect.in is not flexible, so be careful
        when modifying it by hand. The first word is mostly parsed. The number
        of words for each tag and the sequence of information is assumed to be
        fixed. See manual with care.
        E.g.,
            Irreducible element: Mg1
               Equivalent atoms: 1..32
         Fractional coordinates: 0.0000000 0.0000000 0.0000000
              Electronegativity: 1.31
                Oxidation state: 2
        """
        structure = Structure.from_file(poscar)
        irreducible_sites = []
        antisite_configs = []
        interstitial_coords = []
        included = []
        excluded = []
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
                    element = ''. \
                        join([i for i in irreducible_name if not i.isdigit()])
                    first_index, last_index = \
                        [int(i) for i in di.readline().split()[2].split("..")]
                    representative_coords = \
                        [float(i) for i in di.readline().split()[2:]]
                    irreducible_sites.append(
                        IrreducibleSite(irreducible_name,
                                        element,
                                        first_index,
                                        last_index,
                                        representative_coords))
                    electronegativity[element] = float(di.readline().split()[1])
                    oxidation_states[element] = int(di.readline().split()[2])

                elif line[0] == "Interstitial":
                    # "Interstitial coordinates: 0 0 0
                    # "Interstitial coordinates: 0.25 0.25 0.25"
                    if len(line) == 5:
                        interstitial_coords.append(
                            [float(line[i]) for i in range(2, len(line))])
                    else:
                        print("The number of interstitial coordinates is not a "
                              "multiple of 3.")
                        interstitial_coords = []

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

                elif line[0] == "Exceptionally" and line[1] == "included:":
                    included = line[2:]
                elif line[0] == "Exceptionally" and line[1] == "excluded:":
                    excluded = line[2:]

                else:
                    raise NotSupportedFlagError(l + " is not supported!")

        return cls(structure, irreducible_sites, dopant_configs,
                   antisite_configs, interstitial_coords, included, excluded,
                   distance, cutoff, symprec, oxidation_states,
                   electronegativity)

    @classmethod
    def from_basic_settings(cls, poscar, dopants=None,
                            flattened_interstitial_coords=None,
                            is_antisite=True, en_diff=_EN_DIFF, included="",
                            excluded="", distance=_DISTANCE, cutoff=_CUTOFF,
                            symprec=_SYMPREC):
        """
        Generates DefectInitialSetting object with default settings.

        Args:
            poscar (str): POSCAR type file name used for supercell defect
                          calculations
            dopants (list): dopant element names, e.g., ["Al", "N"]
            flattened_interstitial_coords (list):
                Coordinates for interstitial sites.
                The number of elements needs to be
                divided by 3.
                e.g., [0, 0, 0, 0.1, 0.1, 0.1]
            is_antisite (bool): Whether to consider antisite defects.
            en_diff (float): Electronegativity (EN) difference for determining
                             sets of antisites and dopant sites.
            included (list): Exceptionally added defects with charges
                            e.g., ["Va_O1_-1", "Va_O1_-2"]
            excluded (list): Exceptionally removed defects with charges.
                            If they don't exist, this flag does nothing.
                            e.g., ["Va_O1_1", "Va_O1_2"]
            distance (float): Maximum displacement distance in angstrom. 0 means
                              that random displacement is not considered.
            cutoff (float): Cutoff radius in which atoms are displaced.
            symprec (float): Precision used for symmetry analysis.
        """

        if dopants is None:
            dopants = []
        s = Structure.from_file(poscar).get_sorted_structure()
        symmetrized_structure = \
            SpacegroupAnalyzer(s, symprec=symprec).get_symmetrized_structure()

        # Electronegativity and oxidation states for constituents and dopants
        electronegativity = {}
        oxidation_states = {}

        symbol_set = s.symbol_set
        dopant_symbol_set = tuple(dopants)
        element_set = symbol_set + dopant_symbol_set

        for s in element_set:
            electronegativity[s] = get_electronegativity(s)
            oxidation_states[s] = get_oxidation_state(s)

        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irreducible_sites = {}

        # irreducible_sites (list): a set of IrreducibleSite class objects
        irreducible_sites = []

        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        equiv_sites = symmetrized_structure.equivalent_sites
        last_index = 0

        for i, equiv_site in enumerate(equiv_sites):
            # set element name of equivalent site
            element = equiv_site[0].species_string

            # increment number of inequivalent sites for element
            if element not in num_irreducible_sites.keys():
                num_irreducible_sites[element] = 1
            else:
                num_irreducible_sites[element] += 1

            first_index = last_index + 1
            last_index = last_index + len(equiv_site)
            # the following np.array type must be converted to list
            # to keep the consistency of the IrreducibleSite object.
            repr_coords = list(equiv_site[0].frac_coords)
            irreducible_name = element + str(num_irreducible_sites[element])

            irreducible_sites.append(IrreducibleSite(irreducible_name,
                                                     element,
                                                     first_index,
                                                     last_index,
                                                     repr_coords))

        # E.g., antisite_configs = [["Mg, "O"], ...]
        antisite_configs = []
        if is_antisite is True:
            for elem1 in symbol_set:
                for elem2 in symbol_set:
                    if elem1 == elem2:
                        continue
                    elif elem1 in electronegativity and \
                            elem2 in electronegativity:
                        if abs(electronegativity[elem1] -
                               electronegativity[elem2]) < en_diff:
                            antisite_configs.append([elem1, elem2])
                    else:
                        cls.electronegativity_not_defined(elem1, elem2)

        # E.g., dopant_configs = [["Al", "Mg"], ...]
        dopant_configs = []
        for dopant in dopants:
            if dopant in symbol_set:
                warnings.warn("Dopant " + dopant + " is a constituent of host.")
                continue
            for elem in symbol_set:
                if elem in electronegativity and dopant in electronegativity:
                    if abs(electronegativity[elem] -
                           electronegativity[dopant]) < en_diff:
                        dopant_configs.append([dopant, elem])
                else:
                    cls.electronegativity_not_defined(dopant, elem)

        # E.g., interstitial_coords = [[0, 0, 0], [0.1, 0.1, 0.1]]
        if flattened_interstitial_coords:
            if int(len(flattened_interstitial_coords)) % 3 == 0:
                interstitial_coords = \
                    [flattened_interstitial_coords[3*i: 3*i + 3]
                     for i in
                     range(int(len(flattened_interstitial_coords) / 3))]
            else:
                raise ValueError("The number of interstitial coordinates is not"
                                 " a multiple of 3.")
        else:
            interstitial_coords = []

        return cls(symmetrized_structure, irreducible_sites, dopant_configs,
                   antisite_configs, interstitial_coords, included, excluded,
                   distance, cutoff, symprec, oxidation_states,
                   electronegativity)

    def as_dict(self):
        """
        Dictionary representation of DefectInitialSetting class object.
        """
        d = {"structure":           self._structure,
             "irreducible_sites":   [i.as_dict()
                                     for i in self._irreducible_sites],
             "dopant_configs":      self._dopant_configs,
             "antisite_configs":    self._antisite_configs,
             "interstitial_coords": self._interstitial_coords,
             "included":            self._included,
             "excluded":            self._excluded,
             "distance":            self._distance,
             "cutoff":              self._cutoff,
             "symprec":             self._symprec,
             "oxidation_states":    self._oxidation_states,
             "electronegativity":   self._electronegativity}

        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def to(self, defect_in_file="defect.in", poscar_file="DPOSCAR"):
        """
        Prints readable defect.in file.
        """
        self._write_defect_in(defect_in_file)
        self._structure.to(fmt="poscar", filename=poscar_file)

    def make_defect_name_set(self):
        """
        Returns a set of defect names by default.
        """
        name_set = []

        # Vacancies
        # Charges are set from 0 to minus of the oxidation state.
        for i in self._irreducible_sites:
            os = self._oxidation_states[i.element]
            for o in extended_range(-os):
                name_set.append("Va_" + i.irreducible_name + "_" + str(o))

        # Interstitials
        # Both intrinsic elements and dopants are considered.
        # Charges are set from 0 to the oxidation state.
        inserted_elements = \
            tuple(self._structure.symbol_set) + tuple(self._dopants)
        for e in inserted_elements:
            os = self._oxidation_states[e]
            for j in range(len(self._interstitial_coords)):
                for o in extended_range(os):
                    name_set.append(e + "_i" + str(j + 1) + "_" + str(o))

        # Antisites + Substituted dopants
        for in_elem, out_elem in self._antisite_configs + self._dopant_configs:
            for i in self._irreducible_sites:
                if out_elem == i.element:
                    os_diff = self._oxidation_states[in_elem] - \
                              self._oxidation_states[out_elem]
                    for o in extended_range(os_diff):
                        name_set.append(
                            in_elem + "_" + i.irreducible_name + "_" + str(o))

        for i in self._included:
            if i not in name_set:
                name_set.append(i)
            else:
                print("{} is set to be included, but exists.".format(e))

        for e in self._excluded:
            if e in name_set:
                name_set.remove(e)
            else:
                print("{} is set to be excluded, but does not exist.".format(e))

        return name_set

    def _write_defect_in(self, defect_in_file="defect.in"):
        """
        Helper function to write down defect.in file.
        Args:
            defect_in_file (str): Name of defect.in type file.
        """
        with open(defect_in_file, 'w') as di:

            for e in self._irreducible_sites:
                di.write("   Irreducible element: {}\n".format(
                    e.irreducible_name))
                di.write("      Equivalent atoms: {}\n".format(
                    str(e.first_index) + ".." + str(e.last_index)))
                di.write("Fractional coordinates: %9.7f  %9.7f  %9.7f\n" %
                         tuple(e.repr_coords))
                di.write("     Electronegativity: {}\n".format(
                    self._electronegativity[e.element]))
                di.write("       Oxidation state: {}\n\n".format(
                    self._oxidation_states[e.element]))

            if self._interstitial_coords:
                for i in self._interstitial_coords:
                    di.write("Interstitial coordinates: %6.4f  %6.4f  %6.4f\n" %
                             tuple(i))
            else:
                di.write("Interstitial coordinates: \n")

            if self._antisite_configs is not []:
                di.write("Antisite defects: ")
                di.write(' '.join(i[0] + "_" + i[1]
                                  for i in self._antisite_configs) + "\n\n")

            for d in self._dopants:
                if d not in self._structure.symbol_set:
                    di.write("   Dopant element: {}\n".format(d))
                    di.write("Electronegativity: {}\n".format(
                        self._electronegativity[d]))
                    di.write("  Oxidation state: {}\n\n".format(
                        self._oxidation_states[d]))

            di.write("Substituted defects: ")
            di.write(' '.join(i[0] + "_" + i[1]
                              for i in self._dopant_configs) + "\n\n")

            di.write("Maximum Displacement: {}\n\n".format(self._distance))
            di.write("Exceptionally included: ")
            di.write(' '.join(i for i in self._included) + "\n")
            di.write("Exceptionally excluded: ")
            di.write(' '.join(i for i in self._excluded) + "\n")
            di.write(
                "Cutoff region of atoms perturbed: {}\n".format(self._cutoff))
            di.write("Symprec: {}\n".format(self._symprec))

    @staticmethod
    def electronegativity_not_defined(a, b):
        print("Electronegativity of {} and/or {} is not defined".format(a, b))

    @property
    def structure(self):
        return self._structure

    @property
    def irreducible_sites(self):
        return self._irreducible_sites

    @property
    def dopant_configs(self):
        return self._dopant_configs

    @property
    def dopants(self):
        return self._dopants

    @property
    def antisite_configs(self):
        return self._antisite_configs

    @property
    def interstitial_coords(self):
        return self._interstitial_coords

    @property
    def included(self):
        return self._included

    @property
    def excluded(self):
        return self._excluded

    @property
    def distance(self):
        return self._distance

    @property
    def cutoff(self):
        return self._cutoff

    @property
    def symprec(self):
        return self._symprec

    @property
    def oxidation_states(self):
        return self._oxidation_states

    @property
    def electronegativity(self):
        return self._electronegativity


class NotSupportedFlagError(Exception):
    pass


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
        defect_setting = DefectInitialSetting.from_basic_settings(
            opts.poscar, opts.dopants, opts.interstitial_coords,
            opts.is_antisite, opts.en_diff, opts.included, opts.excluded,
            opts.distance, opts.cutoff, opts.symprec)
        defect_setting.to()


if __name__ == "__main__":
    main()
