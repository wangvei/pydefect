#!/usr/bin/env python

import warnings
import json
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import atom
from defect import IrreducibleSite
from monty.json import MontyEncoder
from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

class DefectSetting():
    """
    This class object holds full information on the setting of the point 
    defect calculations.

    Args:
        structure: pmg Structure/IStructure class object
        irreducible_sites: IrreducibleSite class object
        dopant_configs (array): dopant configurations,
            e.g., ["Al_Mg", "N_O"].
        antisite_configs (array): antisite configuraions, 
            e.g., ["Mg_O", "O_Mg"].
        interstitial_coords (3x1 array): coordinations of interstitial sites,
            e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ...].
        included (array): exceptionally added defect type with a charge state.
            e.g., ["Va_O1_-1", "Va_O1_-2"]                     
        excluded (array): exceptionally removed defect type with a charge state.
                         In case some of them don't exist, they'll be ignored.
            e.g., ["Va_O1_1", "Va_O1_2"]                     
        distance (float): Maximum displacement distance in angstrom.
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
        # because irreducible_sites indicate just pointers.
        return self.as_dict() == other.as_dict()

    @classmethod
    def from_dict(cls, d):
        """
        Construct a DefectSetting class object from a dictionary.
        """
        irreducible_sites = []
        for i in d["irreducible_sites"]:
            irreducible_sites.append(IrreducibleSite.from_dict(i))

        #        return cls({k: v for k, v in d.items()})
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
        #        print(loadfn(filename))
        return cls.from_dict(loadfn(filename))

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """
        Construct DefectSetting class object from a defect.in file.
        Currently, format of the defect.in file is not flexible,
        so be careful for manipulating it.
        """
        structure = Structure.from_file(poscar)
        irreducible_sites = []
        electronegativity = {}
        oxidation_states = {}
        dopant_configs = []
        distance = None
        cutoff = None
        symprec = None

        with open(defect_in_file) as defect_in:
            for l in defect_in:
                line = l.split()

                if line == []:
                    continue

                elif line[0] == "Irreducible":
                    irreducible_name = line[2]
                    # remove index number from irreducible_name
                    # e.g., "Mg1" --> "Mg"
                    element = ''.join(
                               [i for i in irreducible_name if not i.isdigit()])
                    # Skip one line which shows the representative index.
                    # First_index atom is assumed to represent equivalent atoms.
                    defect_in.readline()
                    # Representative atom index
                    first_index, last_index = [int(i)
                           for i in defect_in.readline().split()[2].split("..")]
                    repr_coords = \
                           [float(i) for i in defect_in.readline().split()[2:]]
                    irreducible_sites.append(
                        IrreducibleSite(irreducible_name, element, first_index,
                                        last_index, repr_coords))
                    electronegativity[element] = \
                                          float(defect_in.readline().split()[1])
                    oxidation_states[element] = \
                                            int(defect_in.readline().split()[2])

                elif line[0] == "Interstital":
                    # "Interstital coordinates: 0 0 0" -> [0, 0, 0]
                    b = [float(''.join( i for i in line[i] if i.isdigit()
                                     or i == '.')) for i in range(2, len(line))]
                    # TODO: If the numbers for interstitial cannot be divided
                    #       by 3. Return Error.
                    interstitial_coords = [b[i:i + 3] for i in range(int(len(b) / 3))]

                elif line[0] == "Antisite":
                    antisite_configs = line[2:]
                    if antisite_configs is []:
                        antisite_configs = False

                elif line[0] == "Dopant":
                    d = line[2]
                    electronegativity[d] = \
                                          float(defect_in.readline().split()[1])
                    oxidation_states[d] = int(defect_in.readline().split()[2])

                elif line[0] == "Substituted":
                    dopant_configs = line[2:]
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
                    raise NotSupportedFlagError(
                        line[0] + " flag is not supported!")

        return cls(structure, irreducible_sites, dopant_configs,
                   antisite_configs, interstitial_coords, included, excluded,
                   distance, cutoff, symprec, oxidation_states,
                   electronegativity)

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


class DefectInMaker():
    """
    This class generates defect.in file based on some default setting.

    Args:
        dopants (array): dopant element names, e.g., ["Al", "N"]
        interstitial_coords (3x1 array): coordinations for interstitial sites,
            e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ...].
        is_antisite (bool): Whether to consider antisite defects.
        EN_diff (float): Electronegativity (EN) difference for determining sets 
                            of antisites and dopant sites. 
        included (array): exceptionally added defect type with a charge state.
            e.g., ["Va_O1_-1", "Va_O1_-2"]                     
        excluded (array): exceptionally removed defect type with a charge state.
            e.g., ["Va_O1_1", "Va_O1_2"]                     
        distance (float): Maximum displacement distance in angstrom.
        cutoff (float): Cutoff radius for detemining atoms displaced.
        symprec (float): Precision used for symmetry analysis.
    """

    def __init__(self, structure, dopants, interstitial_coords, is_antisite,
                 EN_diff, included="", excluded="", distance=0.2, cutoff=3.0,
                 symprec=0.01):

        self.dopants = dopants

        if interstitial_coords:
            if not len(interstitial_coords) % 3 == 0:
                comment = "The number of intersitial coordinates cannot be \
                         divided by 3."
                raise ValueError(comment)
        self.interstitial_coords = interstitial_coords
        self.is_antisite = is_antisite
        self.EN_diff = EN_diff
        self.included = included
        self.excluded = excluded
        self.distance = distance
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
        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irreducible_sites = {}
        # irreducible_sites (aray): a set of IrreducibleSite class objects
        self.irreducible_sites = []
        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        last = 0
        equiv_sites = self.symmetrized_structure.equivalent_sites
        for i, e in enumerate(equiv_sites):
            element = e[0].species_string
            if element not in num_irreducible_sites.keys():
                num_irreducible_sites[element] = 1
            else:
                # increment number of inequiv site for element
                num_irreducible_sites[element] += 1
            first = last + 1
            last = last + len(e)
            repr_coords = e[0].frac_coords
            irreducible_name = element + str(num_irreducible_sites[element])
            self.irreducible_sites.append(
                IrreducibleSite(irreducible_name, element, first, last, repr_coords))

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
                               self.electronegativity[s2]) < EN_diff:
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
                               self.electronegativity[d]) < EN_diff:
                            self.dopant_configs.append([d, s1])
                    else:
                        self.electronegativity_not_defined(d, s1)

        self.setting = DefectSetting(
            self.symmetrized_structure, self.irreducible_sites,
            self.dopant_configs, self.antisite_configs,
            self.interstitial_coords, self.included, self.excluded,
            self.distance, self.cutoff,
            self.symprec, self.oxidation_states,
            self.electronegativity)

    def electronegativity_not_defined(self, element1, element2):
        print("Electronegativity of {} and/or {} is not defined". \
              format(element1, element2))

    @classmethod
    def from_structure_file(cls, poscar, dopants=[], interstitial_coords=False,
                            is_antisite=False, EN_diff=1.0, included=None,
                            excluded=None, distance=0.2, cutoff=3.0, symprec=0.01):
        """
        Construct DefectInMaker class object from a POSCAR file.
        VERY IMPORTANT: Some parameters are set by default.
        """
        structure = Structure.from_file(poscar)

        return cls(structure, dopants, interstitial_coords, is_antisite,
                   EN_diff, included, excluded, distance, cutoff,
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

            for e in self.irreducible_sites:
                fw.write("  Irreducible element: {}\n".format(e.irreducible_name))
                fw.write("  Representative atom: {}\n".format(e.first_index))
                fw.write("     Equivalent atoms: {}\n".format(
                    str(e.first_index) + ".." + str(e.last_index)))
                fw.write("Factional coordinates: %9.7f %9.7f %9.7f\n" % tuple(e.repr_coords))
                fw.write("    Electronegativity: {}\n".format(
                    self.electronegativity[e.element]))
                fw.write("      Oxidation state: {}\n\n".format(
                    self.oxidation_states[e.element]))
            fw.write("Interstital coordinates: ")
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
                    if not d in self.symmetrized_structure.symbol_set:
                        fw.write("   Dopant element: {}\n".format(d))
                        fw.write("Electronegativity: {}\n".format(self.electronegativity[d]))
                        fw.write("  Oxidation state: {}\n\n".format(self.oxidation_states[d]))

                fw.write("Substituted defects: ")
                fw.write(' '.join(i[0] + "_" + i[1]
                                  for i in self.dopant_configs) + "\n\n")
            fw.write("Maximum Displacement: {}\n\n".format(self.distance))
            fw.write("Exceptionally included: {}\n".format(self.included))
            fw.write("Exceptionally excluded: {}\n".format(self.excluded))
            fw.write("Cutoff region of atoms perturbed: {}\n".format(self.cutoff))
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

        print("   Dopant element: {}".format(dopant))
        print("Electronegativity: {}".format(electronegativity))
        print("  Oxidation state: {}".format(oxidation_states))


class NotSupportedFlagError(Exception):
    pass


def main():
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str, help="POSCAR name.")
    parser.add_argument("-d", "--dopants", dest="dopants", default="", nargs="+",
                        type=str, help="Dopant elements. Eg. Al Ga In.")
    parser.add_argument("-i", dest="interstitial_coords", nargs="+", default=None,
                        type=float, help="Inetrstitials. Eg. 0.5 0.5 0.5.")
    parser.add_argument("-a", "--antisite", dest="is_antisite",
                        action="store_false",
                        help="Set if antisites are not considered.")
    parser.add_argument("-e", "--EN_diff", dest="EN_diff", type=float,
                        default=1.0, help="Criterion of the electronegativity \
                        difference for constructing antisite and impurities.")
    parser.add_argument("--included", dest="included", type=str, default="",
                        help="Exceptionally included defects. E.g. Va_O2_-1.")
    parser.add_argument("--excluded", dest="excluded", type=str, default="",
                        help="Exceptionally excluded defects. E.g. Va_O2_0.")
    parser.add_argument("--distance", dest="distance", type=float, default=0.2,
                        help="Displacement distance.")
    parser.add_argument("--cutoff", dest="cutoff", type=float, default=3.0,
                        help="Set the cutoff [A].")
    parser.add_argument("--symprec", dest="symprec", type=float, default=0.01,
                        help="Set the symprec [A].")
    parser.add_argument("--print_dopant", dest="print_dopant", default=None,
                        type=str, help="Print Dopant information.")

    opts = parser.parse_args()

    if opts.print_dopant:
        DefectInMaker.print_dopant_info(opts.print_dopant)
    else:
        defect_in = DefectInMaker.from_structure_file(opts.poscar, opts.dopants,
                                                      opts.interstitial_coords, opts.is_antisite,
                                                      opts.EN_diff, opts.included, opts.excluded,
                                                      opts.distance, opts.cutoff, opts.symprec)
        defect_in.to()


if __name__ == "__main__": main()
