# -*- coding: utf-8 -*-
from collections import defaultdict
import json
import warnings

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from monty.json import MontyEncoder
from monty.serialization import loadfn

from obadb.util.structure_handler import get_symmetry_dataset, \
    get_point_group_from_dataset, get_coordination_distances

from pydefect.database.atom import electronegativity, charge
from pydefect.util.utils import get_logger
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.core.config \
    import ELECTRONEGATIVITY_DIFFERENCE, DISPLACEMENT_DISTANCE, CUTOFF_RADIUS, \
    SYMMETRY_TOLERANCE, ANGLE_TOL

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


logger = get_logger(__name__)


def electronegativity_not_defined(a, b):
    logger.warning("Electronegativity of {} and/or {} not defined".format(a, b))


def extended_range(i):
    """
    Extended range method used for positive/negative input value.
    The input value is included for both positive and negative numbers.
    E.g., extended_range(3) = [0, 1, 2, 3]
          extended_range(-3) = [-3, -2, -1, 0]

    Args:
        i (int): an integer
    """
    if isinstance(i, int):
        return range(i + 1) if i >= 0 else range(i, 1)
    else:
        raise AttributeError


def get_electronegativity(element):
    try:
        return electronegativity[element]
    except KeyError:
        logger.warning("Electronegativity of " + element + " is unavailable, " +
                       "so set to 0.")
        return 0


def get_oxidation_state(element):
    try:
        return charge[element]
    except KeyError:
        logger.warning("Oxidation state of " + element + " is unavailable," +
                       "so set to 0.")
        return 0


def print_dopant_info(dopant):
    """
    This method is used to add dopant information a posteriori.
    """
    if Element.is_valid_symbol(dopant):
        electronegativity = get_electronegativity(dopant)
        oxidation_state = get_oxidation_state(dopant)

        print("   Dopant element: {}".format(dopant))
        print("Electronegativity: {}".format(electronegativity))
        print("  Oxidation state: {}".format(oxidation_state))
    else:
        logger.warnings(dopant + " is not a proper element name.")


# def element_set(defect_initial_setting):
#     """
#     Args: defect_initial_setting (DefectInitialSetting)
#     """
#     e_set = set()

    # for i in defect_initial_setting.irreducible_sites:
    #     e_set.add(i.element)

    # for a in defect_initial_setting.antisite_configs:
    #     e_set.add(a[0])

    # for d in defect_initial_setting.dopant_configs:
    #     e_set.add(d[0])

    # return e_set


class DefectInitialSetting:
    """
    This class object holds full information used for setting of a series of
    point defect calculations for a particular material.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        space_group_symbol (str):
            space group symbol
        irreducible_sites (list):
            List of IrreducibleSite class objects
        dopant_configs (Nx2 list):
            Dopant configurations, e.g., [["Al", Mg"], ["N", "O"]]
        antisite_configs (Nx2 list):
            Antisite configurations, e.g., [["Mg","O"], ["O", "Mg"]]
        interstitial_coords (Nx3 list):
            Coordinates for interstitial sites, e.g., [[0, 0, 0], [0, 0, 0.1]]
        included (list):
            Exceptionally added defects with charges,
            e.g., ["Va_O1_-1", "Va_O1_-2"]
        excluded (list):
            Exceptionally removed defects with charges. If they don't exist,
            this flag does nothing. e.g., ["Va_O1_1", "Va_O1_2"]
        distance (float):
            Maximum displacement in angstrom. 0 means random displacement is
            switched off.
        cutoff (float):
            Cutoff radius in which atoms are displaced.
        symprec (float):
            Precision used for symmetry analysis.
        oxidation_states (dict):
            Oxidation states for relevant elements. Used to determine the
            default defect charges.
        electronegativity (dict):
            Electronegativity for relevant elements. Used to determine the
            substitutional defects.
    """

    def __init__(self, structure, space_group_symbol, irreducible_sites,
                 dopant_configs, antisite_configs, interstitial_coords,
                 included, excluded, distance, cutoff, symprec, angle_tolerance,
                 oxidation_states, electronegativity):

        self.structure = structure
        self.space_group_symbol= space_group_symbol
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
        self.angle_tolerance = angle_tolerance
        self.oxidation_states = oxidation_states
        self.electronegativity = electronegativity

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        # "return self.__dict__ == other.__dict__" is inapplicable,
        # because irreducible_sites returns pointers.
        # print(self.as_dict())
        # print(other.as_dict())
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

        return cls(d["structure"], d["space_group_symbol"], irreducible_sites,
                   d["dopant_configs"], d["antisite_configs"],
                   d["interstitial_coords"], d["included"], d["excluded"],
                   d["distance"], d["cutoff"], d["symprec"],
                   d["angle_tolerance"], d["oxidation_states"],
                   d["electronegativity"])

    @classmethod
    def load_json(cls, filename):
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
        space_group_symbol = None
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
                elif line[0] == "Space":
                    space_group_symbol = line[2]
                elif line[0] == "Irreducible":
                    irreducible_name = line[2]
                    # remove index from irreducible_name, e.g., "Mg1" --> "Mg"
                    element = ''. \
                        join([i for i in irreducible_name if not i.isdigit()])
                    wyckoff = di.readline().split()[2]
                    site_symmetry = di.readline().split()[2]

                    def get_distances(string):
                        distances = {}
                        for i in string:
                            if i[-1] == ":":
                                distances[i[:-1]] = []
                                key = i[:-1]
                            else:
                                distances[key].append(float(i))
                        return distances
                    coordination_distances = \
                        get_distances(di.readline().split()[1:])
                    first_index, last_index = \
                        [int(i) for i in di.readline().split()[2].split("..")]
                    representative_coords = \
                        [float(i) for i in di.readline().split()[2:]]
                    irreducible_sites.append(
                        IrreducibleSite(irreducible_name,
                                        element,
                                        first_index,
                                        last_index,
                                        representative_coords,
                                        wyckoff,
                                        site_symmetry,
                                        coordination_distances))
                    electronegativity[element] = float(di.readline().split()[1])
                    oxidation_states[element] = int(di.readline().split()[2])

                elif line[0] == "Interstitial":
                    # "Interstitial coordinates: 0 0 0
                    # "Interstitial coordinates: 0.25 0.25 0.25"
                    if len(line) == 5:
                        interstitial_coords.append(
                            [float(line[i]) for i in range(2, len(line))])
                    elif len(line) != 2:
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

                elif line[0] == "Angle":
                    angle_tolerance = float(line[2])

                elif line[0] == "Exceptionally" and line[1] == "included:":
                    included = line[2:]
                elif line[0] == "Exceptionally" and line[1] == "excluded:":
                    excluded = line[2:]

                else:
                    raise NotSupportedFlagError(l + " is not supported!")

        return cls(structure, space_group_symbol, irreducible_sites,
                   dopant_configs, antisite_configs, interstitial_coords,
                   included, excluded, distance, cutoff, symprec,
                   angle_tolerance, oxidation_states, electronegativity)

    @classmethod
    def from_basic_settings(cls,
                            structure,
                            dopants=None,
                            flattened_interstitial_coords=None,
                            is_antisite=True,
                            en_diff=ELECTRONEGATIVITY_DIFFERENCE,
                            included="",
                            excluded="",
                            distance=DISPLACEMENT_DISTANCE,
                            cutoff=CUTOFF_RADIUS,
                            symprec=SYMMETRY_TOLERANCE,
                            angle_tolerance=ANGLE_TOL):
        """
        Generates DefectInitialSetting object with some default settings.

        Args:
            structure (Structure/IStructure):
                Structure used for supercell defect calculations
            dopants (list):
                Dopant element names, e.g., ["Al", "N"]
            flattened_interstitial_coords (list):
                Coordinates for interstitial sites. The number of elements needs
                to be divided by 3. e.g., [0, 0, 0, 0.1, 0.1, 0.1]
            is_antisite (bool):
                Whether to consider antisite defects.
            en_diff (float):
                Electronegativity (EN) difference for determining sets of
                antisites and dopant sites.
            included (list):
                Exceptionally added defects with charges,
                e.g., ["Va_O1_-1", "Va_O1_-2"]
            excluded (list):
                Exceptionally removed defects with charges. If they don't exist,
                this flag does nothing. e.g., ["Va_O1_1", "Va_O1_2"]
            distance (float):
                Maximum displacement distance in angstrom. 0 means that random
                displacement is not considered.
            cutoff (float):
                Cutoff radius in which atoms are displaced.
            symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
        """

        if dopants is None:
            dopants = []
        s = structure.get_sorted_structure()
        space_group_analyzer = \
            SpacegroupAnalyzer(s, symprec=symprec,
                               angle_tolerance=angle_tolerance)
        symmetrized_structure = space_group_analyzer.get_symmetrized_structure()

        space_group_symbol = space_group_analyzer.get_space_group_symbol()

        symbol_set = s.symbol_set
        dopant_symbol_set = tuple(dopants)
        element_set = symbol_set + dopant_symbol_set

        # Electronegativity and oxidation states for constituents and dopants
        electronegativity = {}
        oxidation_states = {}
        for element in element_set:
            electronegativity[element] = get_electronegativity(element)
            oxidation_states[element] = get_oxidation_state(element)

        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irreducible_sites = defaultdict(int)

        # irreducible_sites (list): a set of IrreducibleSite class objects
        irreducible_sites = []

        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        equiv_sites = symmetrized_structure.equivalent_sites
        last_index = 0

        lattice = symmetrized_structure.lattice.matrix
        sym_dataset = get_symmetry_dataset(symmetrized_structure, symprec, angle_tolerance=angle_tolerance)

        for i, equiv_site in enumerate(equiv_sites):
            # set element name of equivalent site
            element = equiv_site[0].species_string

            # increment number of inequivalent sites for element
            num_irreducible_sites[element] += 1

            first_index = last_index + 1
            last_index = last_index + len(equiv_site)
            # the following np.array type must be converted to list
            # to keep the consistency of the IrreducibleSite object.
            representative_coords = list(equiv_site[0].frac_coords)
            irreducible_name = element + str(num_irreducible_sites[element])
            wyckoff = sym_dataset["wyckoffs"][first_index]

            simplified_point_group = \
                get_point_group_from_dataset(sym_dataset, representative_coords, lattice,
                                             symprec)[0]
            coordination_distances = \
                get_coordination_distances(symmetrized_structure, first_index)

            irreducible_sites.append(IrreducibleSite(irreducible_name,
                                                     element,
                                                     first_index,
                                                     last_index,
                                                     representative_coords,
                                                     wyckoff,
                                                     simplified_point_group,
                                                     coordination_distances))

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
                        electronegativity_not_defined(elem1, elem2)

        # E.g., dopant_configs = [["Al", "Mg"], ...]
        dopant_configs = []
        for dopant in dopants:
            if dopant in symbol_set:
                logger.warning("Dopant " + dopant + " constitutes host.")
                continue
            for elem in symbol_set:
                if elem in electronegativity and dopant in electronegativity:
                    if abs(electronegativity[elem] -
                           electronegativity[dopant]) < en_diff:
                        dopant_configs.append([dopant, elem])
                else:
                    electronegativity_not_defined(dopant, elem)

        # E.g., interstitial_coords = [[0, 0, 0], [0.1, 0.1, 0.1]]
        if flattened_interstitial_coords:
            if int(len(flattened_interstitial_coords)) % 3 == 0:
                interstitial_coords = \
                    [flattened_interstitial_coords[3*i: 3*i + 3]
                     for i in
                     range(int(len(flattened_interstitial_coords) / 3))]
            else:
                raise ValueError("Number of interstitial coordinates must be"
                                 " a multiple of 3.")
        else:
            interstitial_coords = []

        return cls(symmetrized_structure, space_group_symbol, irreducible_sites,
                   dopant_configs, antisite_configs, interstitial_coords,
                   included, excluded, distance, cutoff, symprec,
                   oxidation_states, electronegativity)

    def as_dict(self):
        """
        Dictionary representation of DefectInitialSetting class object.
        """
        d = {"structure":           self.structure,
             "space_group_symbol":  self.space_group_symbol,
             "irreducible_sites":   [i.as_dict()
                                     for i in self.irreducible_sites],
             "dopant_configs":      self.dopant_configs,
             "antisite_configs":    self.antisite_configs,
             "interstitial_coords": self.interstitial_coords,
             "included":            self.included,
             "excluded":            self.excluded,
             "distance":            self.distance,
             "cutoff":              self.cutoff,
             "symprec":             self.symprec,
             "angle_tolerance":     self.angle_tolerance,
             "oxidation_states":    self.oxidation_states,
             "electronegativity":   self.electronegativity}

        return d

    def to_yaml_file(self, filename="defect.yaml"):
        """
        Returns a yaml file.
        """
        from monty.serialization import dumpfn
        d = self.as_dict()
        print(d)
        dumpfn(d, filename)

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
        self.structure.to(fmt="poscar", filename=poscar_file)

    def make_defect_name_set(self):
        """
        Returns a set of defect names by default.
        """
        name_set = []

        # Vacancies
        # Charges are set from 0 to minus of the oxidation state.
        for i in self.irreducible_sites:
            os = self.oxidation_states[i.element]
            for o in extended_range(-os):
                name_set.append("Va_" + i.irreducible_name + "_" + str(o))

        # Interstitials
        # Both intrinsic elements and dopants are considered.
        # Charges are set from 0 to the oxidation state.
        inserted_elements = \
            tuple(self.structure.symbol_set) + tuple(self.dopants)
        for e in inserted_elements:
            os = self.oxidation_states[e]
            for j in range(len(self.interstitial_coords)):
                for o in extended_range(os):
                    name_set.append(e + "_i" + str(j + 1) + "_" + str(o))

        # Antisites + Substituted dopants
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
                logger.warning("{} is set to be included, but exists.".format(e))

        for e in self.excluded:
            if e in name_set:
                name_set.remove(e)
            else:
                logger.warning("{} is set to be excluded, but does not exist.".format(e))

        return name_set

    def _write_defect_in(self, defect_in_file="defect.in"):
        """
        Helper method to write down defect.in file.
        Args:
            defect_in_file (str):
                Name of defect.in type file.
        """
        with open(defect_in_file, 'w') as di:
            di.write("   Space group: {}\n\n".format(self.space_group_symbol))

            for e in self.irreducible_sites:
                di.write("   Irreducible element: {}\n".format(
                    e.irreducible_name))
                di.write("        Wyckoff letter: {}\n".format(
                    e.wyckoff))
                di.write("         Site symmetry: {}\n".format(
                    e.site_symmetry))
                c = ""
                for k, v in e.coordination_distances.items():
                    c += k + ": " + " ".join([str(round(i, 2)) for i in v]) + " "
                di.write("          Coordination: {}\n".format(c))
                di.write("      Equivalent atoms: {}\n".format(
                    str(e.first_index) + ".." + str(e.last_index)))
                di.write("Fractional coordinates: %9.7f  %9.7f  %9.7f\n" %
                         tuple(e.representative_coords))
                di.write("     Electronegativity: {}\n".format(
                    self.electronegativity[e.element]))
                di.write("       Oxidation state: {}\n\n".format(
                    self.oxidation_states[e.element]))

            if self.interstitial_coords:
                for i in self.interstitial_coords:
                    di.write("Interstitial coordinates: %6.4f  %6.4f  %6.4f\n" %
                             tuple(i))
            else:
                di.write("Interstitial coordinates: \n")

            if self.antisite_configs is not []:
                di.write("Antisite defects: ")
                di.write(' '.join(i[0] + "_" + i[1]
                                  for i in self.antisite_configs) + "\n\n")

            for d in self.dopants:
                if d not in self.structure.symbol_set:
                    di.write("   Dopant element: {}\n".format(d))
                    di.write("Electronegativity: {}\n".format(
                        self.electronegativity[d]))
                    di.write("  Oxidation state: {}\n\n".format(
                        self.oxidation_states[d]))

            di.write("Substituted defects: ")
            di.write(' '.join(i[0] + "_" + i[1]
                              for i in self.dopant_configs) + "\n\n")

            di.write("Maximum Displacement: {}\n\n".format(self.distance))
            di.write("Exceptionally included: ")
            di.write(' '.join(i for i in self.included) + "\n")
            di.write("Exceptionally excluded: ")
            di.write(' '.join(i for i in self.excluded) + "\n")
            di.write(
                "Cutoff region of atoms perturbed: {}\n".format(self.cutoff))
            di.write("Symprec: {}\n".format(self.symprec))
            di.write("Angle tolerance: {}\n".format(self.angle_tolerance))


class NotSupportedFlagError(Exception):
    pass


