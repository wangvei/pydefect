# -*- coding: utf-8 -*-
from collections import defaultdict
from itertools import permutations
import json
from typing import Union

from pydefect.core.defect_name import SimpleDefectName
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn, dumpfn

from obadb.util.structure_handler import get_point_group_from_dataset, \
    get_coordination_distances

from pydefect.core.error_classes import InvalidFileError
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.database.atom import electronegativity_list, oxidation_state_dict
from pydefect.util.logger import get_logger
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.core.config \
    import ELECTRONEGATIVITY_DIFFERENCE, DISPLACEMENT_DISTANCE, \
    CUTOFF_RADIUS, SYMMETRY_TOLERANCE, ANGLE_TOL

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def electronegativity_not_defined(a, b):
    logger.warning(
        "Electronegativity of {} and/or {} not defined".format(a, b))


def candidate_charge_set(i: int) -> list:
    """ Candidate charge set for defect charge states.

    The input value is included for both positive and negative numbers, and
    -1 (1) is included for positive (negative) odd number.
    E.g., potential_charge_set(3) = [-1, 0, 1, 2, 3]
          potential_charge_set(-3) = [-3, -2, -1, 0, 1]
          potential_charge_set(2) = [0, 1, 2]
          potential_charge_set(-4) = [-4, -3, -2, -1, 0]

    Args:
        i (int): an integer
    """
    if i % 2 == 0:
        return list(range(i + 1) if i >= 0 else range(i, 1))
    else:
        return list(range(-1, i + 1) if i >= 0 else range(i, 2))


def get_electronegativity(element) -> Union[float, None]:
    try:
        return electronegativity_list[element]
    except KeyError:
        logger.warning(
            f"Electronegativity of {element} is unavailable, so set to 0.")
        return


def get_oxidation_state(element) -> Union[int, None]:
    try:
        return oxidation_state_dict[element]
    except KeyError:
        logger.warning(
            f"Oxidation state of {element} is unavailable, so set to 0.")
        return


def dopant_info(dopant) -> Union[str, None]:
    """ Print dopant info.

    This method is used to add dopant info a posteriori.
    Args:
        dopant (str): Dopant element name e.g., Mg
    """
    if Element.is_valid_symbol(dopant):
        electronegativity = get_electronegativity(dopant)
        oxidation_state = get_oxidation_state(dopant)

        out = ["   Dopant element: {}".format(dopant),
               "Electronegativity: {}".format(electronegativity),
               "  Oxidation state: {}".format(oxidation_state)]
        return "\n".join(out)
    else:
        logger.warnings(dopant + " is not a proper element name.")
        return


def get_distances_from_string(string: list) -> dict:
    """
    Arg:
        string (list): Split string e.g. "Mg: 2.1 2.2 O: 2.3 2.4".split()

    Return:
        dict: e.g. {"Mg": [2,1, 2.2], "O": [2.3, 2.4]
    """
    distances = {}
    key = None
    for i in string:
        if i[-1] == ":":
            distances[i[:-1]] = []
            key = i[:-1]
        else:
            if key:
                distances[key].append(float(i))
            else:
                raise ValueError("Invalid string {} to create distance list".
                                 format(string))
    return distances


class DefectInitialSetting(MSONable):
    """ Holds full information for creating a series of DefectEntry objects."""

    def __init__(self,
                 structure: Structure,
                 space_group_symbol: str,
                 transformation_matrix: list,
                 cell_multiplicity: int,
                 irreducible_sites: list,
                 dopant_configs: list,
                 antisite_configs: list,
                 interstitial_site_names: list,
                 included: list,
                 excluded: list,
                 displacement_distance: float,
                 cutoff: float,
                 symprec: float,
                 angle_tolerance: float,
                 oxidation_states: dict,
                 electronegativity: dict,
                 interstitials_yaml: str = "interstitials.yaml"):
        """
        Args:
            structure (Structure):
                pmg Structure class object for perfect supercell
            space_group_symbol (str):
                space group symbol
            transformation_matrix (list):
                The matrix used for expanding the unitcell.
            cell_multiplicity (int):
                How much is the supercell larger than the *primitive* cell.
                This is used for calculating the defect concentration.
                Note that (multiplicity of the symmetry) * cell_multiplicity
                is the number of irreducible sites in the supercell.
            irreducible_sites (list):
                List of IrreducibleSite class objects
            dopant_configs (Nx2 list):
                Dopant configurations, e.g., [["Al", Mg"], ["N", "O"]]
            antisite_configs (Nx2 list):
                Antisite configurations, e.g., [["Mg","O"], ["O", "Mg"]]
            interstitial_site_names (list):
                Interstitial site indices written in interstitital.in file
            included (list):
                Exceptionally added defects with charges,
                e.g., ["Va_O1_-1", "Va_O1_-2"]
            excluded (list):
                Exceptionally removed defects with charges. If they don't exist,
                this flag does nothing. e.g., ["Va_O1_1", "Va_O1_2"]
            displacement_distance (float):
                Maximum displacement in angstrom. 0 means random displacement is
                switched off.
            cutoff (float):
                Cutoff radius in which atoms are displaced.
            symprec (float):
                Precision used for symmetry analysis.
            angle_tolerance (float):
                Angle tolerance for symmetry analysis in degree
            oxidation_states (dict):
                Oxidation states for relevant elements. Used to determine the
                default defect charges.
            electronegativity (dict):
                Electronegativity for relevant elements. Used to determine the
                substitutional defects.
            interstitials_yaml (str):
                Interstitial yaml file name
        """

        self.structure = structure
        self.space_group_symbol = space_group_symbol
        self.transformation_matrix = transformation_matrix
        self.cell_multiplicity = cell_multiplicity
        self.irreducible_sites = list(irreducible_sites)
        self.dopant_configs = list(dopant_configs)
        # dopant element names are derived from in_name of dopant_configs.
        self.dopants = set([d[0] for d in dopant_configs])
        self.antisite_configs = list(antisite_configs)
        self.interstitial_site_names = list(interstitial_site_names)
        self.included = list(included) if included else []
        self.excluded = list(excluded) if excluded else []
        self.displacement_distance = displacement_distance
        self.cutoff = cutoff
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.oxidation_states = oxidation_states
        self.electronegativity = electronegativity

        if self.interstitial_site_names:
            try:
                t = InterstitialSiteSet.from_files(self.structure,
                                                   interstitials_yaml)
            except FileNotFoundError:
                logger.error("interstitial.yaml file is needed.")
                raise

            sites = t.interstitial_sites
            if self.interstitial_site_names[0] == "all":
                self.interstitials = {k: v for k, v in sites.items()}
            else:
                for site_name in self.interstitial_site_names:
                    if site_name not in sites.keys():
                        raise ValueError(
                            f"Interstitial site name {site_name} "
                            f"do not exist in {interstitials_yaml}")
                self.interstitials = {k: v for k, v in sites.items()}
        else:
            self.interstitials = None

    @classmethod
    def from_dict(cls, d):
        # Expansion of irreducible_sites is necessary.
        irreducible_sites = []
        for i in d["irreducible_sites"]:
            irreducible_sites.append(IrreducibleSite.from_dict(i))

        structure = d["structure"]
        if isinstance(structure, dict):
            structure = Structure.from_dict(structure)

        return cls(structure=structure,
                   space_group_symbol=d["space_group_symbol"],
                   transformation_matrix=d["transformation_matrix"],
                   cell_multiplicity=d["cell_multiplicity"],
                   irreducible_sites=irreducible_sites,
                   dopant_configs=d["dopant_configs"],
                   antisite_configs=d["antisite_configs"],
                   interstitial_site_names=d["interstitial_site_names"],
                   included=d["included"],
                   excluded=d["excluded"],
                   displacement_distance=d["displacement_distance"],
                   cutoff=d["cutoff"],
                   symprec=d["symprec"],
                   angle_tolerance=d["angle_tolerance"],
                   oxidation_states=d["oxidation_states"],
                   electronegativity=d["electronegativity"])

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    @classmethod
    def from_defect_in(cls, poscar="DPOSCAR", defect_in_file="defect.in"):
        """ Constructor with defect.in file.

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
        interstitial_site_names = []
        included = None
        excluded = None
        electronegativity = {}
        oxidation_states = {}
        dopant_configs = []
        displacement_distance = None
        cutoff = None
        symprec = None

        with open(defect_in_file) as di:
            for l in di:
                line = l.split()

                if not line:
                    continue
                elif line[0] == "Space":
                    space_group_symbol = line[2]
                elif line[0] == "Transformation":
                    transformation_matrix = [int(i) for i in line[2:]]
                elif line[0] == "Cell":
                    cell_multiplicity = int(line[2])
                elif line[0] == "Irreducible":
                    irreducible_name = line[2]
                    # remove index from irreducible_name, e.g., "Mg1" --> "Mg"
                    element = ''. \
                        join([i for i in irreducible_name if not i.isdigit()])
                    wyckoff = di.readline().split()[2]
                    site_symmetry = di.readline().split()[2]

                    coordination_distances = \
                        get_distances_from_string(di.readline().split()[1:])
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
                    electronegativity[element] = \
                        float(di.readline().split()[1])
                    oxidation_states[element] = int(di.readline().split()[2])

                elif line[0] == "Interstitials:":
                    interstitial_site_names = list(line[1:])

                elif line[0] == "Antisite":
                    antisite_configs = [i.split("_") for i in line[2:]]

                elif line[0] == "Dopant":
                    d = line[2]
                    electronegativity[d] = float(di.readline().split()[1])
                    oxidation_states[d] = int(di.readline().split()[2])

                elif line[0] == "Substituted":
                    dopant_configs = [i.split("_") for i in line[2:]]

                elif line[0] == "Maximum":
                    displacement_distance = float(line[2])

                elif line[0] == "Cutoff":
                    cutoff = float(line[5])

                elif line[0] == "Symprec:":
                    symprec = float(line[1])

                elif line[0] == "Angle":
                    angle_tolerance = float(line[2])

                elif line[0] == "Exceptionally":
                    if line[1] == "included:":
                        included = line[2:]
                    elif line[1] == "excluded:":
                        excluded = line[2:]

                else:
                    raise InvalidFileError(
                        f"{line} is not supported in defect.in!")

        return cls(structure=structure,
                   space_group_symbol=space_group_symbol,
                   transformation_matrix=transformation_matrix,
                   cell_multiplicity=cell_multiplicity,
                   irreducible_sites=irreducible_sites,
                   dopant_configs=dopant_configs,
                   antisite_configs=antisite_configs,
                   interstitial_site_names=interstitial_site_names,
                   included=included,
                   excluded=excluded,
                   displacement_distance=displacement_distance,
                   cutoff=cutoff,
                   symprec=symprec,
                   angle_tolerance=angle_tolerance,
                   oxidation_states=oxidation_states,
                   electronegativity=electronegativity)

    @property
    def are_atoms_perturbed(self) -> bool:
        return False if self.cutoff < 1e-6 else True

    @classmethod
    def from_basic_settings(cls,
                            structure: Structure,
                            transformation_matrix: list,
                            cell_multiplicity: int,
                            oxidation_states: dict = None,
                            dopants: list = None,
                            is_antisite: bool = True,
                            en_diff: float = ELECTRONEGATIVITY_DIFFERENCE,
                            included: list = None,
                            excluded: list = None,
                            displacement_distance: float
                            = DISPLACEMENT_DISTANCE,
                            cutoff: float = CUTOFF_RADIUS,
                            symprec: float = SYMMETRY_TOLERANCE,
                            angle_tolerance: float = ANGLE_TOL,
                            interstitial_site_names: list = None):
        """ Generates object with some default settings.

        Args:
            structure (Structure/IStructure):
                Structure used for supercell defect calculations
            transformation_matrix (list):
                Diagonal component of transformation matrix.
            cell_multiplicity (int):
                How much is the supercell larger than the *primitive* cell.
            oxidation_states (dict):
                Oxidation states. Values are int.
            dopants (list):
                Dopant element names, e.g., ["Al", "N"]
            is_antisite (bool):
                Whether to consider antisite defects.
            en_diff (float):
                Electronegativity (EN) difference for determining sets of
                antisites and dopant sites.
            included (list):
                Exceptionally added defects with charges,
                e.g., ["Va_O1_-1", "Va_O1_-2"]
            excluded (list):
                this flag does nothing. e.g., ["Va_O1_1", "Va_O1_2"]
            displacement_distance (float):
                Maximum displacement distance in angstrom. 0 means that random
                displacement is not considered.
            cutoff (float):
                Cutoff radius in which atoms are displaced.
            symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
            interstitial_site_names (list):
                Names of interstitial sites written in interstitial.in file.
        """
        dopants = [] if dopants is None else dopants
        interstitial_site_names = [] \
            if interstitial_site_names is None else interstitial_site_names

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
        for element in element_set:
            electronegativity[element] = get_electronegativity(element)

        if oxidation_states is None:
            primitive = space_group_analyzer.find_primitive()
            oxi_guess = primitive.composition.oxi_state_guesses()[0]
            oxidation_states = {k: round(v) for k, v in oxi_guess.items()}
            for d in dopants:
                oxidation_states[d] = get_oxidation_state(d)
        else:
            for e in s.composition.elements:
                if e not in oxidation_states.keys():
                    raise ValueError("{} is not included in oxidation states".
                                     format(e))

        symmetrized_structure.add_oxidation_state_by_element(oxidation_states)

        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irreducible_sites = defaultdict(int)

        # irreducible_sites (list): a set of IrreducibleSite class objects
        irreducible_sites = []

        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        equiv_sites = symmetrized_structure.equivalent_sites

        lattice = symmetrized_structure.lattice.matrix

        sym_dataset = space_group_analyzer.get_symmetry_dataset()

        # Initialize the last index for iteration.
        last_index = 0
        for i, equiv_site in enumerate(equiv_sites):
            # set element name of equivalent site
            element = equiv_site[0].species_string
            # need to omit sign and numbers, eg Zn2+ -> Zn
            element = ''.join([i for i in element
                               if not (i.isdigit() or i == "+" or i == "-")])

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
                get_point_group_from_dataset(sym_dataset,
                                             representative_coords,
                                             lattice, symprec)[0]
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
            for elem1, elem2 in permutations(symbol_set, 2):
                if elem1 == elem2:
                    continue
                elif elem1 in electronegativity and elem2 in electronegativity:
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
                if elem and dopant in electronegativity:
                    if abs(electronegativity[elem] -
                           electronegativity[dopant]) < en_diff:
                        dopant_configs.append([dopant, elem])
                else:
                    electronegativity_not_defined(dopant, elem)

        return cls(structure=structure,
                   space_group_symbol=space_group_symbol,
                   transformation_matrix=transformation_matrix,
                   cell_multiplicity=cell_multiplicity,
                   irreducible_sites=irreducible_sites,
                   dopant_configs=dopant_configs,
                   antisite_configs=antisite_configs,
                   interstitial_site_names=interstitial_site_names,
                   included=included,
                   excluded=excluded,
                   displacement_distance=displacement_distance,
                   cutoff=cutoff,
                   symprec=symprec,
                   angle_tolerance=angle_tolerance,
                   oxidation_states=oxidation_states,
                   electronegativity=electronegativity)

    def as_dict(self):

        irreducible_sites = [i.as_dict() for i in self.irreducible_sites]

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "structure": self.structure,
             "space_group_symbol": self.space_group_symbol,
             "transformation_matrix": self.transformation_matrix,
             "cell_multiplicity": self.cell_multiplicity,
             "irreducible_sites": irreducible_sites,
             "dopant_configs": self.dopant_configs,
             "antisite_configs": self.antisite_configs,
             "interstitial_site_names": self.interstitial_site_names,
             "included": self.included,
             "excluded": self.excluded,
             "displacement_distance": self.displacement_distance,
             "cutoff": self.cutoff,
             "symprec": self.symprec,
             "angle_tolerance": self.angle_tolerance,
             "oxidation_states": self.oxidation_states,
             "electronegativity": self.electronegativity}

        return d

    def to_yaml_file(self, filename="defect.yaml"):
        dumpfn(self.as_dict(), filename)

    def to_json_file(self, filename):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def to(self, defect_in_file="defect.in", poscar_file="DPOSCAR"):
        """ Prints readable defect.in file and corresponding DPOSCAR file. """
        self._write_defect_in(defect_in_file)
        self.structure.to(fmt="poscar", filename=poscar_file)

    def make_defect_name_set(self) -> list:
        """ Return defect name list based on DefectInitialSetting object. """
        name_set = list()

        # Vacancies
        for i in self.irreducible_sites:
            # Vacancy charge is minus of element charge.
            charge = -self.oxidation_states[i.element]
            for c in candidate_charge_set(charge):
                in_atom = None
                out_site = i.irreducible_name
                name_set.append(SimpleDefectName(in_atom, out_site, c))

        inserted_elements = \
            tuple(self.structure.symbol_set) + tuple(self.dopants)
        for element in inserted_elements:
            if self.interstitials is not None:
                for name in self.interstitials.keys():
                    for charge in candidate_charge_set(
                            self.oxidation_states[element]):
                        in_atom = element
                        out_site = name
                        name_set.append(
                            SimpleDefectName(in_atom, out_site, charge))

        # Antisites + Substituted dopants
        for in_atom, out_elem in self.antisite_configs + self.dopant_configs:
            for i in self.irreducible_sites:
                if out_elem == i.element:
                    os_diff = self.oxidation_states[in_atom] - \
                              self.oxidation_states[out_elem]
                    for charge in candidate_charge_set(os_diff):
                        out_site = i.irreducible_name
                        name_set.append(
                            SimpleDefectName(in_atom, out_site, charge))

        for i in self.included:
            defect_name = SimpleDefectName.from_str(i)
            if defect_name not in name_set:
                name_set.append(defect_name)
            else:
                logger.warning(f"{defect_name} included, but already exists.")

        for e in self.excluded:
            defect_name = SimpleDefectName.from_str(e)
            if defect_name in name_set:
                name_set.remove(defect_name)
            else:
                logger.warning(f"{defect_name} excluded, but does not exist.")

        return name_set

    def _write_defect_in(self, defect_in_file="defect.in"):
        """ Helper method to write down defect.in file.

        Args:
            defect_in_file (str):
                Name of defect.in type file.
        """
        lines = list()
        lines.append(f"  Space group: {self.space_group_symbol}\n")

        lines.append("Transformation matrix: {0[0]:2d} {0[1]:2d} {0[2]:2d}".
                     format(self.transformation_matrix))

        lines.append(f"Cell multiplicity: {self.cell_multiplicity}\n")

        for site in self.irreducible_sites:
            lines.append(f"   Irreducible element: {site.irreducible_name}")
            lines.append(f"        Wyckoff letter: {site.wyckoff}")
            lines.append(f"         Site symmetry: {site.site_symmetry}")

            neighbor_atom_distances = []
            for k, v in site.coordination_distances.items():
                neighbor_atom_distances.append(
                   k + ": " + " ".join([str(round(i, 2)) for i in v]))
            neighbor_atom_distances = " ".join(neighbor_atom_distances)

            lines.append(f"          Coordination: {neighbor_atom_distances}")

            equiv_atoms = str(site.first_index) + ".." + str(site.last_index)
            lines.append(f"      Equivalent atoms: {equiv_atoms}")

            lines.append("Fractional coordinates: "
                         "{0[0]:9.7f}  {0[1]:9.7f}  {0[2]:9.7f}".
                         format(site.representative_coords))

            electronegativity = self.electronegativity[site.element]
            lines.append(f"     Electronegativity: {electronegativity}")

            oxidation_state = self.oxidation_states[site.element]
            lines.append(f"       Oxidation state: {oxidation_state}")

            # Append "" for one blank line.
            lines.append("")

        interstitial_site_names = " ".join(self.interstitial_site_names)
        lines.append(f"Interstitials: {interstitial_site_names}")

        # [["Ga", "N], ["N", "Ga"]] -> Ga_N N_Ga
        antisite_configs = \
            " ".join(["_".join(i) for i in self.antisite_configs])

        lines.append(f"Antisite defects: {antisite_configs}")
        lines.append("")

        for dopant in self.dopants:
            if dopant not in self.structure.symbol_set:
                electronegativity = self.electronegativity[dopant]
                oxidation_state = self.oxidation_states[dopant]

                lines.append(f"   Dopant element: {dopant}")
                lines.append(f"Electronegativity: {electronegativity}")
                lines.append(f"  Oxidation state: {oxidation_state}")
                lines.append("")

        substituted = " ".join("_".join(i) for i in self.dopant_configs)

        lines.append(f"Substituted defects: {substituted}")
        lines.append("")

        lines.append(f"Maximum Displacement: {self.displacement_distance}")
        lines.append("")

        lines.append(f"Exceptionally included: {' '.join(self.included)}")
        lines.append(f"Exceptionally excluded: {' '.join(self.excluded)}")
        lines.append("")

        lines.append(f"Cutoff region of neighboring atoms: {self.cutoff}")
        lines.append(f"Symprec: {self.symprec}")
        lines.append(f"Angle tolerance: {self.angle_tolerance}")
        lines.append("")

        with open(defect_in_file, 'w') as defect_in:
            defect_in.write("\n".join(lines))
