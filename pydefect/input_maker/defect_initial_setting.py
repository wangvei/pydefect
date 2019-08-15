# -*- coding: utf-8 -*-
import json
from collections import defaultdict
from itertools import permutations
from typing import Union, List, Optional, Tuple

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn, dumpfn
from obadb.util.structure_handler import (
    get_point_group_from_dataset, get_coordination_distances)
from pydefect.core.config import (
    ELECTRONEGATIVITY_DIFFERENCE, DISPLACEMENT_DISTANCE, CUTOFF_RADIUS,
    SYMMETRY_TOLERANCE, ANGLE_TOL)
from pydefect.core.defect_entry import DefectType, DefectEntry
from pydefect.core.defect_name import DefectName
from pydefect.core.error_classes import InvalidFileError
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.database.atom import electronegativity_list, oxidation_state_dict
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import first_appearance_index
from pydefect.util.structure_tools import perturb_neighboring_atoms
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def candidate_charge_set(i: int) -> list:
    """ Candidate charge set for defect charge states.

    The input value is included for both positive and negative numbers, and
    -1 (1) is included for positive (negative) odd number.
    E.g., candidate_charge_set(3) = [-1, 0, 1, 2, 3]
          candidate_charge_set(-3) = [-3, -2, -1, 0, 1]
          candidate_charge_set(2) = [0, 1, 2]
          candidate_charge_set(-4) = [-4, -3, -2, -1, 0]

    Args:
        i (int): an integer
    """
    if i >= 0:
        charge_set = [i for i in range(i + 1)]
        if i % 2 == 1:
            charge_set.insert(0, -1)
    else:
        charge_set = [i for i in range(i, 1)]
        if i % 2 == 1:
            charge_set.append(1)

    return charge_set


def get_electronegativity(element: Union[str, Element]) -> float:
    """get a Pauling electronegativity obtained from wikipedia"""
    try:
        return electronegativity_list[str(element)]
    except KeyError:
        logger.warning(
            f"Electronegativity of {element} is unavailable, so set to 0.")
        return 0.0


def get_oxidation_state(element: Union[str, Element]) -> int:
    """get a typical oxidation stat """
    try:
        return oxidation_state_dict[str(element)]
    except KeyError:
        logger.warning(
            f"Oxidation state of {element} is unavailable, so set to 0.")
        return 0


def dopant_info(dopant: Union[str, Element]) -> Union[str, None]:
    """ Print dopant info.

    This method is used to add dopant info a posteriori.
    Args:
        dopant (str): Dopant element name e.g., Mg
    """
    dopant = str(dopant)
    if Element.is_valid_symbol(dopant):
        electronegativity = get_electronegativity(dopant)
        oxidation_state = get_oxidation_state(dopant)

        out = [f"   Dopant element: {dopant}",
               f"Electronegativity: {electronegativity}",
               f"  Oxidation state: {oxidation_state}"]
        return "\n".join(out)
    else:
        logger.warning(f"{dopant} is not an element name.")
        return


def get_distances_from_string(string: list) -> dict:
    """ Parse a string including distances to neighboring atoms
    Arg:
        string (list):
            List of string composed as "Mg: 2.1 2.2 O: 2.3 2.4".split()

    Return:
        dict: e.g. {"Mg": [2,1, 2.2], "O": [2.3, 2.4]
    """
    distances = {}
    key = None
    for i in string:
        if i[-1] == ":":
            distances[i[:-1]] = list()
            key = i[:-1]
        else:
            try:
                distances[key].append(float(i))
            except (UnboundLocalError, ValueError, NameError):
                print(f"Invalid string {string} for distance list.")
                raise
    return distances


def insert_atoms(structure: Structure, atoms: dict) -> Tuple[Structure, list]:
    """  Return structure with atoms and the indices atoms are inserted.
    """
    inserted_structure = structure.copy()
    inserted_atoms = []
    for k, v in atoms.items():
        index = first_appearance_index(inserted_structure, k)
        inserted_structure.insert(index, k, v)
        # The atom indices locating after k need to be incremented.
        for i in inserted_atoms:
            if i["index"] >= index:
                i["index"] += 1
        inserted_atoms.append({"element": k, "index": index, "coord": v})

    return inserted_structure, inserted_atoms


def select_defects(defect_set: List[dict],
                   keywords: Union[str, list, None] = None,
                   specified_defects: Optional[List[str]] = None) -> List[dict]:
    """ Returns names including one of keywords.

    Args:
        defect_set (list):
            A set of defect dict with "name" keys.
        keywords (str/list):
            Keywords determining if name is selected or not.
        specified_defects (list):
                Specifies particular defect to be considered.

    Return:
         list of defects.
    """

    if specified_defects:
        new_list = []
        for sd in specified_defects:
            name = DefectName.from_str(sd)
            for d in defect_set:
                print(name)
                if name.is_name_matched(d["name"]):
                    dd = d.copy()
                    dd["charge"] = name.charge
                    new_list.append(dd)
                    break
            else:
                raise ValueError(f"Specified defect {sd} is invalid.")
    elif keywords:
        new_list = []
        for defect in defect_set:
            name = DefectName(name=defect["name"], charge=defect["charge"])
            if name.is_name_matched(keywords):
                new_list.append(defect)
    else:
        return defect_set

    return new_list


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
                 interstitial_sites: Union[list, str],
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
                Structure for perfect supercell
            space_group_symbol (str):
                space group symbol, e.g., "Pm-3m".
            transformation_matrix (list):
                Matrix used for expanding the unitcell.
            cell_multiplicity (int):
                How much is the supercell larger than the *primitive* cell.
                This is used for calculating the defect concentration.
                (multiplicity of the symmetry) * cell_multiplicity corresponds
                to the number of irreducible sites in the supercell.
            irreducible_sites (list):
                List of IrreducibleSite objects
            dopant_configs (Nx2 list):
                Dopant configurations, e.g., [["Al", Mg"], ["N", "O"]]
            antisite_configs (Nx2 list):
                Antisite configurations, e.g., [["Mg","O"], ["O", "Mg"]]
            interstitial_sites (list):
                Interstitial site indices written in interstitial.in file
                "all" means that all the interstitials are considered.
            included (list):
                Exceptionally added defects with charges,
                e.g., ["Va_O1_-1", "Va_O1_-2"]
            excluded (list):
                Exceptionally removed defects with charges. If they don't exist,
                this flag does nothing. e.g., ["Va_O1_1", "Va_O1_2"]
            displacement_distance (float):
                Maximum displacement in angstrom applied to neighbor atoms.
                0 means that no random displacement is applied.
            cutoff (float):
                Cutoff radius in which atoms are displaced in angstrom.
            symprec (float):
                Precision used for symmetry analysis in angstrom.
            angle_tolerance (float):
                Angle tolerance for symmetry analysis in degree
            oxidation_states (dict):
                Oxidation states for relevant elements.
                Used to determine the default defect charges.
            electronegativity (dict):
                Electronegativity for relevant elements.
                Used to determine the substitutional defects.
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
        self.interstitial_sites = list(interstitial_sites)
        self.included = list(included) if included else list()
        self.excluded = list(excluded) if excluded else list()
        self.displacement_distance = displacement_distance
        self.cutoff = cutoff
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.oxidation_states = oxidation_states
        self.electronegativity = electronegativity
        self.interstitials_yaml = interstitials_yaml

        if not self.interstitial_sites:
            self.interstitials = None
        else:
            try:
                interstitial_site_set = \
                    InterstitialSiteSet.from_files(self.structure,
                                                   interstitials_yaml)
                sites = interstitial_site_set.interstitial_sites
            except FileNotFoundError:
                logger.error("interstitial.yaml file is needed.")
                raise

            if self.interstitial_sites[0] == "all":
                self.interstitials = dict(sites)
            else:
                self.interstitials = {}
                for i_site in self.interstitial_sites:
                    if i_site not in sites.keys():
                        raise ValueError(
                            f"Interstitial site name {i_site} "
                            f"do not exist in {interstitials_yaml}")
                    self.interstitials[i_site] = sites[i_site]

        self.defect_entries = None

    # see history at 2019/8/2 if from_dict and as_dict needs to recover

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    @classmethod
    def from_defect_in(cls,
                       poscar: str = "DPOSCAR",
                       defect_in_file: str = "defect.in"):
        """ Class object construction with defect.in file.

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

        Args:
            poscar (str):
                POSCAR type filename generated by from_basic_settings method.
                Usually, DPOSCAR
            defect_in_file (str):
                defect.in type filename.
        """
        structure = Structure.from_file(poscar)
        space_group_symbol = None
        irreducible_sites = list()
        antisite_configs = list()
        interstitial_sites = list()
        included = None
        excluded = None
        electronegativity = dict()
        oxidation_states = dict()
        dopant_configs = list()
        displacement_distance = None
        cutoff = None
        symprec = None

        with open(defect_in_file) as di:
            for l in di:
                line = l.split()

                # Vacant line is skipped.
                if not line:
                    continue
                elif line[0] == "Space":
                    space_group_symbol = line[-1]
                elif line[0] == "Transformation":
                    # Transformation matrix
                    transformation_matrix = [int(i) for i in line[-9:]]
                elif line[0] == "Cell":
                    cell_multiplicity = int(line[-1])
                elif line[0] == "Irreducible":
                    irreducible_name = line[-1]
                    # remove index from irreducible_name, e.g., "Mg1" --> "Mg"
                    element = "".join(
                        [i for i in irreducible_name if not i.isdigit()])

                    # Here, reading defect.in file is put forward.
                    # Wyckoff letter: line
                    split_next_line = di.readline().split()
                    wyckoff = split_next_line[-1]

                    # Site symmetry: line
                    split_next_line = di.readline().split()
                    site_symmetry = split_next_line[-1]

                    # Coordination: line
                    split_next_line = di.readline().split()
                    coordination_distances = \
                        get_distances_from_string(split_next_line[1:])

                    # Equivalent atoms: line
                    split_next_line = di.readline().split()
                    first_index, last_index = \
                        [int(i) for i in split_next_line[2].split("..")]

                    # Fractional coordinates: line
                    split_next_line = di.readline().split()
                    repr_coords = [float(i) for i in split_next_line[2:]]

                    irreducible_sites.append(
                        IrreducibleSite(irreducible_name,
                                        element,
                                        first_index,
                                        last_index,
                                        repr_coords,
                                        wyckoff,
                                        site_symmetry,
                                        coordination_distances))

                    # Electronegativity: line
                    split_next_line = di.readline().split()
                    electronegativity[element] = float(split_next_line[1])

                    # Oxidation state: line
                    split_next_line = di.readline().split()
                    oxidation_states[element] = int(split_next_line[2])

                elif line[0] == "Interstitials:":
                    interstitial_sites = list(line[1:])
                    if "all" in interstitial_sites:
                        interstitial_sites = "all"

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
                    cutoff = float(line[-1])

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
                   interstitial_sites=interstitial_sites,
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
                            interstitial_sites: list = None):
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
            interstitial_sites (list):
                Names of interstitial sites written in interstitial.in file.
        """
        dopants = dopants if dopants else list()
        interstitial_sites = \
            interstitial_sites if interstitial_sites else "all"

        s = structure.get_sorted_structure()

        sga = SpacegroupAnalyzer(structure=s, symprec=symprec,
                                 angle_tolerance=angle_tolerance)
        symmetrized_structure = sga.get_symmetrized_structure()
        space_group_symbol = sga.get_space_group_symbol()

        symbol_set = s.symbol_set
        dopant_symbol_set = tuple(dopants)
        element_set = symbol_set + dopant_symbol_set

        # Electronegativity and oxidation states for constituents and dopants
        electronegativity = {e: get_electronegativity(e) for e in element_set}

        if oxidation_states is None:
            primitive = sga.find_primitive()
            oxi_guess = primitive.composition.oxi_state_guesses()
            if oxi_guess:
                oxidation_states = \
                    {k: round(v) for k, v in oxi_guess[0].items()}
            else:
                oxidation_states = \
                    {str(e): get_oxidation_state(e)
                     for e in s.composition.elements}
                logger.warning(f"Oxidation states set to {oxidation_states}")
            for d in dopants:
                oxidation_states[d] = get_oxidation_state(d)
        else:
            for e in s.composition.elements:
                if e not in oxidation_states.keys():
                    raise ValueError(f"{e} not included in oxidation states.")

        symmetrized_structure.add_oxidation_state_by_element(oxidation_states)

        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites
        num_irreducible_sites = defaultdict(int)

        # irreducible_sites (list): a set of IrreducibleSite class objects
        irreducible_sites = list()

        # equivalent_sites: Equivalent site indices from SpacegroupAnalyzer.
        equiv_sites = symmetrized_structure.equivalent_sites
        lattice = symmetrized_structure.lattice.matrix
        sym_dataset = sga.get_symmetry_dataset()

        # Initialize the last index for iteration.
        last_index = 0
        for i, equiv_site in enumerate(equiv_sites):
            # set element name of equivalent site
            element = equiv_site[0].species_string
            # need to omit sign and numbers, eg Zn2+ -> Zn
            element = \
                ''.join([i for i in element if not (i.isdigit() or i in "+-")])

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

        def electronegativity_not_defined(a, b):
            logger.warning(f"Electronegativity of {a} and/or {b} not defined")

        # E.g., antisite_configs = [["Mg, "O"], ...]
        antisite_configs = list()
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
        dopant_configs = list()
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
                   interstitial_sites=interstitial_sites,
                   included=included,
                   excluded=excluded,
                   displacement_distance=displacement_distance,
                   cutoff=cutoff,
                   symprec=symprec,
                   angle_tolerance=angle_tolerance,
                   oxidation_states=oxidation_states,
                   electronegativity=electronegativity)

    def to_yaml_file(self, filename: str = "defect.yaml") -> None:
        dumpfn(self.as_dict(), filename)

    def to_json_file(self, filename: str) -> None:
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def to(self,
           defect_in_file: str = "defect.in",
           poscar_file: str = "DPOSCAR") -> None:
        """ Prints readable defect.in file and corresponding DPOSCAR file. """
        self._write_defect_in(defect_in_file)
        self.structure.to(fmt="poscar", filename=poscar_file)

    # def complex_set(self,
    #                 name,
    #                 removed_sites,
    #                 inserted_elements: dict,
    #                 extreme_charge) -> dict:
    #     return None

    def _substituted_set(self,
                         removed_sites: List[IrreducibleSite],
                         inserted_element: Optional[str]) -> list:
        defect_set = []
        changes_of_num_elements = defaultdict(int)
        for rs in removed_sites:
            in_name = inserted_element if inserted_element else "Va"
            name = "_".join([in_name, rs.irreducible_name])
            center = rs.representative_coords
            symmetry = rs.site_symmetry
            structure = self.structure.copy()
            structure.remove_sites([rs.first_index])
            removed_atoms = [{"element": rs.element,
                              "index": rs.first_index,
                              "coods": rs.representative_coords}]
            changes_of_num_elements[rs.element] += -1
            num_equiv_sites = rs.num_atoms

            if inserted_element:
                defect_type = DefectType.substituted
                atom = {inserted_element: rs.representative_coords}
                structure, inserted_atoms = insert_atoms(structure, atom)
                changes_of_num_elements[inserted_element] += 1
                inserted_oxi_state = self.oxidation_states[inserted_element]
            else:
                defect_type = DefectType.vacancy
                inserted_atoms = {}
                inserted_oxi_state = 0

            charges = candidate_charge_set(
                inserted_oxi_state - self.oxidation_states[rs.element])

            defect_set.append(
                {"name": name,
                 "defect_type": defect_type,
                 "initial_structure": structure,
                 "removed_atoms": removed_atoms,
                 "inserted_atoms": inserted_atoms,
                 "changes_of_num_elements": dict(changes_of_num_elements),
                 "initial_site_symmetry": symmetry,
                 "charges": charges,
                 "num_equiv_sites": num_equiv_sites,
                 "center": center})

        return defect_set

    def _inserted_set(self, inserted_elements: List[str]) -> list:
        if not self.interstitials:
            return []

        defect_set = []
        for name, i in self.interstitials.items():
            center = i.representative_coords
            symmetry = i.site_symmetry
            num_equiv_sites = i.symmetry_multiplicity

            for e in inserted_elements:
                changes_of_num_elements = defaultdict(int)
                name = "_".join([e, name])
                extreme_charge = self.oxidation_states[e]
                charges = candidate_charge_set(extreme_charge)
                changes_of_num_elements[e] += -1
                atom = {e: i.representative_coords}
                structure, inserted_atoms = \
                    insert_atoms(structure=self.structure, atoms=atom)

                defect_set.append(
                    {"name": name,
                     "defect_type": DefectType.interstitial,
                     "initial_structure": structure,
                     "removed_atoms": list(),
                     "inserted_atoms": inserted_atoms,
                     "changes_of_num_elements": dict(changes_of_num_elements),
                     "initial_site_symmetry": symmetry,
                     "charges": charges,
                     "num_equiv_sites": num_equiv_sites,
                     "center": center})

        return defect_set

    def make_defect_set(self,
                        keywords: Optional[list] = None,
                        specified_defects: Optional[list] = None):
        """ Return defect name list based on DefectInitialSetting object. """
        vacancies = self._substituted_set(removed_sites=self.irreducible_sites,
                                          inserted_element=None)

        substituted = []
        configs = self.antisite_configs + self.dopant_configs
        for in_elem, out_elem in configs:
            removed_sites = \
                [i for i in self.irreducible_sites if out_elem == i.element]
            substituted += self._substituted_set(removed_sites=removed_sites,
                                                 inserted_element=in_elem)
        inserted_elements = self.structure.symbol_set + tuple(self.dopants)
        interstitials = self._inserted_set(inserted_elements=inserted_elements)

        complexes = list()
        defect_set = vacancies + substituted + interstitials + complexes
        defect_set = select_defects(defect_set, keywords, specified_defects)

        self.defect_entries = []
        for defect in defect_set:
            center = defect.pop("center")
            for charge in defect.pop("charges"):
                inserted_indices = \
                    [i["index"] for i in defect["inserted_atoms"]]
                # By default, neighboring atoms are perturbed.
                # If one wants to avoid it, set displacement_distance = 0
                perturbed_structure, neighboring_sites = \
                    perturb_neighboring_atoms(
                        structure=defect["initial_structure"],
                        center=center,
                        cutoff=self.cutoff,
                        distance=self.displacement_distance,
                        inserted_atom_indices=inserted_indices)

                defect["initial_structure"].set_charge(charge)
                perturbed_structure.set_charge(charge)

                self.defect_entries.append(
                    DefectEntry(
                        perturbed_initial_structure=perturbed_structure,
                        cutoff=self.cutoff,
                        charge=charge,
                        neighboring_sites=neighboring_sites,
                        **defect))

    def _write_defect_in(self, defect_in_file: str = "defect.in"):
        """ Helper method to write down defect.in file.

        Args:
            defect_in_file (str):
                Name of defect.in type file.
        """
        lines = list()
        lines.append(f"  Space group: {self.space_group_symbol}\n")

        lines.append(f"Transformation matrix: "
                     f"{' '.join(str(x) for x in self.transformation_matrix)}")

        lines.append(f"Cell multiplicity: {self.cell_multiplicity}\n")

        for site in self.irreducible_sites:
            lines.append(f"   Irreducible element: {site.irreducible_name}")
            lines.append(f"        Wyckoff letter: {site.wyckoff}")
            lines.append(f"         Site symmetry: {site.site_symmetry}")

            neighbor_atom_distances = list()
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

        if self.interstitial_sites == "all":
            sites = "all"
        else:
            sites = ' '.join(self.interstitial_sites)
        lines.append(f"Interstitials: {sites}")

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
