# -*- coding: utf-8 -*-
import json
from collections import defaultdict
from functools import reduce
from itertools import permutations,groupby
import operator
from typing import Union, List, Optional, Tuple, Dict

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn, dumpfn
from pydefect.core.complex_defects import ComplexDefects
from pydefect.core.config import (
    ELECTRONEGATIVITY_DIFFERENCE, DISPLACEMENT_DISTANCE, SYMMETRY_TOLERANCE,
    ANGLE_TOL, CUTOFF_FACTOR)
from pydefect.core.defect_entry import DefectType, DefectEntry
from pydefect.core.defect_name import DefectName
from pydefect.core.error_classes import InvalidFileError
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.core.irreducible_site import IrreducibleSite
from pydefect.database.atom import electronegativity_list, oxidation_state_dict
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import (
    get_min_distance, first_appearing_index, perturb_neighboring_atoms,
    defect_center_from_coords, get_coordination_distances, min_distance_from_coords)
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

logger = get_logger(__name__)


def default_charge_set(i: int) -> set:
    """Set of defect charge states.

    -1 (1) is included for positive (negative) odd number.
    E.g., default_charge_set(3) = [-1, 0, 1, 2, 3]
          default_charge_set(-3) = [-3, -2, -1, 0, 1]
          default_charge_set(2) = [0, 1, 2]
          default_charge_set(-4) = [-4, -3, -2, -1, 0]

    Args:
        i (int): an integer

    Return:
        Set of candidate charges
    """
    if i >= 0:
        charge_set = {i for i in range(i + 1)}
        if i % 2 == 1:
            charge_set.add(-1)
    else:
        charge_set = {i for i in range(i, 1)}
        if i % 2 == 1:
            charge_set.add(1)

    return set(charge_set)


def get_electronegativity(element: Union[str, Element]) -> float:
    """Get a Pauling electronegativity obtained from wikipedia

    If it doesn't exist in the database, 0.0 is returned.

    Args:
         element (str/ Element): Input element

    Return:
         Electronegativity
    """
    try:
        return electronegativity_list[str(element)]
    except KeyError:
        logger.warning(f"Electronegativity: {element} is unavailable. Set 0.")
        return 0.0


def get_oxidation_state(element: Union[str, Element]) -> int:
    """Get a typical oxidation state

    If it doesn't exist in the database, 0 is returned.

    Args:
         element (str/ Element): Input element

    Return:
         Oxidation state
     """
    try:
        return oxidation_state_dict[str(element)]
    except KeyError:
        logger.warning(f"Oxidation state: {element} is unavailable. Set 0.")
        return 0


def get_oxidation_states(dopants: list,
                         oxidation_states: dict,
                         structure: Structure) -> dict:
    """Get oxidation states.

    Args:
        dopants (list): Dopant element list ["Al", "N"]
        oxidation_states (dict): Dopant element list {"Al": 3, "N": -3, ..}
        structure (Structure): Host structure.

     The oxidation states are determined in the following order.
        1. Specified oxidation states in oxidation_states. If an element
           doesn't exist in oxidation_states, use get_oxidation_state.
        2. Pymatgen oxidation guess.
        3. get_oxidation_state function.

    Return
        Dict of oxidation_states.
    """
    host_element_set = structure.composition.elements
    element_set = host_element_set + dopants

    # Electronegativity and oxidation states for constituents and dopants
    if oxidation_states:
        d = {}
        for element in element_set:
            if element in oxidation_states.keys():
                d[element] = oxidation_states[element]
            else:
                d[element] = get_oxidation_state(element)
                logger.warning(f"Oxidation state of {element} is set to "
                               f"{d[element]}.")
    else:
        guess = structure.composition.reduced_composition.oxi_state_guesses()
        if guess:
            d = {k: round(v) for k, v in guess[0].items()}
            for dopant in dopants:
                d[dopant] = get_oxidation_state(dopant)
        else:
            d = {str(e): get_oxidation_state(e) for e in element_set}
            logger.warning(f"Oxidation states set to {d}")

    return d


def dopant_info(dopant: Union[str, Element]) -> Union[str, None]:
    """Print dopant info.

    This method is used to add dopant info a posteriori.

    Args:
        dopant (str): Dopant element name e.g., Mg

    Return:
        String of dopant information.
    """
    dopant = str(dopant)
    if Element.is_valid_symbol(dopant):
        out = [f"   Dopant element: {dopant}",
               f"Electronegativity: {get_electronegativity(dopant)}",
               f"  Oxidation state: {get_oxidation_state(dopant)}"]
        return "\n".join(out)
    else:
        logger.warning(f"{dopant} is not a proper element name.")
        return


def get_distances_from_string(string: list) -> dict:
    """Parse a string including distances to neighboring atoms
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


def insert_atoms(structure: Structure,
                 inserted_atoms: List[dict]
                 ) -> Tuple[Structure, List[dict]]:
    """Return the structure with inserted atoms and its inserted atom indices

    Args:
        structure (Structure):
            Input structure.
        inserted_atoms (list):
            List of dict with "element" and "coords" keys.
            Not that "index" is absent here as it is not determined, yet.

    Return:
        (inserted Structure, list of dict of inserted atom info)
    """
    inserted_structure = structure.copy()
    inserted_atoms_with_indices = []
    for atom in inserted_atoms:
        index = first_appearing_index(inserted_structure, atom["element"])
        inserted_structure.insert(index, atom["element"], atom["coords"])
        # The atom indices locating after k need to be incremented.
        for i in inserted_atoms_with_indices:
            if i["index"] >= index:
                i["index"] += 1

        inserted_atoms_with_indices.append({"element": atom["element"],
                                            "index": index,
                                            "coords": atom["coords"]})

    return inserted_structure, inserted_atoms_with_indices


def select_defects(defect_set: Dict[str, dict],
                   keywords: Union[str, list, None] = None,
                   included: Optional[list] = None,
                   excluded: Optional[list] = None,
                   specified_defects: Optional[List[str]] = None
                   ) -> Dict[str, dict]:
    """Returns defect names that contains at least one of keywords.

    Args:
        defect_set (list):
            A list of defect dict with "name" and "charge" keys.
            e.g.  {"Va_Mg1": {"charges": {-1, 0, 1, 2}},
                   "Va_O1": {"charges": {0}},
                   "O_i1": {"charges": {0, 1}}}
        keywords (str/list):
            Keywords determining if name is selected or not.
        included (list):
            Exceptionally included defects with full names.
        excluded (list):
            Exceptionally excluded defects with full names.
        specified_defects (list):
            Specifies particular defect to be considered.
            In this case, only name needs to exist in the candidate.
            For example, Va_O1_-1 is fine while Va_O2_-1 is not in the above
            defect_set.

    Return:
         list of defect dict. Charges are set
           e.g.  {"Va_Mg1": {"charges": {0, 1, 2}},
                  "Va_O1": {"charges": {0}},..}
    """
    included = included or []
    excluded = excluded or []
    specified_defects = specified_defects or []

    for name, defect in defect_set.items():
        if specified_defects:
            new_charges = set()
            for i in specified_defects:
                defect_name = DefectName.from_str(i)
                if defect_name.is_name_matched(name):
                    new_charges.add(defect_name.charge)
            defect["charges"] = new_charges
            continue

        for e in excluded:
            defect_name = DefectName.from_str(e)
            if defect_name.is_name_matched(name):
                defect["charges"].discard(defect_name.charge)

        for i in included:
            defect_name = DefectName.from_str(i)
            if defect_name.is_name_matched(name):
                defect["charges"].add(defect_name.charge)

        if keywords:
            new_charges = set()
            for charge in defect["charges"]:
                defect_name = DefectName(name=name, charge=charge)
                if defect_name.is_name_matched(keywords):
                    new_charges.add(charge)
            defect["charges"] = new_charges

    return defect_set


class DefectInitialSetting(MSONable):
    """Holds full information for creating a series of DefectEntry objects."""

    def __init__(self,
                 structure: Structure,
                 space_group_symbol: str,
                 transformation_matrix: list,
                 cell_multiplicity: int,
                 irreducible_sites: list,
                 dopant_configs: list,
                 antisite_configs: list,
                 interstitial_site_names: Union[list, str],
                 complex_defect_names: Union[list, str],
                 included: Optional[list],
                 excluded: Optional[list],
                 displacement_distance: float,
                 symprec: float,
                 angle_tolerance: float,
                 oxidation_states: dict,
                 electronegativity: dict,
                 interstitials_yaml: str = "interstitials.yaml",
                 complex_defect_yaml: str = "complex_defects.yaml",
                 defect_entries: List[DefectEntry] = None):
        """
        Args:
            structure (Structure):
                Structure for perfect supercell.
            space_group_symbol (str):
                Space group symbol, e.g., "Pm-3m".
            transformation_matrix (list):
                Matrix used for expanding the *primitive* unitcell.
            cell_multiplicity (int):
                How much is the supercell larger than the *primitive* unitcell.
                This is used for calculating the defect concentration.
                symmetry_multiplicity * cell_multiplicity equals to the number
                of irreducible sites in the supercell.
            irreducible_sites (list):
                List of IrreducibleSite objects
            dopant_configs (Nx2 list):
                Dopant configurations, e.g., [["Al", Mg"], ["N", "O"]] meaning
                Al-on-Mg and N-on-O substitutions.
            antisite_configs (Nx2 list):
                Antisite configurations, e.g., [["Mg","O"], ["O", "Mg"]]
            interstitial_site_names (list/str):
                Interstitial site names written in the interstitial.yaml file
                "all" means all the interstitials are considered.
            complex_defect_names (list/str):
                Complex defect names in complex_defects.yaml file.
            included (list):
                Exceptionally included defects with each charge state,
                e.g., ["Va_O1_-1", "Va_O1_-2"]
            excluded (list):
                Exceptionally removed defects with each charge state.
                If they don't exist in the default set, this flag does nothing.
                e.g., ["Va_O1_1", "Va_O1_2"]
            displacement_distance (float):
                Maximum displacement in angstrom applied to neighbor atoms.
                Note that defect itself, e.g., the Al atom of Al-on-Mg is not
                displaced.
                0 means that no random displacement is applied.
            symprec (float):
                Precision used for symmetry analysis in angstrom.
            angle_tolerance (float):
                Angle tolerance for symmetry analysis in degree.
            oxidation_states (dict):
                Oxidation states for relevant elements. E.g., {"Mg": 2}
                Used to determine the default defect charges.
                If not set, use pymatgen oxidation guess.
            electronegativity (dict):
                Electronegativity for relevant elements.
                Used to determine the substitutional defects.
            interstitials_yaml (str):
                Interstitial yaml file name
            complex_defect_yaml (str):
                Complex defect yaml file name
        """

        self.structure = structure
        self.space_group_symbol = space_group_symbol
        if len(transformation_matrix) == 9:
            self.transformation_matrix = transformation_matrix
        elif len(transformation_matrix) == 3:
            self.transformation_matrix = sum(transformation_matrix, [])
        else:
            raise ValueError(f"transformation matrix {transformation_matrix} "
                             f"is invalid.")

        self.cell_multiplicity = cell_multiplicity
        self.irreducible_sites = irreducible_sites[:]
        self.dopant_configs = dopant_configs[:]
        # dopant element names are derived from in_name of dopant_configs.
        self.dopants = set([d[0] for d in dopant_configs])
        self.antisite_configs = antisite_configs[:]
        if isinstance(interstitial_site_names, list):
            self.interstitial_site_names = interstitial_site_names[:]
        else:
            self.interstitial_site_names = [interstitial_site_names]
        if isinstance(complex_defect_names, list):
            self.complex_defect_names = complex_defect_names[:]
        else:
            self.complex_defect_names = [complex_defect_names]
        self.included = included[:] if included else []
        self.excluded = excluded[:] if excluded else []
        self.displacement_distance = displacement_distance
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance
        self.oxidation_states = oxidation_states
        self.electronegativity = electronegativity
        self.interstitials_yaml = interstitials_yaml
        self.complex_defect_yaml = complex_defect_yaml

        if not self.interstitial_site_names:
            self.interstitials = {}
        else:
            try:
                interstitial_site_set = InterstitialSiteSet.from_files(
                    self.structure, interstitials_yaml).interstitial_sites
            except FileNotFoundError:
                if not self.interstitial_site_names[0] == "all":
                    logger.error("interstitial.yaml file is needed.")
                    raise
                interstitial_site_set = []

            if self.interstitial_site_names[0] == "all":
                self.interstitials = dict(interstitial_site_set)
            else:
                self.interstitials = {}
                for i_site in self.interstitial_site_names:
                    if i_site not in interstitial_site_set.keys():
                        raise ValueError(
                            f"Interstitial site name {i_site} "
                            f"does not exist in {interstitials_yaml}")
                    self.interstitials[i_site] = interstitial_site_set[i_site]

        self.complex_defects = {}
        try:
            complexes = ComplexDefects.from_files(
                self.structure, complex_defect_yaml).complex_defects
        except FileNotFoundError:
            complexes = {}

        for c in self.complex_defect_names:
            if c not in complexes.keys():
                raise ValueError(
                    f"Complex defect name {c} does not exist in "
                    f"{complex_defect_yaml}.")
            self.complex_defects[c] = complexes[c]

        self.defect_entries = defect_entries or None

    # see history at 2019/8/2 if from_dict and as_dict need to be recovered.
    @classmethod
    def load_json(cls, filename: str) -> "DefectInitialSetting":
        return loadfn(filename)

    @classmethod
    def from_defect_in(cls,
                       poscar: str = "DPOSCAR",
                       defect_in_file: str = "defect.in",
                       interstitials_yaml: str = "interstitials.yaml",
                       complex_defect_yaml: str = "complex_defects.yaml",
                       ) -> "DefectInitialSetting":
        """Class object construction with defect.in file.

        Currently, the file format of defect.in is not flexible, so be careful
        when modifying it by hand. The first word is mostly parsed. The number
        of words for each tag and the sequence of information is assumed to be
        fixed. See manual with care. Example of defect.in is as follows.
        ------------------------------------------------------------------------
        Space group: Fm-3m

        Transformation matrix: 2 0 0 0 2 0 0 0 2
        Cell multiplicity: 32

           Irreducible element: Mg1
                Wyckoff letter: a
                 Site symmetry: m-3m
                 Cutoff radius: 2.5
                  Coordination: O: 2.1 2.1 2.1 2.1 2.1 2.1
              Equivalent atoms: 0..31
        Fractional coordinates: 0.0000000  0.0000000  0.0000000
             Electronegativity: 1.31
               Oxidation state: 2

           Irreducible element: O1
                Wyckoff letter: b
                 Site symmetry: m-3m
                 Cutoff radius: 2.5
                  Coordination: Mg: 2.1 2.1 2.1 2.1 2.1 2.1
              Equivalent atoms: 32..63
        Fractional coordinates: 0.2500000  0.0000000  0.0000000
             Electronegativity: 3.44
               Oxidation state: -2

        Interstitials: all
        Complex defects:
        Antisite defects:

        Substituted defects:

        Maximum Displacement: 0.0

        Exceptionally included:
        Exceptionally excluded:

        Symprec: 0.01
        Angle tolerance: 5
        ------------------------------------------------------------------------

        Args:
            poscar (str):
                POSCAR type filename generated by from_basic_settings method.
                Usually, called as DPOSCAR.
            defect_in_file (str):
                defect.in type filename.
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
        symprec = None

        with open(defect_in_file) as defect_in:
            for l in defect_in:
                line = l.split()

                if not line:
                    continue  # Vacant line is skipped.
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
                    split_next_line = defect_in.readline().split()
                    wyckoff = split_next_line[-1]

                    split_next_line = defect_in.readline().split()
                    site_symmetry = split_next_line[-1]

                    split_next_line = defect_in.readline().split()
                    cutoff_radius = float(split_next_line[-1])

                    split_next_line = defect_in.readline().split()
                    coordination_distances = \
                        get_distances_from_string(split_next_line[1:])

                    split_next_line = defect_in.readline().split()
                    first_index, last_index = \
                        [int(i) for i in split_next_line[2].split("..")]

                    split_next_line = defect_in.readline().split()
                    repr_coords = [float(i) for i in split_next_line[2:]]

                    irreducible_sites.append(
                        IrreducibleSite(irreducible_name,
                                        element,
                                        first_index,
                                        last_index,
                                        repr_coords,
                                        wyckoff,
                                        site_symmetry,
                                        cutoff_radius,
                                        coordination_distances))

                    split_next_line = defect_in.readline().split()
                    electronegativity[element] = float(split_next_line[1])

                    split_next_line = defect_in.readline().split()
                    oxidation_states[element] = int(split_next_line[2])

                elif line[0] == "Interstitials:":
                    interstitial_site_names = line[1:]
                    if "all" in interstitial_site_names:
                        interstitial_site_names = "all"

                elif line[0] == "Complex":
                    complex_defect_names = line[2:]

                elif line[0] == "Antisite":
                    antisite_configs = [i.split("_") for i in line[2:]]

                elif line[0] == "Dopant":
                    d = line[2]
                    electronegativity[d] = \
                        float(defect_in.readline().split()[1])
                    oxidation_states[d] = int(defect_in.readline().split()[2])

                elif line[0] == "Substituted":
                    dopant_configs = [i.split("_") for i in line[2:]]

                elif line[0] == "Maximum":
                    displacement_distance = float(line[2])

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
                   complex_defect_names=complex_defect_names,
                   included=included,
                   excluded=excluded,
                   displacement_distance=displacement_distance,
                   symprec=symprec,
                   angle_tolerance=angle_tolerance,
                   oxidation_states=oxidation_states,
                   electronegativity=electronegativity,
                   interstitials_yaml=interstitials_yaml,
                   complex_defect_yaml=complex_defect_yaml)

    # @property
    # def are_atoms_perturbed(self) -> bool:
    #     return self.cutoff > 1e-5

    @classmethod
    def from_basic_settings(cls,
                            structure: Structure,
                            transformation_matrix: list,
                            cell_multiplicity: int,
                            oxidation_states: Optional[dict] = None,
                            dopants: list = None,
                            is_antisite: bool = True,
                            en_diff: float = ELECTRONEGATIVITY_DIFFERENCE,
                            included: Optional[list] = None,
                            excluded: Optional[list] = None,
                            displacement_distance: float
                            = DISPLACEMENT_DISTANCE,
                            symprec: float = SYMMETRY_TOLERANCE,
                            angle_tolerance: float = ANGLE_TOL,
                            interstitial_sites: list = None,
                            complex_defect_names: list = None,
                            interstitials_yaml: str = "interstitials.yaml",
                            complex_defect_yaml: str = "complex_defects.yaml",
                            ) -> "DefectInitialSetting":
        """ Generates object with some default settings.

        Args:
            structure (Structure/IStructure):
                Structure used for supercell defect calculations
            transformation_matrix (list):
                Matrix used for expanding the unitcell.
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
            symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
            interstitial_sites (list):
                Names of interstitial sites written in interstitial.yaml file.
            complex_defect_names (list):
                Names of complex defects written in complex_defect.yaml file.

        """
        dopants = dopants or []
        interstitial_sites = interstitial_sites or "all"
        complex_defect_names = complex_defect_names or []

        # Here, the structure is sorted by elements.
        sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
        symmetrized_structure = sga.get_symmetrized_structure()
        equiv_sites = symmetrized_structure.equivalent_sites
        sorted_structure = \
            Structure.from_sites(reduce(lambda a, b: a + b, equiv_sites))

        sga = SpacegroupAnalyzer(sorted_structure, symprec, angle_tolerance)
        sym_dataset = sga.get_symmetry_dataset()

        species = [str(s) for s in sorted_structure.species]
        # unique_justseen https://docs.python.org/ja/3/library/itertools.html
        # ["H", "H", "O", "O", "H"] -> ['H', 'O', 'H']
        host_element_list = \
            list(map(next, map(operator.itemgetter(1), groupby(species, None))))
        element_list = host_element_list + dopants

        electronegativity = {e: get_electronegativity(e) for e in element_list}
        oxidation_states = \
            get_oxidation_states(dopants, oxidation_states, structure)

        # num_irreducible_sites["Mg"] = 2 means Mg has 2 inequivalent sites.
        num_irreducible_sites = defaultdict(int)

        # construct IrreducibleSite class objects
        irreducible_sites = []
        # Initialize the last index, which must start from -1 to begin with 0.
        last_index = -1
        for i, equiv_site in enumerate(equiv_sites):
            # set element name of equivalent site
            element = equiv_site[0].species_string
            # need to omit sign and numbers, eg Zn2+ -> Zn
            element = ''.join([i for i in element if i.isalpha()])

            # increment number of inequivalent sites for element
            num_irreducible_sites[element] += 1

            first_index = last_index + 1
            last_index = last_index + len(equiv_site)
            # np.array must be converted to list, being consistent with arg type
            repr_coords = list(equiv_site[0].frac_coords)
            irreducible_name = element + str(num_irreducible_sites[element])
            wyckoff = sym_dataset["wyckoffs"][first_index]
            site_symmetry = sym_dataset["site_symmetry_symbols"][first_index]

            # get the minimum distance from the representative coordinates.
            min_dist = min_distance_from_coords(structure, repr_coords)
            cutoff = round(min_dist * CUTOFF_FACTOR, 2)
            coordination_distances = \
                get_coordination_distances(structure=sorted_structure,
                                           atom_index=first_index,
                                           cutoff=cutoff)

            irreducible_sites.append(IrreducibleSite(irreducible_name,
                                                     element,
                                                     first_index,
                                                     last_index,
                                                     repr_coords,
                                                     wyckoff,
                                                     site_symmetry,
                                                     cutoff,
                                                     coordination_distances))

        # List of inserted and removed atoms, e.g., [["Mg, "O"], ...]
        antisite_configs = []
        if is_antisite is True:
            for elem1, elem2 in permutations(host_element_list, 2):
                if elem1 == elem2:
                    continue
                elif elem1 in oxidation_states and elem2 in oxidation_states:
                    if oxidation_states[elem1] * oxidation_states[elem2] > 0:
                        antisite_configs.append([elem1, elem2])
                else:
                    logger.warning(f"Oxidation states of {elem1} and/or " 
                                   f"{elem2} is not defined, so consideration"
                                   f"of antisite of its pair is skipped.")

            # for elem1, elem2 in permutations(host_element_list, 2):
            #     if elem1 == elem2:
            #         continue
            #     elif elem1 in electronegativity and elem2 in electronegativity:
            #         if abs(electronegativity[elem1]
            #                - electronegativity[elem2]) < en_diff:
            #             antisite_configs.append([elem1, elem2])
            #     else:
            #         logger.warning(f"Electronegativity of {elem1} and/or "
            #                        f"{elem2} is not defined")

        # List of inserted and removed atoms, e.g., [["Al", "Mg"], ...]
        dopant_configs = []
        for dopant in dopants:
            if dopant in host_element_list:
                logger.warning(f"Dopant {dopant} constitutes host.")
                continue
            for elem in host_element_list:
                if elem and dopant in electronegativity:
                    if abs(electronegativity[elem]
                           - electronegativity[dopant]) < en_diff:
                        dopant_configs.append([dopant, elem])
                else:
                    logger.warning(f"Electronegativity of {dopant} and/or "
                                   f"{elem} is not defined")

        return cls(structure=sorted_structure,
                   space_group_symbol=sga.get_space_group_symbol(),
                   transformation_matrix=transformation_matrix,
                   cell_multiplicity=cell_multiplicity,
                   irreducible_sites=irreducible_sites,
                   dopant_configs=dopant_configs,
                   antisite_configs=antisite_configs,
                   interstitial_site_names=interstitial_sites,
                   complex_defect_names=complex_defect_names,
                   included=included,
                   excluded=excluded,
                   displacement_distance=displacement_distance,
                   symprec=symprec,
                   angle_tolerance=angle_tolerance,
                   oxidation_states=oxidation_states,
                   electronegativity=electronegativity,
                   interstitials_yaml=interstitials_yaml,
                   complex_defect_yaml=complex_defect_yaml)

    def to_yaml_file(self, filename: str = "defect.yaml") -> None:
        """Dump to yaml file. """
        dumpfn(self.as_dict(), filename)

    def to_json_file(self, filename: str) -> None:
        """Dump to json file. """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def to(self,
           defect_in_file: str = "defect.in",
           poscar_file: str = "DPOSCAR") -> None:
        """Prints readable defect.in file and corresponding DPOSCAR file. """
        self._write_defect_in(defect_in_file)
        self.structure.to(fmt="poscar", filename=poscar_file)

    def _substituted_set(self,
                         removed_sites: List[IrreducibleSite],
                         inserted_element: Optional[str]) -> dict:
        """Helper method to create substituted and vacancy defect set. """
        defect_set = {}

        for rs in removed_sites:
            changes_of_num_elements = defaultdict(int)
            in_name = inserted_element if inserted_element else "Va"
            name = "_".join([in_name, rs.irreducible_name])
            symmetry = rs.site_symmetry
            structure = self.structure.copy()
            structure.remove_sites([rs.first_index])
            removed_atoms = [{"element": rs.element,
                              "index": rs.first_index,
                              "coords": rs.representative_coords}]
            changes_of_num_elements[rs.element] -= 1
            multiplicity = rs.multiplicity
            cutoff = rs.cutoff
            removed_oxi_state = self.oxidation_states[rs.element]

            if inserted_element:
                defect_type = DefectType.substituted
                atom = [{"element": inserted_element,
                         "coords": rs.representative_coords}]
                structure, inserted_atoms = insert_atoms(structure, atom)
                changes_of_num_elements[inserted_element] += 1
                inserted_oxi_state = self.oxidation_states[inserted_element]
            else:
                defect_type = DefectType.vacancy
                inserted_atoms = {}
                inserted_oxi_state = 0

            charges = default_charge_set(inserted_oxi_state - removed_oxi_state)

            defect_set[name] = \
                {"defect_type": defect_type,
                 "initial_structure": structure,
                 "removed_atoms": removed_atoms,
                 "inserted_atoms": inserted_atoms,
                 "changes_of_num_elements": dict(changes_of_num_elements),
                 "initial_site_symmetry": symmetry,
                 "charges": charges,
                 "multiplicity": multiplicity,
                 "cutoff": cutoff,
                 "center": rs.representative_coords}

        return defect_set

    def _inserted_set(self,
                      inserted_elements: Union[List[str], Tuple[str]]) -> dict:
        """Helper method to create interstitial defect set. """
        defect_set = {}
        for interstitial_name, i in self.interstitials.items():
            center = i.representative_coords
            symmetry = i.site_symmetry
            multiplicity = i.multiplicity
            cutoff = i.cutoff

            for ie in inserted_elements:
                changes_of_num_elements = defaultdict(int)
                # e.g., "i1_Mg"
                name = "_".join([ie, interstitial_name])
                extreme_charge = self.oxidation_states[ie]
                charges = default_charge_set(extreme_charge)
                changes_of_num_elements[ie] += 1
                inserted_atom = [{"element": ie,
                                  "coords": i.representative_coords}]
                structure, inserted_atoms = \
                    insert_atoms(self.structure, inserted_atom)

                defect_set[name] = \
                    {"defect_type": DefectType.interstitial,
                     "initial_structure": structure,
                     "removed_atoms": list(),
                     "inserted_atoms": inserted_atoms,
                     "changes_of_num_elements": dict(changes_of_num_elements),
                     "initial_site_symmetry": symmetry,
                     "charges": charges,
                     "multiplicity": multiplicity,
                     "cutoff": cutoff,
                     "center": center}

        return defect_set

    def _complex_set(self) -> dict:
        """Complex defect set. """
        defect_set = {}
        changes_of_num_elements = defaultdict(int)

        for name, complex_defect in self.complex_defects.items():
            structure = self.structure.copy()
            removed_atoms = []
            for index in complex_defect.removed_atom_indices:
                element = self.structure[index].species_string
                coords = list(self.structure[index].frac_coords)
                removed_atoms.append({"element": element,
                                      "index": index,
                                      "coords": coords})
                changes_of_num_elements[element] += -1

            structure.remove_sites(complex_defect.removed_atom_indices)
            removed_structure = structure.copy()

            # inserted_atoms contains indices of the inserted atoms in the
            # returned structure.
            structure, inserted_atoms = \
                insert_atoms(structure, complex_defect.inserted_atoms)
            for i in inserted_atoms:
                changes_of_num_elements[i["element"]] += -1

            coords = [i["coords"] for i in inserted_atoms + removed_atoms]

            # need to evaluate min_dist using the atom removed structure
            center = defect_center_from_coords(coords, removed_structure)
            min_dist = min_distance_from_coords(removed_structure, center)

            cutoff = round(min_dist * CUTOFF_FACTOR, 2)

            charges = default_charge_set(complex_defect.extreme_charge_state)

            defect_set[name] = \
                {"defect_type": DefectType.complex,
                 "initial_structure": structure,
                 "removed_atoms": removed_atoms,
                 "inserted_atoms": inserted_atoms,
                 "changes_of_num_elements": dict(changes_of_num_elements),
                 "initial_site_symmetry": complex_defect.point_group,
                 "charges": charges,
                 "multiplicity": complex_defect.multiplicity,
                 "cutoff": cutoff,
                 "center": center}

        return defect_set

    def make_defect_set(self,
                        keywords: Union[str, list, None] = None,
                        specified_defects: Optional[list] = None) -> None:
        """Return defect name list based on DefectInitialSetting object. """
        defects = {}
        # vacancies
        defects.update(self._substituted_set(
            removed_sites=self.irreducible_sites, inserted_element=None))

        # substituted + antisite
        configs = self.antisite_configs + self.dopant_configs
        for in_elem, out_elem in configs:
            removed_sites = [i for i in self.irreducible_sites
                             if out_elem == i.element]
            substituted = self._substituted_set(removed_sites=removed_sites,
                                                inserted_element=in_elem)
            defects.update(substituted)

        # interstitials
        inserted_elements = self.structure.symbol_set + tuple(self.dopants)
        defects.update(self._inserted_set(inserted_elements=inserted_elements))

        # complex defects
        defects.update(self._complex_set())

        defects = select_defects(defect_set=defects,
                                 keywords=keywords,
                                 included=self.included,
                                 excluded=self.excluded,
                                 specified_defects=specified_defects)

        self.defect_entries = []
        for name, defect in defects.items():
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
                        cutoff=defect["cutoff"],
                        distance=self.displacement_distance,
                        inserted_atom_indices=inserted_indices)

                defect["initial_structure"].set_charge(charge)
                perturbed_structure.set_charge(charge)

                self.defect_entries.append(
                    DefectEntry(name=name,
                                perturbed_initial_structure=perturbed_structure,
                                charge=charge,
                                neighboring_sites=neighboring_sites,
                                **defect))

    def _write_defect_in(self, defect_in_file: str = "defect.in") -> None:
        """Helper method to write down defect.in file.

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
            lines.append(f"         Cutoff radius: {site.cutoff}")

            neighbor_atom_distances = list()
            for k, v in site.coordination_distances.items():
                neighbor_atom_distances.append(
                   f"{k}: {' '.join([str(round(i, 2)) for i in v])}")
            neighbor_atom_distances = " ".join(neighbor_atom_distances)

            lines.append(f"          Coordination: {neighbor_atom_distances}")

            lines.append(f"      Equivalent atoms: "
                         f"{site.first_index}..{site.last_index}")

            lines.append("Fractional coordinates: "
                         "{0[0]:9.7f}  {0[1]:9.7f}  {0[2]:9.7f}".
                         format(site.representative_coords))

            electronegativity = self.electronegativity[site.element]
            lines.append(f"     Electronegativity: {electronegativity}")

            oxidation_state = self.oxidation_states[site.element]
            lines.append(f"       Oxidation state: {oxidation_state}")

            # Append "" for one blank line.
            lines.append("")

        if self.interstitial_site_names == "all":
            sites = "all"
        else:
            sites = ' '.join(self.interstitial_site_names)
        lines.append(f"Interstitials: {sites}")

        complex_defects = ' '.join(self.complex_defect_names)
        lines.append(f"Complex defects: {complex_defects}")

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
        lines.append(f"Symprec: {self.symprec}")
        lines.append(f"Angle tolerance: {self.angle_tolerance}")
        lines.append("")

        with open(defect_in_file, 'w') as defect_in:
            defect_in.write("\n".join(lines))

