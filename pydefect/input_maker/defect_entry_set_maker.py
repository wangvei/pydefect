# -*- coding: utf-8 -*-

import re

from pymatgen.core.periodic_table import Element

from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.core.defect_entry import DefectEntry
from pydefect.util.structure import perturb_neighboring_atoms, \
    defect_center_from_coords
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_int_from_string(x):
    """
    Returns joined integer number from a string.

    Args:
        x (str): a string
    """
    return int(''.join(i for i in x if i.isdigit()))


def parse_defect_name(defect_name):
    """
    Divides defect name by "_" to three part.
    E.g., "Va_Mg1_0" --> in_name="Va", out_name="Mg1", charge=0

    Args:
        defect_name (str): defect name defined in PyDefect, e.g., "Va_Mg1_2"
    """
    try:
        d = defect_name.split("_")
        in_name = d[0]
        out_name = d[1]
        charge = int(d[2])
        name = in_name + "_" + out_name
    except ValueError or IndexError:
        raise ValueError("Defect name {} is not proper.".format(defect_name))

    if not re.match(r'^[a-xA-Z]+[0-9]+$', out_name):
        raise ValueError("Defect name {} is not proper.".format(defect_name))
    return name, in_name, out_name, charge


def log_is_being_removed(name):
    """
    Shows the message.
    Args:
        name (str): a string
    """
    logger.warning("{:>10} is being removed.".format(name))


def log_already_exist(name):
    """
    Shows the message.
    Args:
        name (str): a string
    """
    logger.warning("{:>10} already exists, so nothing is done.".format(name))


def log_is_being_constructed(name):
    """
    Shows the message.
    Args:
        name (str): a string
    """
    logger.warning("{:>10} is being constructed.".format(name))


def is_name_matched(name, keywords):
    """
    Returns if name is matched by the selected_keywords.
     Args:
        name (str): Target name.
        keywords (str/list): Keywords used for checking if name is selected.

    When the following type names are given, constructs a set of defects.
        "Va"    --> A set of all the vacancies.
        "_i"     --> A set of all the interstitials.
        "Va_O"  --> A set of all the oxygen vacancies
        "Va_O[0-9]_0" --> All the oxygen vacancies in neutral charge states
        "Va_O1" --> A set of oxygen vacancies at O1 site
        "Mg_O"  --> A set of all the Mg-on-O antisite pairs.
        "Mg_O1" --> A set of Mg-on-O1 antisite pairs.

    When complete defect_name is given, constructs a particular defect.
        e.g., "Va_O1_2",  "Mg_O1_0"
    """

    try:
        if isinstance(keywords, str):
            keywords = [keywords]
        else:
            keywords = list(keywords)
    except TypeError:
        raise TypeError("The type of keywords {} is not list.".format(keywords))

    return any([re.search(p, name) for p in keywords])


def select_defect_names(name_set, keywords):
    """
    Returns names including one of keywords.
     Args:
        name_set (list): A set of names.
        keywords (str/list): Keywords used for checking if name is selected or not.
    """
    names = []

    for d in name_set:
        if is_name_matched(d, keywords):
            names.append(d)

    return list(set(names))


class DefectEntrySetMaker:
    """ Constructs a set of DefectEntry objects based on defect_initial_setting

    Parameters in use:
        in_pattern (str):
            Pattern for screening in_name
        out_pattern (str):
            Pattern for screening out_name
    """

    def __init__(self,
                 defect_initial_setting: DefectInitialSetting,
                 interstitial_in_file: str = None,
                 keywords: list = None,
                 particular_defects: list = None):
        """
        Args:
            defect_initial_setting (DefectInitialSetting):
                 DefectInitialSetting class object.
            interstitial_in_file (str):
                Name of interstitial.in file.
            keywords (list):
                Specify regular expression to narrow the defects considered.
            particular_defects (list):
                Specifies particular defect to be considered.
        """
        self.perfect_structure = defect_initial_setting.structure
        self.irreducible_sites = defect_initial_setting.irreducible_sites
        self.cutoff = defect_initial_setting.cutoff
        self.displacement_distance = \
            defect_initial_setting.displacement_distance
        self.cell_multiplicity = defect_initial_setting.cell_multiplicity

        if interstitial_in_file:
            from pydefect.core.interstitial_site import InterstitialSites
            self.interstitials = \
                InterstitialSites.from_interstitial_in(interstitial_in_file)
        else:
            self.interstitials = None

        if particular_defects:
            defect_name_set = particular_defects
        else:
            defect_name_set = defect_initial_setting.make_defect_name_set()
            if keywords:
                defect_name_set = select_defect_names(defect_name_set, keywords)

        self.defect_entries = \
            [self.make_defect_entry(d) for d in defect_name_set]

    def make_defect_entry(self, defect_name):
        """
        Constructs a single DefectEntry class object from a given defect_name.

        Args:
            defect_name (str):
                defect name in PyDefect manner, e.g., "Va_Mg2_-2".

        Parameters in use:
            in_name" (str):
                Inserted element name. "Va" is used for vacancies.
            out_name" (str):
                Removed site name. "in", where n is an integer, is used for
                interstitials. E.g., "i1".
            charge (int):
                DefectSupercell charge state
            removed_atom_index (int):
                Removed atom index from the perfect structure.
        """

        defect_structure = self.perfect_structure.copy()
        name, in_name, out_name, charge = parse_defect_name(defect_name)

        changes_of_num_elements = {}
        # -------------------- analyze out_name --------------------------------
        removed_atoms = {}
        # interstitial
        if re.match(r"^i[a-xA-Z0-9]+$", out_name):
            for i in self.interstitials:
                if out_name == i.site_name:
                    defect_coords = i.representative_coords
                    initial_site_symmetry = i.site_symmetry
                    num_equiv_sites = \
                        self.cell_multiplicity * i.symmetry_multiplicity
                    break
            else:
                raise IndexError("{} in {} is not in interstitial.in."
                                 .format(out_name, defect_name))

        else:
            for i in self.irreducible_sites:
                if out_name == i.irreducible_name:
                    changes_of_num_elements[i.element] = -1
                    removed_index = i.first_index - 1
                    defect_coords = i.representative_coords
                    removed_atoms[removed_index] = i.representative_coords
                    initial_site_symmetry = i.site_symmetry
                    num_equiv_sites = i.num_atoms
                    break
            else:
                raise IndexError("{} in {} is improper.".format(out_name,
                                                                defect_name))

            try:
                defect_structure.remove_sites([removed_index])
            except IndexError:
                print("{} in {} is improper.".format(out_name, defect_name))

        # -------------------- analyze in_name ---------------------------------
        inserted_atoms = []

        # This block must be following analyzing out_name because
        # defect coordinates is needed when inserting an in_name atom.
        if in_name == "Va":
            pass
        elif Element.is_valid_symbol(in_name):
            # There may be multiple irreducible sites for inserted element,
            # e.g., Mg1 and Mg2, element of in_name is inserted to just before
            # the same elements, otherwise to the 1st index.
            changes_of_num_elements[in_name] = 1
            if in_name in defect_structure.symbol_set:
                inserted_index = \
                    min(defect_structure.indices_from_symbol(in_name))
                inserted_atoms.append(inserted_index)
            else:
                inserted_index = 0
                inserted_atoms = [0]
            defect_structure.insert(inserted_index, in_name, defect_coords)
        else:
            raise ValueError("{} in {} is improper.".format(out_name,
                                                            defect_name))

        if self.cutoff:
            inserted_atom_coords = list([self.perfect_structure.frac_coords[k]
                                         for k in inserted_atoms])
            removed_atom_coords = list(removed_atoms.values())
            defect_coords = inserted_atom_coords + removed_atom_coords
            center = defect_center_from_coords(defect_coords,
                                               self.perfect_structure)

            perturbed_defect_structure, perturbed_sites = \
                perturb_neighboring_atoms(defect_structure, center,
                                          self.cutoff,
                                          self.displacement_distance)
        else:
            perturbed_defect_structure = defect_structure
            perturbed_sites = []

        defect_structure.set_charge(charge)
        perturbed_defect_structure.set_charge(charge)

        return DefectEntry(name=name,
                           initial_structure=defect_structure,
                           perturbed_initial_structure=perturbed_defect_structure,
                           removed_atoms=removed_atoms,
                           inserted_atoms=inserted_atoms,
                           changes_of_num_elements=changes_of_num_elements,
                           charge=charge,
                           initial_site_symmetry=initial_site_symmetry,
                           perturbed_sites=perturbed_sites,
                           num_equiv_sites=num_equiv_sites)
