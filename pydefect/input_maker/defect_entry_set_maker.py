# -*- coding: utf-8 -*-

from copy import deepcopy
import re

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure

from pydefect.core.defect_entry import DefectEntry
from pydefect.util.structure import perturb_neighbors

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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
    except ValueError:
        print("DefectSupercell {} is improper.".format(defect_name))

    if not re.match(r'^[a-xA-Z]+[0-9]+$', out_name):
        raise ValueError("Defect name {} is not proper.".format(defect_name))
    return in_name, out_name, charge


def print_is_being_removed(name):
    """
    Shows the message.
    Args:
        name (str): a string
    """
    print("{:>10} is being removed.".format(name))


def print_already_exist(name):
    """
    Shows the message.
    Args:
        name (str): a string
    """
    print("{:>10} already exists, so nothing is done.".format(name))


def print_is_being_constructed(name):
    """
    Shows the message.
    Args:
        name (str): a string
    """
    print("{:>10} is being constructed.".format(name))


def is_name_selected(name, keywords):
    """
    Returns if name is selected by selected_keywords.
     Args:
        name (str): Target name.
        keywords (list): Keywords used for checking if name is selected or not.

    When the following type names are given, constructs a set of defects.
        "Va"    --> A set of all the vacancies.
        "_i"     --> A set of all the interstitials.
        "Va_O"  --> A set of all the oxygen vacancies
        "Va_O1" --> A set of oxygen vacancies at O1 site
        "Mg_O"  --> A set of all the Mg-on-O antisite pairs.
        "Mg_O1" --> A set of Mg-on-O1 antisite pairs.

    When complete defect_name is given, constructs a particular defect.
        e.g., "Va_O1_2",  "Mg_O1_0"
    """

    if type(keywords) is not list:
        raise TypeError("The type of keywords is not list.")

    for p in keywords:
        if re.search(p, name):
            return True

    return False


def select_defect_names(name_set, keywords):
    """
    Returns names selected by selected_keywords.
     Args:
        name_set (list): A set of names.
        keywords (list): Keywords used for checking if name is selected or not.
    """
    names = []

    for d in name_set:
        if is_name_selected(d, keywords):
            names.append(d)

    return list(set(names))


class DefectEntrySetMaker:
    """
    Abstract class that is subclassed by a particular first-principles code
    implementation. Constructs a set of DefectSupercell class objects based on
    given oxidation states.

    Args:
        defect_initial_setting (DefectInitialSetting):
            DefectInitialSetting class object.
        keywords (list):
            Specify regular expression to narrow the defects considered.
        particular_defects (list):
            Specifies particular defect to be considered.

    Parameters in use:
        in_pattern (str):
            Pattern for screening in_name
        out_pattern (str):
            Pattern for screening out_name
    """
    def __init__(self, defect_initial_setting, keywords=None,
                 particular_defects=None):

        self._structure = defect_initial_setting.structure
        self._irreducible_sites = defect_initial_setting.irreducible_sites
        self._interstitial_coords = defect_initial_setting.interstitial_coords
        self._cutoff = defect_initial_setting.cutoff
        self._distance = defect_initial_setting.distance

        if particular_defects:
            defect_name_set = particular_defects
        else:
            defect_all_name_set = defect_initial_setting.make_defect_name_set()
            if keywords:
                defect_name_set = \
                    select_defect_names(defect_all_name_set, keywords)
                if "perfect" in keywords:
                    defect_name_set.append("perfect")
            else:
                defect_all_name_set.append("perfect")
                defect_name_set = defect_all_name_set

        self.defect_entries = []
        for d in defect_name_set:
            self.defect_entries.append(self.make_defect_entry(d))

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

        defect_structure = deepcopy(self._structure)
        in_name, out_name, charge = parse_defect_name(defect_name)
        name = in_name + "_" + out_name

        # TODO: add initial_symmetry and multiplicity for interstitials.
        initial_symmetry = "1"
        multiplicity = 1

        changes_of_num_elements = {}
        # -------------------- analyze out_name --------------------------------
        removed_atoms = {}
        # interstitial
        if re.match(r'^i[0-9]+$', out_name):
            interstitial_index = get_int_from_string(out_name)
            try:
                defect_coords = self._interstitial_coords[interstitial_index - 1]
            except ValueError:
                print("#{} interstitial not defined".format(interstitial_index))
        else:
            for i in self._irreducible_sites:
                if out_name == i.irreducible_name:
                    changes_of_num_elements[i.element] = -1
                    removed_index = i.first_index - 1
                    defect_coords = i.repr_coords
                    removed_atoms[removed_index] = i.repr_coords
                    initial_symmetry = i.site_symmetry
                    multiplicity = i.num_atoms
                    break
            try:
                defect_structure.remove_sites([removed_index])
            except ValueError:
                print("{} in {} is improper.".format(out_name, defect_name))

        # -------------------- analyze in_name ---------------------------------
        inserted_atoms = []

        # This block needs to be run after finishing analyze_out_name because
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

        if self._cutoff:
            perturbed_defect_structure, perturbed_sites = \
                perturb_neighbors(defect_structure, defect_coords,
                                  self._cutoff, self._distance)
        else:
            perturbed_defect_structure = defect_structure
            perturbed_sites = []

        defect_structure.set_charge(charge)
        perturbed_defect_structure.set_charge(charge)

        return DefectEntry(name, defect_structure, perturbed_defect_structure,
                           removed_atoms, inserted_atoms,
                           changes_of_num_elements, charge, initial_symmetry,
                           multiplicity, perturbed_sites)

