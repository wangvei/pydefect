# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod
from copy import deepcopy
import re

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure

from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


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

    if not re.match(r'^[a-xA-Z]+[1-9]+$', out_name):
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
        pattern = r"" + re.escape(p)
        if re.search(pattern, name):
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


class DefectMaker:
    """
    Constructs a single DefectEntry class object from a given defect_name.

    Args:
        defect_name (str):
            defect name in PyDefect manner, e.g., "Va_Mg2_-2".
        structure (Structure):
            Structure class object for the *perfect* supercell.
        irreducible_sites ([IrreducibleSite])
        interstitial_coords (Nx3 list): coordinates of interstitial sites,
                                         e.g., [[0, 0, 0], [0.1, 0.1, 0.1], ..]

    Parameters in use:
        in_name" (str): Inserted element name. "Va" is used for vacancies.
        out_name" (str): Removed site name. "in", where n is an integer,
                         is used for interstitials. E.g., "i1".
        charge (int): DefectSupercell charge state
        removed_atom_index (int): removed atom index from the perfect structure.
    """
    def __init__(self, defect_name, structure, irreducible_sites,
                 interstitial_coords):

        # deepcopy is required for modifying the original structure.
        defect_structure = deepcopy(structure)
        in_name, out_name, charge = parse_defect_name(defect_name)
        name = in_name + "_" + out_name

        changes_of_num_elements = {}
        # -------------------- analyze out_name --------------------------------
        removed_atoms = {}
        # interstitial
        if re.match(r'^i[0-9]+$', out_name):
            interstitial_index = get_int_from_string(out_name)
            try:
                defect_coords = interstitial_coords[interstitial_index - 1]
            except ValueError:
                print("#{} interstitial not defined".format(interstitial_index))
        else:
            for i in irreducible_sites:
                if out_name == i.irreducible_name:
                    changes_of_num_elements[i.element] = -1
                    removed_index = i.first_index - 1
                    defect_coords = i.repr_coords
                    removed_atoms[removed_index] = i.repr_coords
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
        self.defect = DefectEntry(name, defect_structure, removed_atoms,
                                  inserted_atoms, changes_of_num_elements,
                                  charge)


class DefectInputSetMaker(metaclass=ABCMeta):
    """
    Abstract class that is subclassed by a particular first-principles code
    implementation. Constructs a set of DefectSupercell class objects based on
    given oxidation states.

    Args:
        defect_initial_setting (DefectInitialSetting):
            DefectInitialSetting class object.
        keywords (list):
            Specify a type of defects.
        particular_defects (list):
            It specifies (a) particular defect(s).

    Parameters in use:
        in_pattern (str): pattern for screening in_name
        out_pattern (str): pattern for screening out_name
    """

    def __init__(self, defect_initial_setting, keywords=None,
                 particular_defects=None, force_overwrite=False):

        self._defect_initial_setting = defect_initial_setting
        self._keywords = keywords
        self._particular_defects = particular_defects
        self._force_overwrite = force_overwrite

        if particular_defects:
            self._defect_name_set = particular_defects
        else:
            defect_all_name_set = defect_initial_setting.make_defect_name_set()
            if keywords:
                self._defect_name_set = \
                    select_defect_names(defect_all_name_set, keywords)
            else:
                defect_all_name_set.append("perfect")
                self._defect_name_set = defect_all_name_set

    def make_input(self):
        # Construct defect input files.
        for d in self._defect_name_set:
            if d == "perfect":
                self._make_perfect_input()
            else:
                self._make_defect_input(d)

    @abstractmethod
    def _make_perfect_input(self):
        pass

    @abstractmethod
    def _make_defect_input(self, defect_name):
        pass
