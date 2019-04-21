# -*- coding: utf-8 -*-

import re

from pymatgen.core.periodic_table import Element

from pydefect.input_maker.defect_initial_setting import DefectInitialSetting, SimpleDefectName
from pydefect.core.defect_entry import DefectEntry
from pydefect.util.structure_tools import perturb_neighboring_atoms, \
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


def select_defect_names(name_set, keywords):
    """
    Returns names including one of keywords.
     Args:
        name_set (list): A set of names.
        keywords (str/list): Keywords determining if name is selected or not.
    """
    names = []

    for name in name_set:
        if name.is_name_matched(keywords):
            names.append(name)

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
                 keywords: list = None,
                 particular_defects: list = None):
        """
        Args:
            defect_initial_setting (DefectInitialSetting):
                 DefectInitialSetting class object.
            keywords (list):
                Specify regular expression to narrow the defects considered.
            particular_defects (list of strings):
                Specifies particular defect to be considered.
        """
        self.perfect_structure = defect_initial_setting.structure
        self.irreducible_sites = defect_initial_setting.irreducible_sites
        self.cutoff = defect_initial_setting.cutoff
        self.displacement_distance = \
            defect_initial_setting.displacement_distance
        self.cell_multiplicity = defect_initial_setting.cell_multiplicity
        self.interstitials = defect_initial_setting.interstitials
        self.is_displaced = defect_initial_setting.are_atoms_perturbed

        if particular_defects:
            defect_name_set = \
                [SimpleDefectName.from_str(i) for i in particular_defects]
        else:
            defect_name_set = defect_initial_setting.make_defect_name_set()
            if keywords:
                defect_name_set = select_defect_names(defect_name_set, keywords)

        self.defect_entries = \
            [self.make_defect_entry(d) for d in defect_name_set]

    def make_defect_entry(self,
                          defect_name: SimpleDefectName):
        """ Constructs a single DefectEntry object from a given defect_name.

        Args:
            defect_name (SimpleDefectName):

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

        changes_of_num_elements = {}
        # -------------------- analyze out_site --------------------------------
        removed_atoms = {}
        # interstitial
        if defect_name.is_interstitial:
            for i in self.interstitials:
                if defect_name.out_site == i.site_name:
                    defect_coords = i.representative_coords
                    initial_site_symmetry = i.site_symmetry
                    num_equiv_sites = i.symmetry_multiplicity
                    break
            else:
                raise IndexError("{} in {} is not in interstitial.in."
                                 .format(defect_name.out_site, defect_name))

        else:
            for i in self.irreducible_sites:
                if defect_name.out_site == i.irreducible_name:
                    changes_of_num_elements[i.element] = -1
                    removed_index = i.first_index - 1
                    defect_coords = i.representative_coords
                    removed_atoms[removed_index] = i.representative_coords
                    initial_site_symmetry = i.site_symmetry
                    num_equiv_sites = i.num_atoms
                    break
            else:
                raise IndexError("{} in {} is improper.".
                                 format(defect_name.out_site, defect_name))

            try:
                defect_structure.remove_sites([removed_index])
            except IndexError:
                print("{} in {} is improper.".format(defect_name.out_site,
                                                     defect_name))

        # -------------------- analyze in_atom ---------------------------------
        inserted_atoms = []

        # This block must be following analyzing out_name because
        # defect coordinates is needed when inserting an in_name atom.
        if defect_name.is_vacancy:
            pass
        elif Element.is_valid_symbol(defect_name.in_atom):
            # There may be multiple irreducible sites for inserted element,
            # e.g., Mg1 and Mg2, element of in_atom is inserted to just before
            # the same elements, otherwise to the 1st index.
            changes_of_num_elements[defect_name.in_atom] = 1
            if defect_name.in_atom in defect_structure.symbol_set:
                inserted_index = min(
                    defect_structure.indices_from_symbol(defect_name.in_atom))
                inserted_atoms.append(inserted_index)
            else:
                inserted_index = 0
                inserted_atoms = [0]
            defect_structure.insert(inserted_index, defect_name.in_atom,
                                    defect_coords)
        else:
            raise ValueError("{} in {} is improper.".format(defect_name.in_atom,
                                                            defect_name))

        if self.is_displaced:
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
            perturbed_defect_structure = defect_structure.copy()
            perturbed_sites = []

        defect_structure.set_charge(defect_name.charge)
        perturbed_defect_structure.set_charge(defect_name.charge)

        return DefectEntry(
            name=defect_name.name_str,
            initial_structure=defect_structure,
            perturbed_initial_structure=perturbed_defect_structure,
            removed_atoms=removed_atoms,
            inserted_atoms=inserted_atoms,
            changes_of_num_elements=changes_of_num_elements,
            charge=defect_name.charge,
            initial_site_symmetry=initial_site_symmetry,
            perturbed_sites=perturbed_sites,
            num_equiv_sites=num_equiv_sites)
