# -*- coding: utf-8 -*-

from collections import defaultdict
from typing import Union, List

from pymatgen.core.periodic_table import Element

from pydefect.core.defect_entry import DefectType, DefectEntry
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.core.defect_name import SimpleDefectName
from pydefect.util.structure_tools import perturb_neighboring_atoms, \
    defect_center_from_coords
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import first_appearance_index


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def get_int_from_string(x: str) -> int:
    """ Returns joined integer number from a string.

    Args:
        x (str): a string
    """
    return int(''.join(i for i in x if i.isdigit()))


def log_is_being_removed(name: str) -> None:
    """ Shows the message.

    Args:
        name (str): a string
    """
    logger.warning("{:>10} is being removed.".format(name))


def log_already_exist(name: str) -> None:
    """ Shows the message.

    Args:
        name (str): a string
    """
    logger.warning("{:>10} already exists, so nothing is done.".format(name))


def log_is_being_constructed(name: str) -> None:
    """ Shows the message.

    Args:
        name (str): a string
    """
    logger.warning("{:>10} is being constructed.".format(name))


def select_defect_names(name_set: List[SimpleDefectName],
                        keywords: Union[str, list],
                        return_str: bool = False) -> List[str]:
    """ Returns names including one of keywords.

    Args:
        name_set (list):
            A set of names (SimpleDefectName).
        keywords (str/list):
            Keywords determining if name is selected or not.
        return_str (bool):
            Whether to return string instead of SimpleDefectName object.

    Return:
         list of names set.
    """
    names = []
    for name in name_set:
        if isinstance(name, str):
            name = SimpleDefectName.from_str(name)

        if name.is_name_matched(keywords):
            if return_str:
                names.append(str(name))
            else:
                names.append(name)

    return names


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
                 keywords: Union[list, None] = None,
                 particular_defects: Union[list, None] = None):
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
        self.is_displaced = defect_initial_setting.are_atoms_perturbed
        self.interstitials = defect_initial_setting.interstitials

        defect_name_set = defect_initial_setting.make_defect_name_set()
        if particular_defects:
            defect_name_set = \
                [SimpleDefectName.from_str(i) for i in particular_defects]
        else:
            if keywords:
                defect_name_set = \
                    select_defect_names(defect_name_set, keywords)

        self.defect_entries = \
            [self.make_defect_entry(d) for d in defect_name_set]

    def make_defect_entry(self,
                          defect_name: SimpleDefectName) -> DefectEntry:
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
        defect_type = DefectType.from_simple_defect_name(defect_name)

        defect_structure = self.perfect_structure.copy()

        changes_of_num_elements = defaultdict(int)
        # -------------------- analyze out_site -------------------------------
        removed_atoms = {}
        # interstitial
        if defect_name.is_interstitial:
            for site_name, i in self.interstitials.items():
                if defect_name.out_site == site_name:
                    defect_coords = i.representative_coords
                    initial_site_symmetry = i.site_symmetry
                    num_equiv_sites = i.symmetry_multiplicity
                    break
            else:
                raise IndexError(f"{defect_name.out_site} in {defect_name} "
                                 f"is not in defect.in.")

        else:
            for i in self.irreducible_sites:
                if defect_name.out_site == i.irreducible_name:
                    changes_of_num_elements[i.element] -= 1
                    removed_index = i.first_index - 1
                    defect_coords = i.representative_coords
                    removed_atoms[removed_index] = i.representative_coords
                    initial_site_symmetry = i.site_symmetry
                    num_equiv_sites = i.num_atoms
                    break
            else:
                raise IndexError(
                    f"{defect_name.out_site} in {defect_name} is improper.")

            try:
                defect_structure.remove_sites([removed_index])
            except IndexError:
                print(f"{defect_name.out_site} in {defect_name} is improper.")

        # -------------------- analyze in_atom --------------------------------
        inserted_atoms = {}
        # This block must be after analyzing out_name as # defect coordinates
        # are needed when inserting an in_name atom.
        if defect_name.is_vacancy:
            pass
        elif Element.is_valid_symbol(defect_name.in_atom):
            # There may be multiple irreducible sites for inserted element,
            # e.g., Mg1 and Mg2, element of in_atom is inserted to just before
            # the same elements, otherwise to the 1st index.
            changes_of_num_elements[defect_name.in_atom] += 1
            inserted_index = first_appearance_index(defect_structure,
                                                    defect_name.in_atom)
            inserted_atoms[inserted_index] = defect_coords
            defect_structure.insert(inserted_index, defect_name.in_atom,
                                    defect_coords)
        else:
            raise ValueError(
                f"{defect_name.in_atom} in {defect_name} is improper.")

        inserted_atom_coords = list([defect_structure.frac_coords[k]
                                     for k in inserted_atoms])
        removed_atom_coords = list(removed_atoms.values())
        all_defect_coords = inserted_atom_coords + removed_atom_coords

        center = defect_center_from_coords(all_defect_coords,
                                           self.perfect_structure)

        # By default, neighboring atoms are perturbed.
        # If one wants to avoid it, set displacement_distance = 0
        perturbed_defect_structure, neighboring_sites = \
            perturb_neighboring_atoms(
                structure=defect_structure,
                center=center,
                cutoff=self.cutoff,
                distance=self.displacement_distance,
                inserted_atom_indices=list(inserted_atoms.keys()))

        defect_structure.set_charge(defect_name.charge)
        perturbed_defect_structure.set_charge(defect_name.charge)

        return DefectEntry(
            name=defect_name.name_str,
            defect_type=defect_type,
            initial_structure=defect_structure,
            perturbed_initial_structure=perturbed_defect_structure,
            removed_atoms=removed_atoms,
            inserted_atoms=inserted_atoms,
            changes_of_num_elements=dict(changes_of_num_elements),
            charge=defect_name.charge,
            initial_site_symmetry=initial_site_symmetry,
            cutoff=self.cutoff,
            neighboring_sites=neighboring_sites,
            num_equiv_sites=num_equiv_sites)
