# -*- coding: utf-8 -*-
from pymatgen import Composition, Structure
from pymatgen.io.vasp import Potcar, Incar

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


def element_diff_from_poscar_files(poscar1, poscar2):
    """
    Returns a dict of change of numbers of elements from poscar2 to poscar1
    For defect calculations, poscar2 should be "perfect".
    """
    c1 = Composition(
        Structure.from_file(poscar1).composition, allow_negative=True)
    c2 = Composition(
        Structure.from_file(poscar2).composition, allow_negative=True)
    c_diff = c1 - c2

    return {str(e): int(c_diff[e]) for e in c_diff}


def get_num_electrons_from_potcar(potcar, nions, charge=0):
    """
    Returns the number of electrons from POTCAR, number of ions, and charge
    state.
    """
    p = Potcar.from_file(potcar)
    # check only the number of ions written in potcar and nions.
    if not len(p) == len(nions):
        raise ValueError("Size of elements in POTCAR file is different")

    return sum([v.nelectrons * nions[i] for i, v in enumerate(p)]) - charge


def get_defect_charge_from_vasp(structure, potcar="POTCAR", incar="INCAR"):
    """
    Returns the defect charge by comparing nion, number of electrons in POTCAR,
    and NELECT in INCAR.
    """
    try:
        num_elect_incar = Incar.from_file(incar)["NELECT"]
    # return 0 if NELECT is absent in INCAR
    except KeyError:
        return 0

    potcar = Potcar.from_file(potcar)

    num_elect_neutral = 0

    for i, e in enumerate(structure.composition):
        if potcar[i].element != str(e):
            raise ValueError("The sequence of elements in POTCAR and Structure "
                             "is different.")

        nelect = potcar[i].nelectrons
        nions = structure.composition[e]
        num_elect_neutral += nelect * nions

    # charge is minus of difference of the electrons
    return int(num_elect_neutral - num_elect_incar)