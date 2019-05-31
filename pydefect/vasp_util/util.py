# -*- coding: utf-8 -*-
import numpy as np
from collections import defaultdict

from pymatgen.core.composition import Composition
from pymatgen.core.sites import Element
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.io.vasp import Potcar, Incar, Procar

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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


def calc_participation_ratio(procar: Procar,
                             spin: Spin,
                             band_index: int,
                             kpoint_index: int,
                             atom_indices: list):
    """ Returns sum of participation ratios at atom_indices sites

    The PROCAR data of the form below. It should VASP uses 1-based indexing,
    but all indices are converted to 0-based here.::

        {
            spin: nd.array accessed with (k-point index, band index,
                                          ion index, orbital index)
        }
    """

    if spin not in procar.data.keys():
        spin = Spin.up

    # sum along k-point and orbital
    projected_to_atoms = \
        np.sum(procar.data[spin][kpoint_index], axis=2)[band_index]

    return np.sum(projected_to_atoms[atom_indices]) / np.sum(projected_to_atoms)


def calc_orbital_character(procar: Procar,
                           structure: Structure,
                           spin: Spin,
                           band_index: int,
                           kpoint_index: int):
    """ Method returning a dictionary of projections on elements.

    Args:
        structure (Structure): Input structure.

    Returns:
        a dictionary in the {Element:{s: value, p: value, d: value}
        where s, p, and d are OrbitalType in pymatgen.
    """
    d = defaultdict(dict)

    for name in structure.symbol_set:
        # get list of index
        indices = \
            [i for i, e in enumerate(structure.species) if e == Element(name)]
        d[name]["s"] = \
            np.sum(procar.data[spin][kpoint_index, band_index, indices, 0])
        d[name]["p"] = \
            np.sum(procar.data[spin][kpoint_index, band_index, indices, 1:4])
        d[name]["d"] = \
            np.sum(procar.data[spin][kpoint_index, band_index, indices, 4:9])
        try:
            d[name]["f"] = \
                np.sum(procar.data[spin][kpoint_index, band_index,
                       indices, 9:16])
        except:
            pass

    return d


def calc_orbital_similarity(orbital_1: defaultdict,
                            orbital_2: defaultdict):
    """
    """
    elements = set(list(orbital_1.keys()) + list(orbital_2.keys()))

    diff = 0
    for e in elements:
        for o in "s", "p", "d", "f":
            diff += abs(orbital_1[e][o] - orbital_2[e][o])

    return diff

