# -*- coding: utf-8 -*-

import numpy as np
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Potcar, Incar, Procar

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def element_diff_from_structures(structure: Structure,
                                 ref_structure: Structure) -> dict:
    """
    Returns dict of change of numbers of elements from structure2 to structure1
    For defect calculations, structure2 should be "perfect".
    """
    to_comp = Composition(structure.composition, allow_negative=True)
    from_comp = Composition(ref_structure.composition, allow_negative=True)
    comp_diff = to_comp - from_comp

    return {str(e): int(diff) for e, diff in comp_diff.items()}


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


def calc_participation_ratio(procar: Procar,
                             spin: Spin,
                             band_index: int,
                             atom_indices: list) -> float:
    """ Returns sum of participation ratios at atom_indices sites

    The PROCAR data of the form below. It should VASP uses 1-based indexing,
    but all indices are converted to 0-based here.::
        { spin: nd.array accessed with (k-point index, band index,
                                          ion index, orbital index) }

    Note that the k-point weight is not considered, so all the k-points are
    treated equally.

    Args:
        procar: Pymatgen Procar class object.
        spin: Spin object
        band_index: Index of band that begins from zero
        atom_indices:
            List of atom indices, for which the participation ratio is
            calculated.

    Return (float):
        float of the participation ratio.
    """
    # sum along k-point and orbital
    sum_per_atom = np.sum(procar.data[spin][:, band_index, :, :], axis=(0, 2))

    return np.sum(sum_per_atom[atom_indices]) / np.sum(sum_per_atom)


def calc_orbital_character(procar: Procar,
                           structure: Structure,
                           spin: Spin,
                           band_index: int,
                           kpoint_index: int):
    """ Method returning a dictionary of projections on elements.

    Args:
        procar (Procar):
            Procar object to be parsed.
        structure (Structure):
            Input structure used for extracting symbol_set
        spin (Spin):
        band_index (int):
        kpoint_index (int):

    Returns (dict):
        a dictionary in the {Element:{s: value, p: value, d: value}
        where s, p, and d are OrbitalType in pymatgen.
    """
    orbital_components = {}

    def projection_sum(atom_indices: tuple, first: int, last: int):
        end = last + 1
        procar_sum = np.sum(
            procar.data[spin][kpoint_index, band_index, atom_indices, first:end]
        )
        return float(procar_sum)

    for element in structure.symbol_set:
        # get list of index
        indices = structure.indices_from_symbol(element)
        orbital_components[element] = \
            {"s": round(projection_sum(indices, 0, 0), 3),
             "p": round(projection_sum(indices, 1, 3), 3),
             "d": round(projection_sum(indices, 4, 8), 3)}
        try:
            orbital_components[element]["f"] = \
                round(projection_sum(indices, 9, 16), 3)
        except KeyError:
            pass

    return orbital_components


def calc_orbital_difference(orbital_1: dict, orbital_2: dict) -> float:
    """ Calculate absolute difference of between two orbitals

    If an element exists in one orbital, its difference corresponds to the
    sum of the orbital component.

    Args:
        orbital_1 (dict):
            {"Mg": {"s": 0.1, ...}, "O": {...},..
        orbital_2 (dict):
            {"Al": {"s": 0.02, ...}, "Mg": {"s": 0.1, ...}, "O": {...},..

    Return (float):
        Sum of difference
    """
    element_set = set(list(orbital_1.keys()) + list(orbital_2.keys()))

    difference = 0
    for e in element_set:
        for o in "s", "p", "d", "f":

            try:
                first = orbital_1[e][o]
            except KeyError:
                first = 0

            try:
                second = orbital_2[e][o]
            except KeyError:
                second = 0

            difference += abs(first - second)

    return difference

