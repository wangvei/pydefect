# -*- coding: utf-8 -*-

import numpy as np
from collections import Iterable
from copy import deepcopy
from typing import Union

from pydefect.core.error_classes import TooLargeUnitcellError
from pymatgen.core.structure import Structure

from obadb.util.structure_handler import find_spglib_standard_conventional, \
    find_spglib_standard_primitive

from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def calc_isotropy(structure, multi):
    """
    Isotropy is defined as the mean absolute deviation of the lattice
    constants from the averaged lattice constant.
    Return
        isotropy (float):
    """
    abc = structure.lattice.abc
    super_abc = multi * abc
    average_abc = np.mean(super_abc)
    return round(np.sum(np.abs(super_abc - average_abc) / average_abc) / 3, 4)


class Supercell:
    def __init__(self,
                 structure: Structure,
                 trans_mat: Union[int, np.array],
                 multiplicity: int):
        """ Constructs a supercell based on a given multiplicity.

        Args:
            structure (Structure):
                Original structure to be expanded.
            trans_mat (3x3 np.array, 3 np.array or a scalar):
                The matrix to be used for expanding the unitcell.
            multiplicity (int):
                The size multiplicity of structure wrt the primitive cell.
        """
        trans_mat = \
            list(trans_mat) if isinstance(trans_mat, Iterable) else [trans_mat]

        if len(trans_mat) == 1:
            trans_mat_str = str(trans_mat)
        elif len(trans_mat) == 9:
            trans_mat_str = ' '.join([str(i) for i in trans_mat])
            trans_mat = np.reshape(trans_mat, (3, 3))
        elif len(trans_mat) == 3:
            trans_mat = np.array(trans_mat)
            if isinstance(trans_mat[0], list) and len(trans_mat[0]) == 3:
                trans_mat_str = ' '.join([str(int(i)) for i in trans_mat.flatten()])
            else:
                trans_mat_str = ' '.join([str(int(i)) for i in trans_mat])
        else:
            raise ValueError("Translation matrix: {} is not proper. 1, 3, or 9 "
                             "numbers  are accepted.")

        self.base_structure = structure
        s = structure * trans_mat
        self.structure = s.get_sorted_structure()
        self.multi = trans_mat
        self.isotropy = calc_isotropy(structure, trans_mat)
        self.num_atoms = self.structure.num_sites

        self.comment = 'trans_mat: ' + trans_mat_str + ', multi: ' + str(multiplicity) + ', isotropy: ' + \
                       str(self.isotropy) + '\n'

    def to_poscar(self, poscar_filename, uposcar_filename):
        self.base_structure.to(filename=uposcar_filename)

        poscar_str = self.structure.to(fmt="poscar").splitlines(True)
        poscar_str[0] = self.comment
        with open(poscar_filename, 'w') as fw:
            for line in poscar_str:
                fw.write(line)


class Supercells:
    def __init__(self,
                 structure: Structure,
                 is_conventional: bool = True,
                 max_num_atoms: int = 400,
                 min_num_atoms: int = 50,
                 isotropy_criterion: float = 0.12):
        """ Constructs a set of supercells satisfying a criterion.

        Args:
            structure (pmg structure class object):
                Unitcell structure
            is_conventional (bool):
                Conventional cell is expanded when True, otherwise primitive
                cell is expanded.
            max_num_atoms (int):
                Maximum number of atoms in the supercell.
            min_num_atoms (int):
                Minimum number of atoms in the supercell.
            isotropy_criterion (float):
                Criterion to judge if a supercell is isotropic or not.
                Isotropy is defined as the mean absolute deviation of the
                lattice constants from the averaged lattice constant.
                np.sum(np.abs(abc - average_abc) / average_abc) / 3
        """

        structure = structure.copy()

        self.is_conventional_based = False
        if is_conventional:
            primitive = find_spglib_standard_primitive(structure)[0]
            conventional = find_spglib_standard_conventional(structure)
            if conventional.num_sites == primitive.num_sites:
                logger.info("Primitive cell is same as is_conventional cell.")
            else:
                self.is_conventional_based = True

            unitcell = conventional
        else:
            primitive = find_spglib_standard_primitive(structure)[0]
            unitcell = primitive

        self.unitcell = unitcell.get_sorted_structure()
        abc = np.array(self.unitcell.lattice.abc)
        num_atoms_in_unitcell = self.unitcell.num_sites

        if max_num_atoms < num_atoms_in_unitcell:
            raise TooLargeUnitcellError("Number of atoms in the unitcell is "
                                        "smaller than the maximum number of "
                                        "atoms in the supercell")

        trans_mat = np.ones(3, dtype="int8")
        self.supercells = []

        for i in range(int(max_num_atoms / num_atoms_in_unitcell)):
            num_atoms = trans_mat.prod() * num_atoms_in_unitcell
            if num_atoms > max_num_atoms:
                break

            isotropy = calc_isotropy(self.unitcell, trans_mat)
            if isotropy < isotropy_criterion and num_atoms >= min_num_atoms:
                multi = int(self.unitcell.num_sites / primitive.num_sites * trans_mat.prod())
                self.supercells.append(Supercell(self.unitcell, trans_mat, multi))

            super_abc = trans_mat * abc
            # multi indices within 1.05a, where a is the shortest supercell
            # lattice length, are incremented.
            for j in range(3):
                if super_abc[j] / min(super_abc) < 1.05:
                    trans_mat[j] += 1

    @property
    def create_sorted_supercells_by_num_atoms(self):
        return sorted(deepcopy(self.supercells),
                      key=lambda x: (x.num_atoms, x.isotropy))

    @property
    def create_sorted_supercells_by_isotropy(self):
        return sorted(deepcopy(self.supercells),
                      key=lambda x: (x.isotropy, x.num_atoms))

    @property
    def create_smallest_supercell(self):
        return self.create_sorted_supercells_by_num_atoms[0]

    @property
    def create_most_isotropic_supercell(self):
        return self.create_sorted_supercells_by_isotropy[0]


