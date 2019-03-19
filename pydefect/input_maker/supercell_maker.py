# -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy

from pymatgen.core.structure import Structure

from obadb.util.structure_handler import find_spglib_standard_conventional, \
    find_spglib_standard_primitive

from pydefect.util.utils import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class Supercell:
    def __init__(self, structure, multi, comment=None):
        """
        Constructs a supercell based on a multi 3x3 matrix.
        Args:
            structure (Structure):
                Original structure to be expanded.
            multi (3x3 numpy array ,list, or a scalar):
                The matrix to be used for expanding the unitcell.
            comment (str):
                Any comment.
        """
        if len(multi) == 9:
            multi = np.reshape(multi, (3, 3))
        elif len(multi) == 3:
            multi = np.array(multi)

        s = structure * multi
        self._multi = multi
        self._supercell_structure = s.get_sorted_structure()
        self._isotropy = self.calc_supercell_isotropy(structure, multi)
        self._num_atoms = self._supercell_structure.num_sites

        if comment is None:
            self._comment = \
                self.supercell_comment(multi, self._isotropy)
        else:
            self._comment = comment

    @classmethod
    def from_poscar(cls, poscar, multi):
        structure = Structure.from_file(poscar)
        return cls(structure, multi)

    def to_poscar(self, filename):
        poscar_str = self.supercell_structure.to(fmt="poscar").splitlines(True)
        poscar_str[0] = self.comment
        with open(filename, 'w') as fw:
            for line in poscar_str:
                fw.write(line)

    @staticmethod
    def calc_isotropy(abc):
        """
        Isotropy is defined as the mean absolute deviation of the lattice
        constants from the averaged lattice constant.
        """
        average_abc = np.mean(abc)
        return np.sum(np.abs(abc - average_abc) / average_abc) / 3

    @classmethod
    def calc_supercell_isotropy(cls, structure, multi):
        abc = structure.lattice.abc
        super_abc = multi * abc
        return cls.calc_isotropy(super_abc)

    @staticmethod
    def supercell_comment(multi, isotropy_value):

        if multi.shape == (3, 3):
            multi_str = ' '.join([str(int(i)) for i in multi.flatten()])
        elif multi.shape == (3,):
            multi_str = ' '.join([str(int(i)) for i in multi])
        else:
            multi_str = str(multi)

        return 'multi: ' + multi_str + ', ' + 'isotropy: ' + \
               str(round(isotropy_value, 3)) + '\n'

    @property
    def multi(self):
        return self._multi

    @property
    def supercell_structure(self):
        return self._supercell_structure

    @property
    def num_atoms(self):
        return self._num_atoms

    @property
    def isotropy(self):
        return self._isotropy

    @property
    def comment(self):
        return self._comment


class Supercells:
    def __init__(self, structure, conventional=True, max_num_atoms=400,
                 min_num_atoms=50, isotropy_criterion=0.12):
        """
        Constructs a set of supercells satisfying a criterion.
        Args:
            structure (pmg structure class object):
                Unitcell structure
            conventional (bool):
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

        if conventional:
            conventional = find_spglib_standard_conventional(structure)
            primitive = find_spglib_standard_primitive(structure)[0]
            if conventional.num_sites == primitive.num_sites:
                logger.info("Primitive cell is same as conventional cell.")
                self._is_conventional_based = False
            else:
                self._is_conventional_based = True
            unitcell = conventional
        else:
            unitcell = find_spglib_standard_primitive(structure)[0]
            self._is_conventional_based = False

        self._sorted_unitcell = unitcell.get_sorted_structure()
        abc = np.array(self._sorted_unitcell.lattice.abc)
        num_atoms_in_unitcell = self._sorted_unitcell.num_sites

        if max_num_atoms < num_atoms_in_unitcell:
            raise TooLargeUnitcellError("Number of atoms in the unitcell is "
                                        "smaller than the maximum number of "
                                        "atoms in the supercell")

        multi = np.ones(3, dtype="int8")
        self._supercells = []

        for i in range(int(max_num_atoms / num_atoms_in_unitcell)):
            num_atoms = multi.prod() * num_atoms_in_unitcell
            if num_atoms > max_num_atoms:
                break

            isotropy = Supercell.calc_supercell_isotropy(self._sorted_unitcell,
                                                         multi)
            if isotropy < isotropy_criterion and num_atoms >= min_num_atoms:
                self._supercells.append(Supercell(self._sorted_unitcell, multi))

            super_abc = multi * abc
            # multi indices within 1.05a, where a is the shortest supercell
            # lattice length, are incremented.
            for j in range(3):
                if super_abc[j] / min(super_abc) < 1.05:
                    multi[j] += 1

        self._are_supercells = True if self._supercells else False

    def sorted_by_num_atoms(self):
        return sorted(deepcopy(self._supercells),
                      key=lambda x: (x.num_atoms, x.isotropy))

    def sorted_by_isotropy(self):
        return sorted(deepcopy(self._supercells),
                      key=lambda x: (x.isotropy, x.num_atoms))

    @property
    def are_supercells(self):
        return self._are_supercells

    @property
    def is_conventional_based(self):
        return self._is_conventional_based

    @property
    def unitcell_structure(self):
        return self._sorted_unitcell


class TooLargeUnitcellError(Exception):
    pass
