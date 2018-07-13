# -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy

from pymatgen.core.structure import Structure
from pydefect.util.structure import find_spglib_standard_conventional, \
    find_spglib_standard_primitive

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


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
    def __init__(self, structure, primitive=False, max_num_atoms=400,
                 min_num_atoms=50, isotropy_criterion=0.12):
        """
        Constructs a set of supercells satisfying a criterion.
        Args:
            structure (pmg structure class object):
                Unitcell structure
            primitive (bool):
                True: Only the conventional cell is expanded.
                False: Both conventional and primitive cells are expanded.
            max_num_atoms (int):
                The maximum number of atoms in the supercell.
            min_num_atoms (int):
                The minimum number of atoms in the supercell.
            isotropy_criterion (float):
                The criterion to judge if a supercell is isotropic or not.
        Return:
            supercell structure (Supercell):
                Supercell class object
            unitcell structure (Structure):
                pmg Structure class object for the based unitcell.
            multi (3x3 numpy array):
                Shape of the supercell
            calc_isotropy (float):
                Isotropic value.
        """
        if primitive is False:
            uc_structure = find_spglib_standard_conventional(structure)
            # check if the conventional cell is same as the primitive cell.
            primitive = find_spglib_standard_primitive(structure)
            if uc_structure.num_sites == primitive.num_sites:
                print("The conventional cell is same as the primitive cell.")
                self._is_primitive = True
            else:
                self._is_primitive = False
        else:
            uc_structure = find_spglib_standard_primitive(structure)
            self._is_primitive = True

        self._uc_structure = uc_structure.get_sorted_structure()
        multi = np.ones(3, dtype="int8")
        abc = np.array(self._uc_structure.lattice.abc)
        num_atoms_in_unitcell = self._uc_structure.num_sites

        if max_num_atoms < num_atoms_in_unitcell:
            raise SupercellSizeError("The number of atoms in the unitcell is "
                                     "smaller than the maximum number of "
                                     "atoms in the supercell")
        self._supercells = []

        for i in range(int(max_num_atoms / num_atoms_in_unitcell)):
            num_atoms = multi.prod() * num_atoms_in_unitcell
            if num_atoms > max_num_atoms:
                break

            # The supercell indices within 1.05a, where a is the shortest
            # supercell lattice length, are incremented.
            isotropy = Supercell.calc_supercell_isotropy(uc_structure, multi)
            if isotropy < isotropy_criterion and num_atoms >= min_num_atoms:
                self._supercells.append(Supercell(uc_structure, multi))

            super_abc = multi * abc
            for j in range(3):
                if super_abc[j] / min(super_abc) < 1.05:
                    multi[j] += 1

        if self._supercells:
            self._converged = True
        else:
            self._converged = False

    def sorted_by_atoms(self):
        return sorted(deepcopy(self._supercells),
                      key=lambda x: (x.num_atoms, x.isotropy))

    def sorted_by_isotropy(self):
        return sorted(deepcopy(self._supercells),
                      key=lambda x: (x.isotropy, x.num_atoms))

    @property
    def supercells(self):
        return self._supercells

    @property
    def converged(self):
        return self._converged

    @property
    def is_primitive(self):
        return self._is_primitive

    @property
    def unitcell_structure(self):
        return self._uc_structure


class SupercellSizeError(Exception):
    pass

