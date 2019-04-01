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
        Constructs a supercell based on a given multiplicity.
        Args:
            structure (Structure):
                Original structure to be expanded.
            multi (3x3 list, list or a scalar):
                The matrix to be used for expanding the unitcell.
            comment (str):
                Any comment.
        """
        print(multi)
        if len(multi) == 1:
            multi_str = str(multi)
        elif len(multi) == 9:
            multi = np.reshape(multi, (3, 3))
            multi_str = ' '.join([str(int(i)) for i in multi])
        elif len(multi) == 3:
            multi = np.array(multi)
            if isinstance(multi[0], list) and len(multi[0]) == 3:
                multi_str = ' '.join([str(int(i)) for i in multi.flatten()])
            else:
                multi_str = ' '.join([str(int(i)) for i in multi])
        else:
             ValueError("Multiplicity: {} is not proper. 1, 3, or 9 numbers"
                       " are accepted.")

        s = structure * multi
        self.multi = multi
        print(s)
        self.structure = s.get_sorted_structure()
        self.isotropy = calc_isotropy(self.structure, multi)
        self.num_atoms = self.structure.num_sites

        if comment is None:
            self.comment = 'multi: ' + multi_str + ', ' + 'isotropy: ' + \
                           str(round(self.isotropy, 4)) + '\n'
        else:
            self.comment = comment

    def to_poscar(self, filename):
        poscar_str = self.structure.to(fmt="poscar").splitlines(True)
        poscar_str[0] = self.comment
        with open(filename, 'w') as fw:
            for line in poscar_str:
                fw.write(line)


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
    return np.sum(np.abs(super_abc - average_abc) / average_abc) / 3


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
                self.is_conventional_based = False
            else:
                self.is_conventional_based = True
            unitcell = conventional
        else:
            unitcell = find_spglib_standard_primitive(structure)[0]
            self.is_conventional_based = False

        self._sorted_unitcell = unitcell.get_sorted_structure()
        abc = np.array(self._sorted_unitcell.lattice.abc)
        num_atoms_in_unitcell = self._sorted_unitcell.num_sites

        if max_num_atoms < num_atoms_in_unitcell:
            raise TooLargeUnitcellError("Number of atoms in the unitcell is "
                                        "smaller than the maximum number of "
                                        "atoms in the supercell")

        multi = np.ones(3, dtype="int8")
        self.supercells = []

        for i in range(int(max_num_atoms / num_atoms_in_unitcell)):
            num_atoms = multi.prod() * num_atoms_in_unitcell
            if num_atoms > max_num_atoms:
                break

            isotropy = calc_isotropy(self._sorted_unitcell, multi)
            m = deepcopy(multi).tolist()
            if isotropy < isotropy_criterion and num_atoms >= min_num_atoms:
                self.supercells.append(Supercell(self._sorted_unitcell, m))

            super_abc = multi * abc
            # multi indices within 1.05a, where a is the shortest supercell
            # lattice length, are incremented.
            for j in range(3):
                if super_abc[j] / min(super_abc) < 1.05:
                    multi[j] += 1

        self._are_supercells = True if self.supercells else False

    def sorted_by_num_atoms(self):
        return sorted(deepcopy(self.supercells),
                      key=lambda x: (x.num_atoms, round(x.isotropy, 4)))

    def sorted_by_isotropy(self):
        return sorted(deepcopy(self.supercells),
                      key=lambda x: (round(x.isotropy, 4), x.num_atoms))

    @property
    def unitcell_structure(self):
        return self._sorted_unitcell


class TooLargeUnitcellError(Exception):
    pass
