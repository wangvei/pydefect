# -*- coding: utf-8 -*-

import numpy as np

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
        Constructs supercell based on a multi 3x3 matrix.
        Args:
            structure (pmg structure class object):
            multi (3x3 numpy array ,list, or a scalar):
            comment (str):
        """
        if len(multi) == 9:
            multi = np.reshape(multi, (3, 3))
        elif len(multi) == 3:
            multi = np.array(multi)

        s = structure * multi
        self.structure = s.get_sorted_structure()
        super_abc = multi * np.array(structure.lattice.abc)
        if comment is None:
            self.comment = \
                supercell_comment(multi, supercell_isotropy(super_abc, multi))
        else:
            self.comment = comment

    @classmethod
    def from_poscar(cls, poscar, multi):
            structure = Structure.from_file(poscar)
            return cls(structure, multi)

    @classmethod
    def recommended_supercell(cls, structure, to_conventional, max_num_atoms,
                              min_num_atoms, isotropy_criterion, make_forcibly):
        """
        Constructs a recommended supercell.
        Note:
        The supercell indices for the axes of which lengths are within 1.05a,
        where a is the shortest supercell lattice length, are incremented.
        Args:
            structure (pmg structure class object):
            to_conventional (bool):
            max_num_atoms (int):
            min_num_atoms (int):
            isotropy_criterion (float):
            make_forcibly (bool):
        """

        if to_conventional:
            structure = find_spglib_standard_conventional(structure)
        else:
            structure = find_spglib_standard_primitive(structure)

        multi = np.ones(3)
        abc = np.array(structure.lattice.abc)
        num_atoms = structure.num_sites
        max_iter = int(max_num_atoms / num_atoms)

        lowest_iso = float(cls.supercell_isotropy(structure, multi))
        candidate_multi = np.copy(multi)

        for i in range(max_iter):
            super_abc = multi * abc
            reduced_abc = super_abc / min(super_abc)

            iso = cls.supercell_isotropy(structure, multi)

            if multi.prod() * num_atoms > max_num_atoms:
                break

            elif multi.prod() * num_atoms > min_num_atoms and \
                    iso < isotropy_criterion:
                comment = cls.supercell_comment(multi, iso, is_converged=True)
                return cls(structure, multi, comment)

            else:
                for j in range(3):
                    if reduced_abc[j] < 1.05:
                        multi[j] += 1
                if iso < lowest_iso:
                    lowest_iso = iso
                    candidate_multi = np.copy(multi)

        if make_forcibly is True:
            comment = cls.supercell_comment(candidate_multi, lowest_iso,
                                        is_converged=False)
            return cls(structure, candidate_multi, comment)

    def to_poscar(self, filename):
        poscar_str = self.structure.to(fmt="poscar").splitlines(True)
        poscar_str[0] = self.comment
        with open(filename, 'w') as fw:
            for line in poscar_str:
                fw.write(line)

    @staticmethod
    def isotropy(abc):
        average_abc = np.mean(abc)
        return np.sum(np.abs(abc - average_abc) / average_abc) / 3

    @classmethod
    def supercell_isotropy(cls, structure, multi):
        abc = structure.lattice.abc
        super_abc = multi * abc
        return cls.isotropy(super_abc)

    @staticmethod
    def supercell_comment(multi, isotropy_value, is_converged=None):

        if is_converged is True:
            converged = "Converged:"
        elif is_converged is False:
            converged = "Not Converged:"
        else:
            converged = ":"

        if multi.shape == (3, 3):
            multi_str = ' '.join([str(int(i)) for i in multi.flatten()])
        elif multi.shape == (3,):
            multi_str = ' '.join([str(int(i)) for i in multi])
        else:
            multi_str = str(multi)

        return 'multi: ' + multi_str + ', ' + \
               'supercell_isotropy ' + converged + str(isotropy_value) + '\n'


