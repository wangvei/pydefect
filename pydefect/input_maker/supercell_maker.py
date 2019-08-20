# -*- coding: utf-8 -*-

from collections import Iterable
from copy import deepcopy
from typing import Union

import numpy as np
from vise.util.structure_handler import (
    find_spglib_standard_conventional, find_spglib_standard_primitive)
from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.error_classes import CellSizeError
from pydefect.util.logger import get_logger
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def calc_isotropy(structure: Structure,
                  trans_mat: np.array) -> tuple:
    """ Return mean absolute deviation of lattice constants from their average

    Args:
        structure (Structure):
            Original structure to be expanded.
        trans_mat (3x3 np.array):
           The matrix to be used for expanding the unitcell.

    Return
        isotropy (float):
           Defined as the mean absolute deviation of the lattice constants
           from the averaged lattice constant
           np.sum(np.abs(abc - average_abc) / average_abc) / 3
        alpha (float):
            Angle alpha of lattice in degree.
    """
    new_structure = structure * trans_mat
    alpha = round(new_structure.lattice.alpha, 2)
    super_abc = new_structure.lattice.abc
    average_abc = np.mean(super_abc)

    isotropy = np.sum(np.abs(super_abc - average_abc) / average_abc) / 3
    return round(isotropy, 4), alpha


class Supercell:
    def __init__(self,
                 structure: Structure,
                 trans_mat: Union[int, np.array],
                 multiplicity: int):
        """ Supercell class constructed based on a given multiplicity.

        Args:
            structure (Structure):
                Original structure to be expanded.
            trans_mat (3x3 np.array, 3 np.array or a scalar):
                The matrix to be used for expanding the unitcell.
            multiplicity (int):
                The size multiplicity of structure wrt the primitive cell.
        """
        if isinstance(trans_mat, Iterable):
            trans_mat = list(trans_mat)
        else:
            trans_mat = [trans_mat]

        if len(trans_mat) == 1:
            trans_mat_str = str(trans_mat)
        elif len(trans_mat) == 9:
            trans_mat_str = ' '.join([str(i) for i in trans_mat])
            trans_mat = np.reshape(trans_mat, (3, 3))
        elif len(trans_mat) == 3:
            if isinstance(trans_mat[0], Iterable) and len(trans_mat[0]) == 3:
                trans_mat = np.array(trans_mat)
                trans_mat_str = \
                    ' '.join([str(int(i)) for i in trans_mat.flatten()])
            else:
                trans_mat = np.array(trans_mat)
                trans_mat_str = ' '.join([str(int(i)) for i in trans_mat])
        else:
            raise ValueError(f"Translation matrix: {trans_mat} is not proper. "
                             f"1, 3, or 9 components are accepted.")

        self.base_structure = structure
        s = structure * trans_mat
        self.structure = s.get_sorted_structure()
        self.trans_mat = trans_mat
        self.multiplicity = multiplicity
        self.isotropy = calc_isotropy(structure, trans_mat)
        self.num_atoms = self.structure.num_sites

        self.comment = f"trans_mat: {trans_mat_str}, multi: {multiplicity}, " \
                       f"isotropy: {self.isotropy[0]}\n"

    def to_poscar(self, poscar_filename: str) -> None:
        poscar_str = self.structure.to(fmt="poscar").splitlines(True)
        poscar_str[0] = self.comment

        with open(poscar_filename, 'w') as fw:
            for line in poscar_str:
                fw.write(line)

    def to_uposcar(self, uposcar_filename: str) -> None:
        self.base_structure.to(filename=uposcar_filename)


class Supercells:
    def __init__(self,
                 structure: Structure,
                 conventional_base: bool = True,
                 max_num_atoms: int = 400,
                 min_num_atoms: int = 50,
                 criterion: float = 0.12,
                 rhombohedral_angle: float = 70,
                 symprec: float = SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL):
        """ Constructs a set of supercells satisfying an isotropic criterion.

        Args:
            structure (pmg structure class object):
                Unitcell structure
            conventional_base (bool):
                Conventional cell is expanded when True, otherwise primitive
                cell is expanded.
            max_num_atoms (int):
                Maximum number of atoms in the supercell.
            min_num_atoms (int):
                Minimum number of atoms in the supercell.
            criterion (float):
                Criterion to judge if a supercell is isotropic or not.
            rhombohedral_angle (float):
                Rhombohedral primitive cells may have very small or very large
                lattice angles not suited for first-principles calculations. If
                rhombohedral_angle is set, only the supercells with
                rhombohedral_angle <= alpha <= 180 - rhombohedral_angle are
                returned. Then, the new supercells are iteratively
                created by multiplying [[1, 1, -1], [-1, 1, 1], [1, -1, 1]] or
                [[1, 1, 0], [0, 1, 1], [1, 0, 1]].
            symprec (float):
                Precision used for symmetry analysis in angstrom.
            angle_tolerance (float):
                Angle tolerance for symmetry analysis in degree
        """
        primitive_cell, _ = \
            find_spglib_standard_primitive(structure, symprec, angle_tolerance)
        sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
        symmetry_dataset = sga.get_symmetry_dataset()
        logger.info(f"Space group: {symmetry_dataset['international']}")

        if conventional_base:
            unitcell = find_spglib_standard_conventional(structure,
                                                         symprec,
                                                         angle_tolerance)
            rhombohedral = False
            self.conventional_base = True
        else:
            unitcell = primitive_cell
            rhombohedral = sga.get_lattice_type() == "rhombohedral"
            self.conventional_base = False

        self.unitcell = unitcell.get_sorted_structure()
        unitcell_mul = int(len(self.unitcell) / len(primitive_cell))

        if max_num_atoms < len(primitive_cell):
            raise CellSizeError("Number of atoms in unitcell is too large.")

        self.supercells = []
        # Only rhombohedral primitive is not an identical matrix.
        based_trans_mat = np.identity(3, dtype="int8")
        # Normal incremented matrix one by one
        incremented_mat = np.identity(3, dtype="int8")
        trans_mat = np.identity(3, dtype="int8")

        for i in range(int(max_num_atoms / len(primitive_cell))):
            isotropy, angle = calc_isotropy(self.unitcell, trans_mat)
            multiplicity = int(unitcell_mul * round(np.linalg.det(trans_mat)))
            num_atoms = int(multiplicity * len(primitive_cell))
            if num_atoms > max_num_atoms:
                break

            rhombohedral_shape = None
            if rhombohedral_angle and rhombohedral:
                if angle < rhombohedral_angle:
                    rhombohedral_shape = "sharp"
                elif angle > 180 - rhombohedral_angle:
                    rhombohedral_shape = "blunt"

            if isotropy < criterion and num_atoms >= min_num_atoms:
                self.supercells.append(
                    Supercell(self.unitcell, trans_mat, multiplicity))

            if rhombohedral_shape == "sharp":
                multiplied_matrix = np.array([[ 1,  1, -1],
                                              [-1,  1,  1],
                                              [ 1, -1,  1]])
                based_trans_mat = np.dot(multiplied_matrix, based_trans_mat)
            elif rhombohedral_shape == "blunt":
                multiplied_matrix = np.array([[1, 1, 0],
                                              [0, 1, 1],
                                              [1, 0, 1]])
                based_trans_mat = np.dot(multiplied_matrix, based_trans_mat)
            else:
                super_abc = (self.unitcell * trans_mat).lattice.abc
                # multi indices within 1.05a, where a is the shortest supercell
                # lattice length, are incremented at the same time.
                for j in range(3):
                    if super_abc[j] / min(super_abc) < 1.05:
                        incremented_mat[j, j] += 1

            trans_mat = np.dot(incremented_mat, based_trans_mat)

    @property
    def sorted_supercells_by_num_atoms(self) -> list:
        return sorted(deepcopy(self.supercells),
                      key=lambda x: (x.num_atoms, x.isotropy))

    @property
    def sorted_supercells_by_isotropy(self) -> list:
        return sorted(deepcopy(self.supercells),
                      key=lambda x: (x.isotropy, x.num_atoms))

    @property
    def smallest_supercell(self) -> Supercell:
        return self.sorted_supercells_by_num_atoms[0]

    @property
    def most_isotropic_supercell(self) -> Supercell:
        return self.sorted_supercells_by_isotropy[0]

