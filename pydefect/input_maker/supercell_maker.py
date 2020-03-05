# -*- coding: utf-8 -*-

from collections import Iterable
from copy import deepcopy
from typing import Union, Optional, List

import numpy as np
from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.error_classes import CellSizeError
from pydefect.util.logger import get_logger
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pydefect.database.symmetry import transmat_primitive2standard

from vise.util.structure_handler import (
    find_spglib_primitive, get_symmetry_dataset)

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def calc_isotropy(structure: Structure,
                  trans_mat: np.array) -> tuple:
    """Return mean absolute deviation of lattice constants from their average.

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
            Average of lattice angle in degree.
    """
    new_structure = structure * trans_mat
    angle = sum(new_structure.lattice.angles) / 3
    super_abc = new_structure.lattice.abc
    average_abc = np.mean(super_abc)

    isotropy = np.sum(np.abs(super_abc - average_abc) / average_abc) / 3
    return round(isotropy, 4), angle


def sanitize_matrix(matrix: list) -> np.ndarray:
    """Sanitize the matrix component to 3x3 np array with int type.

    Example1:  [1, 0, 0, 0, 1, 0, 0, 0, 1] -> [[1, 0, 0],
                                               [0, 1, 0],
                                               [0, 0, 1]]
    Example2:  [1, 2, 3] -> [[1, 0, 0],
                             [0, 2, 0],
                             [0, 0, 3]]

    Args:
        matrix (list or np.array):
           The matrix to be used for expanding the structure.

    Return:
        np.ndarray (3x3)
    """
    if len(matrix) == 1:
        sanitized_matrix = np.eye(3, dtype=int)
        for i in range(3):
            sanitized_matrix[i, i] = matrix[0]
    elif len(matrix) == 9:
        sanitized_matrix = np.reshape(matrix, (3, 3))
    elif len(matrix) == 3:
        if isinstance(matrix[0], Iterable) and len(matrix[0]) == 3:
            sanitized_matrix = np.array(matrix, dtype=int)
        else:
            sanitized_matrix = np.eye(3, dtype=int)
            for i in range(3):
                sanitized_matrix[i, i] = matrix[i]
    else:
        raise ValueError(f"Transformation matrix {matrix} is not proper."
                         f"Only 1, 3, or 9 components are accepted.")

    return sanitized_matrix


class Supercell:
    def __init__(self,
                 structure: Structure,
                 trans_mat: Union[int, np.array, List[List]],
                 multiplicity: Optional[int] = None,
                 check_unitcell: bool = False,
                 symprec: float = SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL):
        """Supercell class constructed based on a given multiplicity.

        Args:
            structure (Structure):
                Primitive ell structure to be expanded.
            trans_mat (3x3 np.array, 3 np.array or a scalar):
                The matrix to be used for expanding the structure.
            multiplicity (int):
                The size multiplicity of structure wrt the primitive cell.
        """
        trans_mat = sanitize_matrix(trans_mat)

        self.is_structure_changed = False
        if check_unitcell:
            primitive, self.is_structure_changed = \
                find_spglib_primitive(structure, symprec, angle_tolerance)
            if self.is_structure_changed:
                logger.warning(f"Structure is changed to primitive cell.")
                sym_dataset = get_symmetry_dataset(structure, symprec,
                                                   angle_tolerance)
                trans_mat = sym_dataset["transformation_matrix"] * trans_mat
                origin_shift = sym_dataset["origin_shift"]
                if np.linalg.norm(origin_shift) > 1e-5:
                    logger.warning(f"Origin is shifted by {origin_shift}.")

            structure = primitive

        s = structure * trans_mat
        self.structure = s.get_sorted_structure()
        self.trans_mat = trans_mat
        self.multiplicity = multiplicity or int(round(np.linalg.det(trans_mat)))
        self.isotropy = calc_isotropy(structure, trans_mat)
        self.num_atoms = self.structure.num_sites
        flatten_trans_mat = sum(self.trans_mat.tolist(), [])
        self.comment = ", ".join(["trans_mat: " + str(flatten_trans_mat),
                                  "multi: " + str(self.multiplicity),
                                  "isotropy: " + str(self.isotropy)])

    def to(self, poscar: str, uposcar: Optional[str] = "UPOSCAR") -> None:
        """Write out the Supercell POSCAR and unitcell UPOSCAR.

        Args:
            poscar (str): POSCAR file name.
            uposcar (str): UPOSCAR file name.

        Return: None
        """
        self.structure.to(filename=poscar, comment=self.comment)
        if self.is_structure_changed:
            self.structure.to(filename=uposcar)


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
        """Constructs a set of supercells satisfying an isotropic criterion.

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
                lattice angles not suited for first-principles calculations.
                Therefore, only the supercells with
                rhombohedral_angle <= lattice angle <= 180 - rhombohedral_angle
                are returned. Then, the new supercells are iteratively
                created by multiplying [[1, 1, -1], [-1, 1, 1], [1, -1, 1]] or
                [[1, 1, 0], [0, 1, 1], [1, 0, 1]].
            symprec (float):
                Precision used for symmetry analysis in angstrom.
            angle_tolerance (float):
                Angle tolerance for symmetry analysis in degree.
        """
        primitive_cell, _ = \
            find_spglib_primitive(structure, symprec, angle_tolerance)
        if max_num_atoms < len(primitive_cell):
            raise CellSizeError("Number of atoms in unitcell is too large.")

        self.unitcell = primitive_cell.get_sorted_structure()
        sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
        symmetry_dataset = sga.get_symmetry_dataset()
        logger.info(f"Space group: {symmetry_dataset['international']}")

        self.conventional_base = conventional_base
        lattice = sga.get_lattice_type()
        centering = symmetry_dataset["international"][0]
        if conventional_base:
            to_unitcell_mat = transmat_primitive2standard(centering)
            rhombohedral = False
            tetragonal = lattice == "tetragonal"
        else:
            to_unitcell_mat = np.eye(3, dtype=int)
            rhombohedral = lattice == "rhombohedral"
            if lattice == "tetragonal" and centering == "P":
                tetragonal = True
            else:
                tetragonal = False

        self.supercells = []
        # Isotropically incremented matrix one by one
        trans_mat = to_unitcell_mat
        incremented_mat = np.eye(3, dtype=int)
        rotation_matrix = np.eye(3, dtype=int)

        for i in range(int(max_num_atoms / len(primitive_cell))):
            isotropy, angle = calc_isotropy(self.unitcell, trans_mat)
            # int is needed when numpy.float is rounded.
            multiplicity = int(round(np.linalg.det(trans_mat)))
            num_atoms = multiplicity * len(primitive_cell)
            if num_atoms > max_num_atoms:
                break

            if isotropy < criterion and num_atoms >= min_num_atoms:
                self.supercells.append(
                    Supercell(self.unitcell, trans_mat, multiplicity))

            rhombohedral_shape = None
            if rhombohedral:
                if angle < rhombohedral_angle:
                    rhombohedral_shape = "sharp"
                elif angle > 180 - rhombohedral_angle:
                    rhombohedral_shape = "blunt"

            if rhombohedral_shape == "sharp":
                m = np.array([[1,  1, -1], [-1,  1,  1], [1, -1,  1]])
                rotation_matrix = np.dot(rotation_matrix, m)
            elif rhombohedral_shape == "blunt":
                m = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
                rotation_matrix = np.dot(rotation_matrix, m)
            else:
                if tetragonal:
                    s = self.unitcell * trans_mat
                    if s.lattice.a < s.lattice.c:
                        if rotation_matrix[0, 1] == 0:
                            a = incremented_mat[0, 0]
                            if (a + 1) ** 2 < a ** 2 * 2:
                                incremented_mat[0, 0] += 1
                                incremented_mat[1, 1] += 1
                            else:
                                rotation_matrix = \
                                    np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]])
                        else:
                            incremented_mat[0, 0] += 1
                            incremented_mat[1, 1] += 1
                            rotation_matrix = np.eye(3, dtype=int)

                        incremented_mat *= rotation_matrix
                    else:
                        incremented_mat[2, 2] += 1
                else:
                    super_abc = (self.unitcell * trans_mat).lattice.abc
                    # multi indices within 1.05a, where a is the shortest
                    # supercell lattice length, are incremented simultaneously.
                    for j in range(3):
                        if super_abc[j] / min(super_abc) < 1.05:
                            incremented_mat[j, j] += 1

            trans_mat = np.dot(np.dot(incremented_mat, rotation_matrix),
                               to_unitcell_mat)

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

    def to_uposcar(self, uposcar: Optional[str] = "UPOSCAR") -> None:
        """ Write out the unitcell UPOSCAR file.

        Args:
            uposcar (str): UPOSCAR file name.

        Return: None
        """
        self.unitcell.to(filename=uposcar)

