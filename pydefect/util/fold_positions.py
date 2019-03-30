import sys
import math

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure


def fold_float(x):
    return x - math.floor(x)


def fold_positions(structure):
    """
    Fold atomic positions with fractional coords less than 0 or larger than 1
    into from 0 to 1.
    For example, coords of site changes from [-0.3, 1.9, 0.5] to [0.7, 0.9, 0.5]

    Args:
        structure(Structure):

    Returns:
        Structure
    """
    for i, site in enumerate(structure):
        modification_vector = [-math.floor(v) for v in site.frac_coords]
        structure.translate_sites(i, modification_vector)
    return structure


def fold_positions_in_poscar(poscar):
    """

    Args:
        poscar (Poscar):

    Returns:
        Poscar:

    """
    s = poscar.structure
    fold_positions(s)
    return Poscar(s)


if __name__ == "__main__":
    path = sys.argv[1]
    p = Poscar.from_file(path)
    fold_positions_in_poscar(p)
    print(p)

