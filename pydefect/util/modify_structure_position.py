import sys
import math

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure


def modify_float(x):
    return x - math.floor(x)


def modify_position_structure(structure):
    """
    modify positions which is out of box (fractional coords x<0 or x>=1)
    into box (0 <= x < 1)
    For example, coords of site changes from [-0.3, 1.9, 0.5] to [0.7, 0.9, 0.5]

    Args:
        structure(Structure):

    Returns:
        None

    """
    for i, site in enumerate(structure):
        modification_vector = [-math.floor(v) for v in site.frac_coords]
        structure.translate_sites(i, modification_vector)
    return structure


def modify_position_poscar(poscar):
    """

    Args:
        poscar(Poscar):

    Returns:

    """
    s = poscar.structure
    modify_position_structure(s)
    return Poscar(s)


if __name__ == "__main__":
    path = sys.argv[1]
    p = Poscar.from_file(path)
    modify_position_poscar(p)
    print(p)

