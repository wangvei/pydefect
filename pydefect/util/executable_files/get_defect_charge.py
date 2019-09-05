#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pymatgen.io.vasp.inputs import Poscar


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

def get_defect_charge_from_vasp(structure, potcar="POTCAR", incar="INCAR"):
    """
    Returns the defect charge by comparing nion, number of electrons in POTCAR,
    and NELECT in INCAR.
    """
    try:
        num_elect_incar = Incar.from_file(incar)["NELECT"]
    # return 0 if NELECT is absent in INCAR
    except KeyError:
        return 0

    potcar = Potcar.from_file(potcar)

    num_elect_neutral = 0

    for i, e in enumerate(structure.composition):
        if potcar[i].element != str(e):
            raise ValueError("The sequence of elements in POTCAR and Structure "
                             "is different.")

        nelect = potcar[i].nelectrons
        nions = structure.composition[e]
        num_elect_neutral += nelect * nions

    # charge is minus of difference of the electrons
    return int(num_elect_neutral - num_elect_incar)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str,
                        help="POSCAR-type file name for the unitcell.",
                        metavar="FILE")
    parser.add_argument("--potcar", dest="potcar", default="POTCAR", type=str,
                        help="POTCAR-type file name for the unitcell.",
                        metavar="FILE")
    parser.add_argument("--incar", dest="incar", default="INCAR", type=str,
                        help="INCAR-type file name.", metavar="FILE")

    opts = parser.parse_args()

    pos = Poscar.from_file(opts.poscar)
    nions = pos.natoms
    charge = get_defect_charge_from_vasp(nions, opts.potcar, opts.incar)

    print("composition:", pos.structure.composition)
    print("defect charge:", charge)


if __name__ == "__main__":
    main()
