#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pymatgen.io.vasp.inputs import Poscar

from pydefect.vasp_util.util import get_num_electrons_from_potcar

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str,
                        help="POSCAR-type file name for the unitcell.",
                        metavar="FILE")
    parser.add_argument("--potcar", dest="potcar", default="POTCAR", type=str,
                        help="POTCAR-type file name for the unitcell.",
                        metavar="FILE")
    parser.add_argument("-c", "--charge", dest="charge", default="0",
                        type=float, help="Charge.")

    opts = parser.parse_args()

    pos = Poscar.from_file(opts.poscar)
    nions = pos.natoms
    nelect = get_num_electrons_from_potcar(opts.potcar, nions, opts.charge)

    print("NELECT =", nelect)


if __name__ == "__main__":
    main()
