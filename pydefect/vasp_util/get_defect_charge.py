#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pymatgen.io.vasp.inputs import Poscar

from pydefect.vasp_util.util import get_defect_charge_from_vasp

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


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
