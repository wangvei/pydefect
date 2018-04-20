#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import warnings

from pydefect.input_maker.defect_initial_setting \
    import print_dopant_info, DefectInitialSetting
from pydefect.input_maker.vasp_input_maker \
    import make_incar, make_kpoints, VaspDefectInputSetMaker
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults


__version__ = "0.0.1"
__date__ = "19.4.2018"

# Following defaults determine the condition of automatic defect calculations.
# electronegativity difference for antisites and substitutional impurities
_EN_DIFF = 1.0
# Maximum displacement distance
_DISTANCE = 0.2
# Cutoff radius in which atoms are perturbed.
_CUTOFF = 3.0
_SYMPREC = 0.01


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package to perform various things related to 
    first-principles point defect calculations. It allows us to construct input
    files, parse first-principles calculation results, and analyze data.""",
        epilog="""                                 
    Author: Yu Kumagai, Akira Takahashi
    Version: {}                                                                 
    Last updated: {}""".format(__version__, __date__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#        allow_abbrev=False)

    subparsers = parser.add_subparsers()

    # -- initial_setting -------------------------------------------------------
    parser_initial = subparsers.add_parser(
        name="initial_setting",
        description="Tools for configuring initial settings for a set of "
                    "defect calculations.",
        aliases=['initial', 'in'])

    parser_initial.add_argument(
        "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR-type file name for the unitcell.")
    parser_initial.add_argument(
        "-d", "--dopants", dest="dopants", default="", nargs="+", type=str,
        help="Dopant elements, e.g., Ga In.")
    parser_initial.add_argument(
        "-i", dest="interstitial_coords", nargs="+", default=None, type=float,
        help="Interstitial coordinates. Eg., 0.5 0.5 0.5.")
    parser_initial.add_argument(
        "-a", "--antisite", dest="is_antisite", action="store_false",
        help="Set if antisite defects are not considered.")
    parser_initial.add_argument(
        "-e", dest="en_diff", type=float, default=_EN_DIFF,
        help="Criterion of the electronegativity difference that determines "
             "antisites and/or substituted impurities.")
    parser_initial.add_argument(
        "--included", dest="included", type=str, default="", nargs="+",
        help="Exceptionally included defects. E.g., Va_O2_-1.")
    parser_initial.add_argument(
        "--excluded", dest="excluded", type=str, default="", nargs="+",
        help="Exceptionally excluded defects. E.g., Va_O2_0.")
    parser_initial.add_argument(
        "--distance", dest="distance", type=float, default=_DISTANCE,
        help="Displacement distance. 0 means that random displacement is not "
             "considered.")
    parser_initial.add_argument(
        "--cutoff", dest="cutoff", type=float, default=_CUTOFF,
        help="Set the cutoff radius [A] in which atoms are displaced.")
    parser_initial.add_argument(
        "--symprec", dest="symprec", type=float, default=_SYMPREC,
        help="Set precision used for symmetry analysis [A].")
    parser_initial.add_argument(
        "--print_dopant", dest="print_dopant", default=None, type=str,
        help="Print dopant information that can be added a posteriori.")

    #    groups = parser_input.add_mutually_exclusive_group(required=True)
    parser_initial.set_defaults(func=initial_setting)

    # -- vasp_input_maker ------------------------------------------------------
    parser_vasp_input = subparsers.add_parser(
        name="vasp_input",
        description="Tools for configuring vasp input files for a set of "
                    "defect calculations. One needs to set .pydefect.yaml"
                    "for potcar setup.",
        aliases=['vasp', 'vi'])

    parser_vasp_input.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_vasp_input.add_argument(
        "--dposcar", dest="dposcar", default="DPOSCAR", type=str,
        help="DPOSCAR-type file name.")
    parser_vasp_input.add_argument(
        "--incar", dest="incar", default="INCAR", type=str,
        help="INCAR-type file name.")
    parser_vasp_input.add_argument(
        "--kpoints", dest="kpoints", default="KPOINTS", type=str,
        help="KPOINTS-type file name.")
    parser_vasp_input.add_argument(
        "--filtering", dest="filtering", type=str, default=None, nargs="+",
        help="Filtering kwargs.")
    parser_vasp_input.add_argument(
        "--add", dest="add", type=str, default=None, nargs="+",
        help="Particular defect names to be added.")
    parser_vasp_input.add_argument(
        "--force_overwrite", dest="force_overwrite", default=False,
        help="Set if the existing folders are overwritten.")
    parser_vasp_input.add_argument(
        "--make_incar", dest="make_incar", action="store_true",
        help="Make INCAR file using several default setting.")
    parser_vasp_input.add_argument(
        "--make_kpoints", dest="make_kpoints", action="store_true",
        help="Make KPOINTS file based on the lattice constants.")

    parser_vasp_input.set_defaults(func=vasp_input)

    # -- defect_entry ----------------------------------------------------------
    parser_defect_entry = subparsers.add_parser(
        name="defect_entry",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
        aliases=['entry', 'de'])

    parser_defect_entry.add_argument(
        "--make_defect_entry", dest="make_defect_entry", action="store_true",
        help="Make defect_entry.json from yaml file.")
    parser_defect_entry.add_argument(
        "--yaml", dest="yaml", type=str,
        help="defect_entry.yaml-type file name to be read.")
    parser_defect_entry.add_argument(
        "--json", dest="json", type=str, default="defect_entry.json",
        help="defect_entry.json type file name.")
    parser_defect_entry.add_argument(
        "--print", dest="print", action="store_true",
        help="Print DefectEntry class object information.")

    parser_defect_entry.set_defaults(func=defect_entry)

    # -- supercell_dft_results -------------------------------------------------
    parser_supercell_results = subparsers.add_parser(
        name="supercell_results",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
        aliases=['supercell', 'sr'])

    parser_supercell_results.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Directory names.")
    parser_supercell_results.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and " "perfect directory.")
    parser_supercell_results.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_supercell_results.add_argument(
        "-o", dest="outcar", type=str, default="OUTCAR")
    parser_supercell_results.add_argument(
        "-v", dest="vasprun", type=str, default="vasprun.xml")

    parser_supercell_results.set_defaults(func=supercell_results)

    # -- unitcell_dft_results -------------------------------------------------
    parser_unitcell_results = subparsers.add_parser(
        name="unitcell_results",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
        aliases=['unitcell', 'ur'])

    parser_unitcell_results.add_argument(
        "--json_file", dest="json_file", default=None, type=str)
    parser_unitcell_results.add_argument(
        "--band_edge_dir", dest="band_edge_dir", default=None, type=str)
    parser_unitcell_results.add_argument(
        "--static_diele_dir", dest="static_diele_dir", default=None, type=str)
    parser_unitcell_results.add_argument(
        "--ionic_diele_dir", dest="ionic_diele_dir", default=None, type=str)
    parser_unitcell_results.add_argument(
        "--total_dos_dir", dest="total_dos_dir", default=None, type=str)

    parser_unitcell_results.set_defaults(func=unitcell_results)

    args = parser.parse_args()
    args.func(args)


def initial_setting(args):

    if args.print_dopant:
        print_dopant_info(args.print_dopant)
    else:
        defect_setting = DefectInitialSetting.from_basic_settings(
            args.poscar, args.dopants, args.interstitial_coords,
            args.is_antisite, args.en_diff, args.included, args.excluded,
            args.distance, args.cutoff, args.symprec)
        defect_setting.to()


def vasp_input(args):

    if args.make_incar:
        make_incar(defect_in=args.defect_in)
    elif args.make_kpoints:
        make_kpoints(poscar=args.dposcar)
    else:
        defect_initial_setting = DefectInitialSetting.\
            from_defect_in(poscar=args.dposcar, defect_in_file=args.defect_in)

        VaspDefectInputSetMaker(defect_initial_setting=defect_initial_setting,
                                filtering_words=args.filtering,
                                particular_defects=args.add,
                                incar=args.incar,
                                kpoints=args.kpoints,
                                force_overwrite=args.force_overwrite)


def defect_entry(args):

    if args.make_defect_entry:
        defect_entry_from_yaml = DefectEntry.from_yaml(args.yaml)
        defect_entry_from_yaml.to_json_file("defect_entry.json")
    elif args.print:
        print(DefectEntry.json_load(args.json))


def supercell_results(args):
    if args.dir_all:
        from glob import glob
        dirs = glob('*[0-9]/')
        dirs.append("perfect")
    else:
        dirs = args.dirs

    for d in dirs:
        print(d)
        if os.path.isdir(d):
            try:
                dft_results = SupercellDftResults. \
                    from_vasp_files(d,
                                    contcar_name=args.poscar,
                                    outcar_name=args.outcar,
                                    vasprun_name=args.vasprun)
                dft_results.to_json_file(
                    filename=os.path.join(d, "dft_results.json"))
            except:
                warnings.warn(message="Parsing data in " + d + " is failed.")
        else:
            warnings.warn(message=d + " does not exist, so nothing is done.")


def unitcell_results(args):

    # if opts.json_file:
    #     try:
    #         dft_results = UnitcellDftResults.json_load(filename=opts.json_file)
    #     except IOError:
    #         print(opts.json_file, "does not exist.")
    # else:
    #     dft_results = UnitcellDftResults()

    dft_results = UnitcellDftResults()

    if args.band_edge_dir:
        try:
            dft_results.set_band_edge_from_vasp(args.band_edge_dir)
        except IOError:
            print(args.band_edge_dir, "is not appropriate.")

    if args.static_diele_dir:
        try:
            dft_results.\
                set_static_dielectric_tensor_from_vasp(args.static_diele_dir)
        except IOError:
            print(args.static_diele_dir, "is not appropriate.")

    if args.ionic_diele_dir:
        try:
            dft_results.\
                set_ionic_dielectric_tensor_from_vasp(args.ionic_diele_dir)
        except IOError:
            print(args.ionic_diele_dir, "is not appropriate.")

    if args.total_dos_dir:
        try:
            dft_results.set_total_dos_from_vasp(args.total_dos_dir)
        except IOError:
            print(args.total_dos_dir, "is not appropriate.")

    dft_results.to_json_file(args.json_file)


if __name__ == "__main__":
    main()
