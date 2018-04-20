#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import warnings

from pydefect.input_maker.defect_initial_setting \
    import print_dopant_info, DefectInitialSetting
from pydefect.input_maker.vasp_input_maker \
    import make_incar, make_kpoints, VaspDefectInputSetMaker
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.core.correction import Ewald, Correction


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
        "--json_file", dest="json_file", default="unitcell.json", type=str)
    parser_unitcell_results.add_argument(
        "--band_edge_dir", dest="band_edge_dir", default=None, type=str,
        help="Set band edge from a vasprun.xml file")

    parser_unitcell_results.add_argument(
        "--static_diele", dest="static_diele", default=None, type=float,
        nargs="+",
        help="Set static dielectric constant")
    parser_unitcell_results.add_argument(
        "--ionic_diele", dest="ionic_diele", default=None, type=float,
        nargs="+",
        help="Set ionic dielectric constant")

    parser_unitcell_results.add_argument(
        "--static_diele_dir", dest="static_diele_dir", default=None, type=str,
        help="Set static dielectric constant from an OUTCAR file")
    parser_unitcell_results.add_argument(
        "--ionic_diele_dir", dest="ionic_diele_dir", default=None, type=str,
        help="Set ionic dielectric constant from an OUTCAR file")
    parser_unitcell_results.add_argument(
        "--total_dos_dir", dest="total_dos_dir", default=None, type=str,
        help="Set total density of states from a vasprun.xml file")
    parser_unitcell_results.add_argument(
        "-o", dest="outcar", type=str, default="OUTCAR")
    parser_unitcell_results.add_argument(
        "-v", dest="vasprun", type=str, default="vasprun.xml")
    parser_unitcell_results.add_argument(
        "--print", dest="print", action="store_true",
        help="Print Unitcell class object information.")

    parser_unitcell_results.set_defaults(func=unitcell_results)

    # -- correction -------------------------------------------------
    parser_correction = subparsers.add_parser(
        name="correction",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
        aliases=['c'])

    # needed files
    parser_correction.add_argument(
        "--unitcell_json", dest="unitcell_json",
        default="unitcell.json", type=str)
    parser_correction.add_argument(
        "--perfect_json", dest="perfect_json",
        default="perfect/dft_results.json", type=str)

    # ewald
    parser_correction.add_argument(
        "--read_ewald_json", dest="read_ewald_json",
        default="ewald.json", type=str)
    parser_correction.add_argument(
        "--dump_ewald_json", dest="dump_ewald_json",
        default="ewald.json", type=str)
    parser_correction.add_argument(
        "--ewald_init", dest="ewald_init",
        default=None)
    parser_correction.add_argument(
        "--ewald_convergence", dest="ewald_convergence",
        default=None)
    parser_correction.add_argument(
        "--ewald_accuracy", dest="ewald_accuracy",
        default=None)

    # correction
    parser_correction.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Directory names.")
    parser_correction.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and " "perfect directory.")

    parser_correction.set_defaults(func=correction)

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

    try:
        dft_results = UnitcellDftResults.json_load(filename=args.json_file)
    except IOError:
        print(args.json_file, "does not exist.")

    if args.print:
        print(dft_results)
        return None

    if args.band_edge_dir:
        try:
            dft_results.set_band_edge_from_vasp(args.band_edge_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            print(args.band_edge_dir, "is not appropriate.")

    if args.static_diele:
        dft_results.static_dielectric_tensor = args.static_diele
    elif args.static_diele_dir:
        try:
            dft_results.\
                set_static_dielectric_tensor_from_vasp(args.static_diele_dir,
                                                       outcar_name=args.outcar)
        except IOError:
            print(args.static_diele_dir, "is not appropriate.")

    if args.ionic_diele:
        dft_results.ionic_dielectric_tensor = args.ionic_diele
    elif args.ionic_diele_dir:
        try:
            dft_results.\
                set_ionic_dielectric_tensor_from_vasp(args.ionic_diele_dir,
                                                      outcar_name=args.outcar)
        except IOError:
            print(args.ionic_diele_dir, "is not appropriate.")

    if args.total_dos_dir:
        try:
            dft_results.set_total_dos_from_vasp(args.total_dos_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            print(args.total_dos_dir, "is not appropriate.")

    dft_results.to_json_file(args.json_file)


def correction(args):
    from glob import glob

    try:
        unitcell_dft_data = UnitcellDftResults.json_load(args.unitcell_json)
    except IOError:
        raise FileNotFoundError("JSON of unitcell was not found.")

    try:
        perfect_dft_data = SupercellDftResults.json_load(args.perfect_json)
    except IOError:
        raise FileNotFoundError("JSON of perfect was not found.")

    if args.dump_ewald_json:
        print("optimizing ewald...")
        ewald_kwargs = {}
        if args.ewald_init:
            ewald_kwargs["initial_value"] = args.ewald_init
        if args.ewald_convergence:
            ewald_kwargs["convergence"] = args.ewald_convergence
        if args.ewald_init:
            ewald_kwargs["prod_cutoff_fwhm"] = args.ewald_accuracy
        ewald_data = \
            Ewald.from_optimization(perfect_dft_data.final_structure,
                                    unitcell_dft_data.total_dielectric_tensor,
                                    **ewald_kwargs)
        ewald_data.to_json_file(args.dump_ewald_json)
    elif args.read_ewald_json:
        try:
            ewald_data = Ewald.load_json(args.read_ewald_json)
        except IOError:
            raise FileNotFoundError("JSON of ewald was not found.")

    if args.dir_all:
        dirs = glob('*[0-9]/')
    else:
        dirs = args.dirs

    for directory in dirs:
        print("correcting {0} ...".format(directory))
        try:
            entry = DefectEntry.json_load(directory + "/defect_entry.json")
            defect_dft_data = \
                SupercellDftResults.json_load(directory + "/dft_results.json")
            c = Correction.compute_alignment_by_extended_fnv(entry,
                                                             defect_dft_data,
                                                             perfect_dft_data,
                                                             unitcell_dft_data,
                                                             ewald_data)
        except Exception as e:
            warnings.warn("Correction for {0} is failed. "
                          "The calculation for {0} is skipped."
                          "Exception type and message is {1}, {2}".
                          format(directory, type(e), e.args))
        c.to_json_file(directory + "/correction.json")


if __name__ == "__main__":
    main()
