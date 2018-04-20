#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pydefect.input_maker.defect_initial_setting \
    import print_dopant_info, DefectInitialSetting
from pydefect.input_maker.vasp_input_maker \
    import make_incar, make_kpoints, VaspDefectInputSetMaker


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


    try:
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError:
        # argcomplete not present.
        pass

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


if __name__ == "__main__":
    main()
