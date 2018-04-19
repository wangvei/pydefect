#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pydefect.input_maker.defect_initial_setting \
    import print_dopant_info, DefectInitialSetting

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

    parser_initial = subparsers.add_parser(
        name="initial_setting",
        description="Tools for configuring initial settings for a set of "
                    "defect calculations. One need to set .pydefect.yaml"
                    "for potcar setup.",
        aliases=['initial', 'in'])

    parser_initial.add_argument(
        "-p", "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR name for the unitcell.")
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


if __name__ == "__main__":
    main()
