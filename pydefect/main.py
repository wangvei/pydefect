#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from glob import glob
from inspect import signature
import os
from os.path import join
import shutil
import warnings

from pymatgen.core.structure import Structure

from pydefect.input_maker.defect_initial_setting \
    import print_dopant_info, DefectInitialSetting
from pydefect.input_maker.supercell_maker import Supercells
from pydefect.core.defect import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.corrections.corrections import Ewald, ExtendedFnvCorrection

from obadb.analyzer.chempotdiag.chem_pot_diag import ChemPotDiag
from obadb.vasp.input_set import ObaSet

# from pydefect.analysis.chempotdiag.chem_pot_diag \
#    import ChemPotDiag
from pydefect.analysis.chempotdiag.make_inputs import make_vasp_inputs_from_mp
from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.analysis.defect_carrier_concentration import DefectConcentration
from pydefect.analysis.defect_energy_plotter import DefectEnergyPlotter
from pydefect.util.logger import get_logger
from pydefect.input_maker.defect_entry_set_maker \
    import DefectEntrySetMaker, log_is_being_removed, log_already_exist, \
    log_is_being_constructed

__author__ = "Yu Kumagai, Akira Takahashi"
__maintainer__ = "Yu Kumagai, Akira Takahashi"

logger = get_logger(__name__)


def overwrite_default_args(class_method, main_args):
    """ Use the defaults in class.classmethod

    Args:
        class_method (classmethod): classmethod. When using __init__, class is fine.
        main_args (dict): Args set by main

    Return:
        args (dict): Overwritten args by options
    """
    from inspect import _empty

    args_with_default = []
    sig = signature(class_method)

    for name, param in sig.parameters.items():
        if param.default != _empty:
            args_with_default.append(name)

    args = {}
    for a in args_with_default:
        if hasattr(main_args, a):
            if getattr(main_args, a) is not None:
                args[a] = getattr(main_args, a)

    return args


def main():
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package for first-principles point defect calculations. It 
    allows us to construct input files, parse first-principles calculation 
    results, and analyze data.""",
        epilog="""                                 
    Author: Yu Kumagai, Akira Takahashi
    Version: {}                                                                 
    Last updated: {}""".format("0.0.1", "will be inserted"),
        #   Last updated: {}""".format(__version__, __date__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #        allow_abbrev=False)

    subparsers = parser.add_subparsers()

    # -- recommend_supercell ---------------------------------------------------
    parser_recommend_supercell = subparsers.add_parser(
        name="recommend_supercell",
        description="Tools for recommendation of an optimal supercell for "
                    "defect calculations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['rs'])
    parser_recommend_supercell.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_recommend_supercell.add_argument(
        "-s", dest="sposcar", type=str, default="SPOSCAR")
    parser_recommend_supercell.add_argument(
        "-u", dest="uposcar", type=str, default="UPOSCAR")
    parser_recommend_supercell.add_argument(
        "-c", "--criterion", dest="isotropy_criterion", type=float,
        help="Isotropy criterion.")
    parser_recommend_supercell.add_argument(
        "--min_num_atoms", dest="min_num_atoms", type=int,
        help="Minimum number of atoms")
    parser_recommend_supercell.add_argument(
        "--max_num_atoms", dest="max_num_atoms", type=int,
        help="Maximum number of atoms")
    parser_recommend_supercell.add_argument(
        "-pr", "--primitive", dest="primitive", action="store_true",
        help="Set when the supercell is expanded based on the primitive cell.")
    parser_recommend_supercell.add_argument(
        "-i", "--most_isotropic", dest="most_isotropic",
        action="store_true",
        help="Output the smallest criterion supercell instead of the smallest "
             "supercell.")
    parser_recommend_supercell.add_argument(
        "-set", dest="set", action="store_true",
        help="Output all the supercells satisfying the criterion.")

    parser_recommend_supercell.set_defaults(func=recommend_supercell)

    # -- initial_setting -------------------------------------------------------
    parser_initial = subparsers.add_parser(
        name="initial_setting",
        description="Tools for configuring initial settings for a set of "
                    "defect calculations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['is'])

    parser_initial.add_argument(
        "-p", "--poscar", dest="poscar", default="SPOSCAR", type=str,
        help="POSCAR-type file name for the supercell.")
    parser_initial.add_argument(
        "-d", "--dopants", dest="dopants", nargs="+", type=str,
        help="Dopant elements, e.g., Ga In.")
    parser_initial.add_argument(
        "-a", "--antisite", dest="is_antisite", action="store_false",
        help="Set if antisite defects are not considered.")
    parser_initial.add_argument(
        "-e", dest="en_diff", type=float,
        help="Criterion of the electronegativity_difference that determines "
             "antisites and/or substituted impurities.")
    parser_initial.add_argument(
        "--included", dest="included", type=str, nargs="+",
        help="Exceptionally included defects. E.g., Va_O2_-1.")
    parser_initial.add_argument(
        "--excluded", dest="excluded", type=str, nargs="+",
        help="Exceptionally excluded defects. E.g., Va_O2_0.")
    parser_initial.add_argument(
        "--displacement_distance", dest="displacement_distance", type=float,
        help="Displacement distance. 0 means that random displacement is not "
             "considered.")
    parser_initial.add_argument(
        "--cutoff", dest="cutoff", type=float,
        help="Set the cutoff radius [A] in which atoms are displaced.")
    parser_initial.add_argument(
        "--symprec", dest="symprec", type=float,
        help="Set precision used for symmetry analysis [A].")
    parser_initial.add_argument(
        "--angle_tol", dest="angle_tolerance", type=float,
        help="Set precision used for symmetry analysis.")
    parser_initial.add_argument(
        "--print_dopant", dest="print_dopant", type=str,
        help="Print dopant information that can be added a posteriori.")
    #    groups = parser_input.add_mutually_exclusive_group(required=True)
    parser_initial.set_defaults(func=initial_setting)

    # -- interstitial ------------------------------------------------------
    parser_interstitial = subparsers.add_parser(
        name="interstitial",
        description="Tools for managing the interstitial sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['i'])

    parser_interstitial.add_argument(
        "--yaml", dest="yaml", type=str, default="interstitials.yaml",
        help="defect_entry.yaml-type file name to be read.")
    parser_interstitial.add_argument(
        "--dposcar", dest="dposcar", default="DPOSCAR", type=str,
        help="DPOSCAR-type file name.")
    parser_interstitial.add_argument(
        "-c", dest="interstitial_coords", nargs="+", type=float,
        help="Interstitial coordinates. Eg., 0.5 0.5 0.5.")
    parser_interstitial.add_argument(
        "--name", dest="site_name", type=str, default=None,
        help="Set the interstitial site name.")
    parser_interstitial.add_argument(
        "--radius", dest="r", type=str,
        help=".")
    parser_interstitial.add_argument(
        "--force_add", dest="force_add", action="store_true",
        help="Set if the interstitial site is added, although it is too close "
             "to other atoms.")

    parser_interstitial.set_defaults(func=interstitial)

    # -- vasp_set_maker -------------------------------------------------
    parser_vasp_set = subparsers.add_parser(
        name="vasp_set",
        description="Tools for configuring vasp defect_set files for a set of "
                    "defect calculations. One needs to set "
                    ".pydefect.yaml for potcar setup.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vs'])

    parser_vasp_set.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_vasp_set.add_argument(
        "--dposcar", dest="dposcar", default="DPOSCAR", type=str,
        help="DPOSCAR-type file name.")
    parser_vasp_set.add_argument(
        "-kw", "--keywords", dest="keywords", type=str, default=None,
        nargs="+", help="Filtering keywords.")
    parser_vasp_set.add_argument(
        "-d", dest="particular_defects", type=str, default=None, nargs="+",
        help="Particular defect names to be added.")
    parser_vasp_set.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Set if the existing folders are overwritten.")
    parser_vasp_set.add_argument(
        "-x", "--xc", dest="xc", default="pbe", type=str,
        help="XC interaction treatment.")
    parser_vasp_set.add_argument(
        "-k", "--kpt_density", dest="kpt_density", default=2.3, type=float,
        help="K-points density.")
    parser_vasp_set.set_defaults(func=vasp_set)

    # -- defect_entry ----------------------------------------------------------
    parser_defect_entry = subparsers.add_parser(
        name="defect_entry",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    parser_defect_entry.add_argument(
        "--make_defect_entry", dest="make_defect_entry", action="store_true",
        help="Make defect_entry.json from yaml file.")
    parser_defect_entry.add_argument(
        "--yaml", dest="yaml", type=str, default="defect_entry.yaml",
        help="defect_entry.yaml-type file name to be read.")
    parser_defect_entry.add_argument(
        "--json", dest="json", type=str, default="defect_entry.json",
        help="defect_entry.json type file name.")
    parser_defect_entry.add_argument(
        "--print", dest="print", action="store_true",
        help="Print DefectEntry class object information.")

    parser_defect_entry.set_defaults(func=defect_entry)

    # -- supercell_calc_results ------------------------------------------------
    parser_supercell_results = subparsers.add_parser(
        name="supercell_results",
        description="Tools for analyzing vasp supercell results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sr'])

    parser_supercell_results.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Directory names.")
    parser_supercell_results.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and " "perfect directory.")
    parser_supercell_results.add_argument(
        "-p", dest="poscar", type=str)
    parser_supercell_results.add_argument(
        "-o", dest="outcar", type=str)
    parser_supercell_results.add_argument(
        "-v", dest="vasprun", type=str)
    parser_supercell_results.add_argument(
        "-pe", dest="perfect_results", type=str, default="perfect")
    parser_supercell_results.add_argument(
        "-de", dest="defect_entry_name", type=str, default="defect_entry.json")
    parser_supercell_results.add_argument(
        "--json", dest="json", type=str, default="dft_results.json",
        help="dft_results.json type file name.")
    parser_supercell_results.add_argument(
        "--print", dest="print", action="store_true",
        help="Print SupercellCalcResults class object information.")

    parser_supercell_results.set_defaults(func=supercell_calc_results)

    # -- unitcell_calc_results -------------------------------------------------
    parser_unitcell_results = subparsers.add_parser(
        name="unitcell_results",
        description="Tools for analyzing vasp unitcell results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ur'])

    parser_unitcell_results.add_argument(
        "--json_file", dest="json_file", default="unitcell.json", type=str,
        help="Json file for the unitcell info.")
    parser_unitcell_results.add_argument(
        "--static_diele", dest="static_diele", default=None, type=float,
        nargs="+",
        help="Set static dielectric constant")
    parser_unitcell_results.add_argument(
        "--ionic_diele", dest="ionic_diele", default=None, type=float,
        nargs="+",
        help="Set ionic dielectric constant")

    parser_unitcell_results.add_argument(
        "--band_edge_dir", dest="band_edge_dir", default=None, type=str,
        help="Set band edge from a vasprun.xml file")

    parser_unitcell_results.add_argument(
        "--static_diele_dir", dest="static_diele_dir", default=None, type=str,
        help="Set static dielectric constant from an OUTCAR file")
    parser_unitcell_results.add_argument(
        "--ionic_diele_dir", dest="ionic_diele_dir", default=None, type=str,
        help="Set ionic dielectric constant from an OUTCAR file")

    parser_unitcell_results.add_argument(
        "--volume_dir", dest="volume_dir", default=None, type=str,
        help="Set volume from a POSCAR file")

    parser_unitcell_results.add_argument(
        "--total_dos_dir", dest="total_dos_dir", default=None, type=str,
        help="Set total density of states from a vasprun.xml file")
    parser_unitcell_results.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_unitcell_results.add_argument(
        "-o", dest="outcar", type=str, default="OUTCAR")
    parser_unitcell_results.add_argument(
        "-v", dest="vasprun", type=str, default="vasprun.xml")
    parser_unitcell_results.add_argument(
        "--print", dest="print", action="store_true",
        help="Print Unitcell class object information.")

    parser_unitcell_results.set_defaults(func=unitcell_calc_results)

    # -- correction ------------------------------------------------------------
    parser_correction = subparsers.add_parser(
        name="correction",
        description="Tools for correction of error of defect formation energy"
                    " due to finite cell size.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['c'])

    # needed files
    parser_correction.add_argument(
        "--unitcell_json", dest="unitcell_json", default="unitcell.json",
        type=str)
    parser_correction.add_argument(
        "--perfect_json", dest="perfect_json",
        default="perfect/dft_results.json", type=str)
    # ewald
    parser_correction.add_argument(
        "--read_ewald_json", dest="read_ewald_json", default="ewald.json",
        type=str)
    parser_correction.add_argument(
        "--dump_ewald_json", dest="dump_ewald_json", default="ewald.json",
        type=str)
    parser_correction.add_argument(
        "--ewald_init", dest="ewald_init", default=None)
    parser_correction.add_argument(
        "--ewald_convergence", dest="ewald_convergence", default=None)
    parser_correction.add_argument(
        "--ewald_accuracy", dest="ewald_accuracy", default=None)
    # correction
    parser_correction.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str, help="Directory names.")
    parser_correction.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and " "perfect directory.")
    parser_correction.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Overwrite already existing correction.json.")

    parser_correction.set_defaults(func=correction)

    # -- chempotdiag -----------------------------------------------------------
    parser_chempotdiag = subparsers.add_parser(
        name="chempotdiag",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cpd'])

    # get poscar from materials project
    parser_chempotdiag.add_argument("-m", "--mat_proj_poscar",
                                    help="",
                                    action="store_true")
    parser_chempotdiag.add_argument("-el", "--elements",
                                    dest="elements", type=str, nargs='+',
                                    default=None,
                                    help="")
    parser_chempotdiag.add_argument("-dp", "--dir_path",
                                    dest="dir_path", type=str,
                                    default=None,
                                    help="")
    parser_chempotdiag.add_argument("-ch", "--criterion_hull",
                                    dest="criterion_hull", type=float,
                                    default=None,
                                    help="Collect materials only if energy above hull is less than criterion_hull. Unit is meV/atom.")
    parser_chempotdiag.add_argument("-k", "--mp_api_key",
                                    help="",
                                    action="store_true")
    parser_chempotdiag.add_argument("-gp", "--gets_poly",
                                    help="",
                                    action="store_true")

    # input
    parser_chempotdiag.add_argument("-e", "--energy", dest="energy_file",
                                    type=str, default=None,
                                    help="Name of text file of "
                                         "energies of compounds")
    parser_chempotdiag.add_argument("-v", "--vasp_dirs",
                                    dest="vasp_dirs", type=str, nargs='+',
                                    default=None,
                                    help="Drawing diagram from specified "
                                         "directories of vasp calculations")
    parser_chempotdiag.add_argument("-p", "--poscar_name",
                                    dest="poscar_name", type=str,
                                    default="POSCAR",
                                    help="Name of POSCAR, like CONTCAR, "
                                         "POSCAR-finish,...")
    parser_chempotdiag.add_argument("-o", "--outcar_name",
                                    dest="outcar_name", type=str,
                                    default="OUTCAR",
                                    help="Name of OUTCAR, like OUTCAR-finish")
    parser_chempotdiag.add_argument("-es", "--energy_shift", type=str,
                                    dest="energy_shift",
                                    nargs='+', default=None,
                                    help="Energy shift, "
                                         "e.g. -es N2/molecule 1 "
                                         "-> make more unstable N2/molecule "
                                         "by 1 eV")

    # thermodynamic status (P and T) input
    parser_chempotdiag.add_argument("-pp", "--partial_pressures",
                                    dest="partial_pressures", type=str,
                                    nargs='+', default=None,
                                    help="partial pressure of system. "
                                         "e.g. -pp O2 1e+5 N2 20000 "
                                         "-> O2: 1e+5(Pa), N2: 20000(Pa)")
    parser_chempotdiag.add_argument("-t", "--temperature",
                                    dest="temperature", type=float,
                                    default=298.15,
                                    help="temperature of system (unit: K)"
                                         "e.g. -t 3000 -> 3000(K)")

    # drawing diagram
    parser_chempotdiag.add_argument("-w", "--without_label",
                                    help="Draw diagram without label.",
                                    action="store_true")
    parser_chempotdiag.add_argument("-c", "--remarked_compound",
                                    dest="remarked_compound", type=str,
                                    default=None,
                                    help="Name of compound you are remarking."
                                         "Outputted equilibrium_points are "
                                         "limited to neighboring that "
                                         "compounds, and those equilibrium "
                                         "points are labeled in "
                                         "chem_pot_diagram.")
    parser_chempotdiag.add_argument("-d", "--draw_range",
                                    dest="draw_range", type=float,
                                    default=None,
                                    help="Drawing range of diagram."
                                         "If range is shallower than the "
                                         "deepest vertex, "
                                         "ValueError will occur")

    # output
    parser_chempotdiag.add_argument("-s", "--save_file",
                                    dest="save_file", type=str,
                                    default=None,
                                    help="File name to save the drawn diagram.")
    parser_chempotdiag.add_argument("-y", "--yaml",
                                    action="store_const", const=True,
                                    default=False,
                                    help="Dumps yaml of remarked_compound")

    parser_chempotdiag.set_defaults(func=chempotdiag)

    # -- plot_energy -----------------------------------------------------------
    parser_plot_energy = subparsers.add_parser(
        name="plot_energy",
        description="Tools for plotting defect formation energies as a "
                    "function of Fermi level",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pe'])

    parser_plot_energy.add_argument(
        "--name", dest="name", type=str, default="",
        help="System name that is written in the title.")
    parser_plot_energy.add_argument(
        "-x", "--xrange", dest="x_range", type=float, nargs='+', default=None,
        help="Two float values for the x-range of the plot.")
    parser_plot_energy.add_argument(
        "-y", "--yrange", dest="y_range", type=float, nargs='+', default=None,
        help="Two float values for the y-range of the plot.")
    parser_plot_energy.add_argument(
        "-s", "--save_file", dest="save_file", type=str, default=None,
        help="PDF file name to save the plot.")
    parser_plot_energy.add_argument(
        "--unitcell", dest="unitcell", type=str, default="unitcell.json",
        help="UnitcellCalcResults class object json file name.")
    parser_plot_energy.add_argument(
        "--perfect", dest="perfect", type=str,
        default="perfect/dft_results.json",
        help="Json file name for the SupercellCalcResults class object of the "
             "perfect supercell result.")
    parser_plot_energy.add_argument(
        "--defect_dirs", dest="defect_dirs", nargs="+", type=str,
        help="Directory names for the defect supercell results. "
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in each directory.")
    parser_plot_energy.add_argument(
        "--chem_pot_yaml", dest="chem_pot_yaml", type=str,
        default="chem_pot.yaml",
        help="Yaml file name for the chemical potential.")
    parser_plot_energy.add_argument(
        "--chem_pot_label", dest="chem_pot_label", type=str, default="A",
        help="Label indicating the equilibrium point in the chemical potential"
             "diagram.")
    parser_plot_energy.add_argument(
        "--filtering", dest="filtering", type=str, help="Filtering word.")
    parser_plot_energy.add_argument(
        "-c", "--concentration", dest="concentration", action="store_true",
        help="Calculate the carrier and defect concentrations.")
    parser_plot_energy.add_argument(
        "-stl", "--show_tl", dest="show_tl", action="store_true",
        help="Show the transition levels.")
    parser_plot_energy.add_argument(
        "-sal", "--show_al", dest="show_al", action="store_true",
        help="Show all the energy lines.")
    parser_plot_energy.add_argument(
        "-t", "--temperature", dest="temperature", nargs="+", type=float,
        help="Temperature for calculating the Fermi level. When two "
             "temperatures are supplied, the first temperature is quenched to "
             "the second temperature.")
    parser_plot_energy.add_argument(
        "-ns", "--num_sites", dest="num_site_file", type=str,
        help="The yaml file name that shows the number of sites. An example is "
             "Va_Mg1: 2")
    # parser_plot_energy.add_argument(
    #     "-u", dest="u", action="store_true",
    #     help="Calculate the U values at given defect name and three charges")

    parser_plot_energy.set_defaults(func=plot_energy)

    args = parser.parse_args()
    args.func(args)


def initial_setting(args):
    if args.print_dopant:
        print_dopant_info(args.print_dopant)
    else:
        structure = Structure.from_file(args.poscar)
        with open(args.poscar) as f:
            first_line = f.readline()
            transformation_matrix = [int(x) for x in
                                     first_line.split(":")[1].split(",")[
                                         0].split()]
            cell_multiplicity = int(
                first_line.split(":")[2].split(",")[0].split()[0])

        overwritten_args = \
            overwrite_default_args(DefectInitialSetting.from_basic_settings,
                                   args)
        defect_setting = \
            DefectInitialSetting.from_basic_settings(structure,
                                                     transformation_matrix,
                                                     cell_multiplicity,
                                                     **overwritten_args)
        defect_setting.to()


def interstitial(args):
    from pydefect.core.interstitial_site import InterstitialSiteSet

    try:
        interstitial_set = InterstitialSiteSet.from_files(args.dposcar,
                                                          args.yaml)
    except FileNotFoundError:
        structure = Structure.from_file(args.dposcar)
        interstitial_set = InterstitialSiteSet(structure=structure)

    coords = args.interstitial_coords
    interstitial_set.add_site(coords, args.site_name)

    interstitial_set.site_set_to_yaml_file(filename=args.yaml)


def vasp_set(args):
    def make_dir(name, obrs):
        if args.force_overwrite:
            if os.path.exists(name):
                log_is_being_removed(name)
                shutil.rmtree(name)

        if os.path.exists(name):
            log_already_exist(name)
        else:
            log_is_being_constructed(name)
            os.makedirs(name)
            obrs.write_input(name)

    dis = DefectInitialSetting.from_defect_in(
        poscar=args.dposcar, defect_in_file=args.defect_in)
    desm = DefectEntrySetMaker(dis, args.keywords, args.particular_defects)

    if not args.particular_defects:
        obrs = ObaSet.make_input(desm.perfect_structure,
                                 standardize_structure=False,
                                 task="defect",
                                 xc=args.xc,
                                 sort_structure=False,
                                 kpt_mode="manual",
                                 kpt_density=args.kpt_density,
                                 user_incar_settings={"ISPIN": 1})

        make_dir("perfect", obrs)

    for de in desm.defect_entries:
        defect_name = "_".join([de.name, str(de.charge)])
        json_file_name = os.path.join(defect_name, "defect_entry.json")

        obrs = ObaSet.make_input(de.perturbed_initial_structure,
                                 charge=de.charge,
                                 standardize_structure=False,
                                 task="defect",
                                 xc=args.xc,
                                 sort_structure=False,
                                 kpt_mode="manual",
                                 kpt_density=args.kpt_density)

        make_dir(defect_name, obrs)
        de.to_json_file(json_file_name)

        if de.perturbed_sites:
            poscar_name = os.path.join(defect_name, "POSCAR")
            with open(poscar_name, "r") as f:
                lines = f.readlines()
                for index, line in enumerate(lines):
                    if index - 8 in de.perturbed_sites:
                        lines[index] = line.strip() + "  Disp\n"

            with open(poscar_name, "w") as f:
                for line in lines:
                    f.write(line)


def recommend_supercell(args):
    structure = Structure.from_file(args.poscar)
    overwritten_args = overwrite_default_args(Supercells, args)

    if args.primitive:
        overwritten_args["is_conventional"] = False

    s = Supercells(structure, **overwritten_args)

    if s.supercells:
        if args.set:
            logger.info("The number of supercells:" + str(len(s.supercells)))
            for supercell in s.supercells:
                cell = "c" if s.is_conventional_based else "p"
                multi_str = cell + "x".join([str(list(supercell.trans_mat)[i])
                                             for i in range(3)])
                name = "_".join([args.sposcar,
                                 multi_str,
                                 str(supercell.num_atoms),
                                 str(supercell.isotropy)])
                supercell.to_poscar(poscar_filename=name)
            s.supercells[0].to_uposcar(uposcar_filename=args.uposcar)

        else:
            if args.most_isotropic:
                supercell = s.create_most_isotropic_supercell
            else:
                supercell = s.create_smallest_supercell
            supercell.to_poscar(poscar_filename=args.sposcar)
            supercell.to_uposcar(uposcar_filename=args.uposcar)

    else:
        logger.warning("Any supercell does not satisfy the criterion. Increase "
                       "the criterion if needed.")


def defect_entry(args):
    if args.make_defect_entry:
        defect_entry_from_yaml = DefectEntry.from_yaml(args.yaml)
        defect_entry_from_yaml.to_json_file("defect_entry.json")
    elif args.print:
        print(DefectEntry.load_json(args.json))


def supercell_calc_results(args):
    if args.print:
        print(SupercellCalcResults.load_json(args.json))

    if args.dir_all:
        dirs = glob('*[0-9]/')
        dirs.append("perfect/")
    else:
        dirs = args.dirs

    for d in dirs:
        overwritten_args = \
            overwrite_default_args(SupercellCalcResults.from_vasp_files, args)

        if os.path.isdir(d):
            logger.info("Parsing data in {}...".format(d))
#            try:
            if d in ["perfect", "perfect/"]:
                dft_results = \
                    SupercellCalcResults.from_vasp_files(d,
                                                         **overwritten_args)
            else:
                filename = join(args.perfect_results, "dft_results.json")
                perfect_results = SupercellCalcResults.load_json(filename)
                de = DefectEntry.load_json(join(d, args.defect_entry_name))

                dft_results = \
                    SupercellCalcResults.from_vasp_files(
                        directory_path=d,
                        defect_entry=de,
                        referenced_dft_results=perfect_results,
                        **overwritten_args)

            dft_results.to_json_file(
                filename=join(d, "dft_results.json"))
#            except ValueError:
#                logger.warning("Parsing data in {} is failed.".format(d))
        else:
            logger.warning("{} does not exist, so nothing is done.".format(d))


def unitcell_calc_results(args):
    try:
        dft_results = UnitcellCalcResults.load_json(filename=args.json_file)
    except:
        dft_results = UnitcellCalcResults()

    if args.print:
        dft_results = UnitcellCalcResults.load_json(filename=args.json_file)
        print(dft_results.total_dielectric_tensor)
        print(dft_results)
        return

    if args.band_edge_dir:
        try:
            dft_results.set_band_edge_from_vasp(args.band_edge_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            raise FileNotFoundError(args.band_edge_dir, "is not appropriate.")

    if args.static_diele:
        dft_results.static_dielectric_tensor = args.static_diele
    elif args.static_diele_dir:
        try:
            dft_results.set_static_dielectric_tensor_from_vasp(
                args.static_diele_dir, outcar_name=args.outcar)
        except IOError:
            raise FileNotFoundError(
                args.static_diele_dir, "is not appropriate.")

    if args.ionic_diele:
        dft_results.ionic_dielectric_tensor = args.ionic_diele
    elif args.ionic_diele_dir:
        try:
            dft_results.set_ionic_dielectric_tensor_from_vasp(
                args.ionic_diele_dir, outcar_name=args.outcar)
        except IOError:
            raise FileNotFoundError(args.ionic_diele_dir, "is not appropriate.")

    if args.volume_dir:
        try:
            dft_results.set_volume_from_vasp(
                args.volume_dir, contcar_name=args.poscar)
        except IOError:
            raise FileNotFoundError(args.volume_dir, "is not appropriate.")

    if args.total_dos_dir:
        try:
            dft_results.set_total_dos_from_vasp(args.total_dos_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            raise FileNotFoundError(args.total_dos_dir, "is not appropriate.")

    dft_results.to_json_file(args.json_file)


def correction(args):
    try:
        unitcell_dft_data = UnitcellCalcResults.load_json(args.unitcell_json)
    except IOError:
        raise FileNotFoundError("JSON for the unitcell info is not found.")

    try:
        perfect_dft_data = SupercellCalcResults.load_json(args.perfect_json)
        print(type(perfect_dft_data))
    except IOError:
        raise FileNotFoundError("JSON for the perfect supercell is not found.")

    if args.dump_ewald_json:
        logger.info("optimizing ewald...")
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
        ewald_filename = args.dump_ewald_json
    elif args.read_ewald_json:
        try:
            ewald_filename = args.read_ewald_json
        except IOError:
            raise FileNotFoundError("JSON for ewald parameter is not found.")
    else:
        raise IOError("Ewald parameter is a mandatory.")

    dirs = glob('*[0-9]/') if args.dir_all else args.dirs

    for directory in dirs:
        json_to_make = join(directory, "correction.json")

        if os.path.exists(json_to_make) and not args.force_overwrite:
            logger.warning("{} exists. ExtendedFnvCorrection was not done."
                           .format(json_to_make))
            continue

        logger.info("correcting {} ...".format(directory))
        entry = DefectEntry.load_json(join(directory, "defect_entry.json"))
        defect_dft_data = SupercellCalcResults.load_json(
            join(directory, "dft_results.json"))

        c = ExtendedFnvCorrection. \
            compute_alignment_by_extended_fnv(entry,
                                              defect_dft_data,
                                              unitcell_dft_data,
                                              ewald_filename)

        c.plot_distance_vs_potential(join(directory, "potential.eps"))
        c.to_json_file(join(directory, "correction.json"))


def chempotdiag(args):
    if args.mat_proj_poscar:
        kwargs_to_make_vasp_inputs = {}
        if args.dir_path:
            kwargs_to_make_vasp_inputs["dir_path"] = args.dir_path
        if args.criterion_hull:
            kwargs_to_make_vasp_inputs["criterion_e_above_hull"] \
                = args.criterion_hull
        if args.mp_api_key:
            kwargs_to_make_vasp_inputs["api_key"] = args.mp_api_key
        if args.gets_poly:
            kwargs_to_make_vasp_inputs["gets_poly"] = True
        make_vasp_inputs_from_mp(args.elements, **kwargs_to_make_vasp_inputs)
    else:
        if args.energy_file and args.vasp_dirs:
            raise ValueError("You can not specify energy_file and vasp_dirs "
                             "simultaneously.")
        if args.energy_file:
            if args.temperature or args.partial_pressures:
                warnings.warn("Now temperature and pressures can not apply when"
                              " reading data from energy_file")
            cp = ChemPotDiag.from_file(args.energy_file)
        if args.vasp_dirs:

            poscar_paths = [d + args.poscar_name for d in args.vasp_dirs]
            outcar_paths = [d + args.outcar_name for d in args.vasp_dirs]

            # pressure and temperature
            partial_pressure_dict = {}
            if args.partial_pressures:
                if len(args.partial_pressures) % 2 != 0:
                    raise ValueError("Invalid partial pressures input {}".
                                     format(args.partial_pressures))
                for i in range(int(len(args.partial_pressures) / 2)):
                    formula = args.partial_pressures[2 * i]
                    pressure = args.partial_pressures[2 * i + 1]
                    partial_pressure_dict[formula] = float(pressure)

            # manually set energy
            energy_shift_dict = {}
            if args.energy_shift:
                if len(args.energy_shift) % 2 != 0:
                    raise ValueError("Invalid energy shift input {}".
                                     format(args.energy_shift))
                for i in range(int(len(args.energy_shift) / 2)):
                    output_name = \
                        args.energy_shift[2 * i] + "/" + args.outcar_name
                    es = args.energy_shift[2 * i + 1]
                    energy_shift_dict[output_name] = float(es)

            cp = ChemPotDiag. \
                from_vasp_calculations_files(poscar_paths,
                                             outcar_paths,
                                             temperature=args.temperature,
                                             pressure=partial_pressure_dict,
                                             energy_shift_dict=energy_shift_dict)
            if args.elements:
                cp.set_elements(args.elements)
        print("Energies of elements ({0}) : {1}"
              .format(cp.elements, cp.element_energy))
        #  Read args of drawing diagram from parser
        if args.remarked_compound:
            try:
                for vertex in cp.get_neighbor_vertices(args.remarked_compound):
                    print(vertex)
            except ValueError:
                print("{0} is unstable. No vertex is labeled."
                      .format(args.remarked_compound))

        kwargs_for_diagram = {}
        if args.remarked_compound:
            kwargs_for_diagram["remarked_compound"] = args.remarked_compound
        if args.save_file:
            kwargs_for_diagram["save_file_name"] = args.save_file
        if args.without_label:
            kwargs_for_diagram["with_label"] = False
        if args.draw_range:
            kwargs_for_diagram["draw_range"] = args.draw_range

        if cp.dim >= 4:
            print("Currently diagram is not available for quaternary or more.")
        else:
            try:
                cp.draw_diagram(**kwargs_for_diagram)
            except ValueError:
                kwargs_for_diagram.pop("remarked_compound")
                cp.draw_diagram(**kwargs_for_diagram)

        if args.yaml:
            if args.remarked_compound is None:
                raise ValueError("remarked_compound is needed to dump yaml")
            cp.dump_vertices_yaml(os.getcwd(), args.remarked_compound)


def plot_energy(args):
    unitcell = UnitcellCalcResults.load_json(args.unitcell)
    perfect = SupercellCalcResults.load_json(args.perfect)

    if not args.defect_dirs:
        defects_dirs = glob('*[0-9]/')
    else:
        defects_dirs = args.defect_dirs

    defects = []
    for d in defects_dirs:
        try:
            e = DefectEntry.load_json(join(d, "defect_entry.json"))
            r = SupercellCalcResults.load_json(join(d, "dft_results.json"))
            c = ExtendedFnvCorrection.load_json(join(d, "correction.json"))
            defects.append(Defect(defect_entry=e, dft_results=r, correction=c))
        except:
            warnings.warn(message="Parsing data in " + d + " is failed.")

    chem_pot = ChemPotDiag.load_vertices_yaml(args.chem_pot_yaml)

    # First construct DefectEnergies class object.
    defect_energies = \
        DefectEnergies.from_files(unitcell=unitcell,
                                  perfect=perfect,
                                  defects=defects,
                                  chem_pot=chem_pot,
                                  chem_pot_label=args.chem_pot_label,
                                  system=args.name)

    if args.concentration:
        defect_concentration = \
            DefectConcentration.from_defect_energies(
                defect_energies=defect_energies,
                temperature=args.temperature[0],
                unitcell=unitcell,
                num_sites_filename=args.num_site_file)

        if len(args.temperature) == 2:
            defect_concentration = \
                DefectConcentration.from_defect_energies(
                    defect_energies=defect_energies,
                    temperature=args.temperature[1],
                    unitcell=unitcell,
                    num_sites_filename=args.num_site_file,
                    previous_concentration=defect_concentration)
    else:
        defect_concentration = None

    plotter = DefectEnergyPlotter(defect_energies, defect_concentration)
    plt = plotter.plot_energy(filtering_words=args.filtering,
                              x_range=args.x_range,
                              y_range=args.y_range,
                              show_fermi_level=args.concentration,
                              show_transition_levels=args.show_tl,
                              show_all_lines=args.show_al)

    if args.save_file:
        plt.savefig(args.save_file, format="pdf")
    else:
        plt.show()


if __name__ == "__main__":
    main()
