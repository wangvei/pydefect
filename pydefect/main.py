#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import warnings
from xml.etree.cElementTree import ParseError

from copy import deepcopy
from glob import glob
from inspect import signature
from itertools import chain
from os.path import join

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

from obadb.analyzer.chempotdiag.chem_pot_diag import ChemPotDiag
from obadb.vasp.input_set import ObaSet
from obadb.vasp.incar import incar_flags

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.analysis.defect_energy_plotter import DefectEnergyPlotter
from pydefect.analysis.defect_eigenvalues import DefectEigenvalue
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.core.prior_info import PriorInfo
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.corrections.corrections \
    import Ewald, NoCorrection, ExtendedFnvCorrection
from pydefect.input_maker.defect_entry_set_maker \
    import DefectEntrySetMaker, log_is_being_removed, log_already_exist, \
    log_is_being_constructed
from pydefect.input_maker.defect_initial_setting \
    import dopant_info, DefectInitialSetting
from pydefect.input_maker.supercell_maker import Supercells
from pydefect.util.logger import get_logger
from pydefect.util.main_tools import get_default_args, list2dict


__author__ = "Yu Kumagai, Akira Takahashi"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package for first-principles point defect calculations. It 
    allows one to construct input files, parse first-principles calculation 
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

    defaults = get_default_args(Supercells)

    parser_recommend_supercell.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_recommend_supercell.add_argument(
        "-s", dest="sposcar", type=str, default="SPOSCAR")
    parser_recommend_supercell.add_argument(
        "-u", dest="uposcar", type=str, default="UPOSCAR")
    parser_recommend_supercell.add_argument(
        "-c", "--criterion", dest="isotropy_criterion", type=float,
        default=defaults["isotropy_criterion"],
        help="Isotropy criterion.")
    parser_recommend_supercell.add_argument(
        "--min_num_atoms", dest="min_num_atoms", type=int,
        default=defaults["min_num_atoms"],
        help="Minimum number of atoms")
    parser_recommend_supercell.add_argument(
        "--max_num_atoms", dest="max_num_atoms", type=int,
        default=defaults["max_num_atoms"],
        help="Maximum number of atoms")
    parser_recommend_supercell.add_argument(
        "-pr", "--primitive", dest="primitive", action="store_true",
        help="Set when the supercell is expanded based on the primitive cell.")
    parser_recommend_supercell.add_argument(
        "-i", "--most_isotropic", dest="most_isotropic", action="store_true",
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

    defaults = get_default_args(DefectInitialSetting.from_basic_settings)

    parser_initial.add_argument(
        "-p", "--poscar", dest="poscar", default="SPOSCAR", type=str,
        help="POSCAR-type file name for the supercell.")
    parser_initial.add_argument(
        "-d", "--dopants", dest="dopants", nargs="+", type=str,
        default=defaults["dopants"],
        help="Dopant elements, e.g., Ga In.")
    parser_initial.add_argument(
        "-a", "--antisite", dest="is_antisite", action="store_false",
        default=defaults["is_antisite"],
        help="Set if antisite defects are not considered.")
    parser_initial.add_argument(
        "-e", dest="en_diff", type=float, default=defaults["en_diff"],
        help="Criterion of the electronegativity_difference that determines "
             "antisites and/or substituted impurities.")
    parser_initial.add_argument(
        "--included", dest="included", type=str, nargs="+",
        default=defaults["included"],
        help="Exceptionally included defects. E.g., Va_O2_-1.")
    parser_initial.add_argument(
        "--excluded", dest="excluded", type=str, nargs="+",
        default=defaults["excluded"],
        help="Exceptionally excluded defects. E.g., Va_O2_0.")
    parser_initial.add_argument(
        "--displacement_distance", dest="displacement_distance", type=float,
        default=defaults["displacement_distance"],
        help="Displacement distance. 0 means that random displacement is not "
             "considered.")
    parser_initial.add_argument(
        "--cutoff", dest="cutoff", type=float,
        default=defaults["cutoff"],
        help="Set the cutoff radius [A] in which atoms are displaced.")
    parser_initial.add_argument(
        "--symprec", dest="symprec", type=float,
        default=defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_initial.add_argument(
        "--angle_tol", dest="angle_tolerance", type=float,
        default=defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
    parser_initial.add_argument(
        "--interstitial_sites", dest="interstitials", type=str, nargs="+",
        default=defaults["interstitial_site_names"],
        help="Interstitial site names.")
    parser_initial.add_argument(
        "--print_dopant", dest="print_dopant", type=str,
        help="Print dopant information that can be added a posteriori.")

    parser_initial.set_defaults(func=initial_setting)

    # -- interstitial ------------------------------------------------------
    parser_interstitial = subparsers.add_parser(
        name="interstitial",
        description="Tools for handling the interstitial sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['i'])

    defaults = get_default_args(InterstitialSiteSet.add_site)
    defaults.update(get_default_args(InterstitialSiteSet.from_files))

    parser_interstitial.add_argument(
        "--yaml", dest="yaml", type=str, default=defaults["filename"],
        help="defect_entry.yaml-type file name to be read.")
    parser_interstitial.add_argument(
        "--dposcar", dest="dposcar", type=str, default=defaults["structure"],
        help="DPOSCAR-type file name.")
    parser_interstitial.add_argument(
        "-c", dest="interstitial_coords", nargs="+", type=float,
        help="Interstitial coordinates. Eg., 0.5 0.5 0.5.")
    parser_interstitial.add_argument(
        "--name", dest="site_name", type=str, default=None,
        help="Set the interstitial site name.")
    parser_interstitial.add_argument(
        "--radius", dest="radius", type=str,
        default=defaults["check_neighbor_radius"],
        help="Radius for checking too closed neighbor atoms.")
    parser_interstitial.add_argument(
        "--force_add", dest="force_add", action="store_true",
        help="Set if the interstitial site is added, although it is too close "
             "to other atoms.")
    parser_interstitial.add_argument(
        "--symprec", dest="symprec", type=float,
        default=defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_interstitial.add_argument(
        "--angle_tol", dest="angle_tolerance", type=float,
        default=defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
    parser_interstitial.add_argument(
        "--method", dest="method", type=str, default=defaults["method"],
        help="Name of method.")

    parser_interstitial.set_defaults(func=interstitial)

    # -- defect_vasp_set_maker -------------------------------------------------
    parser_vasp_set = subparsers.add_parser(
        name="defect_vasp_set",
        description="Tools for configuring vasp defect_set files for a set of "
                    "defect calculations. One needs to set "
                    ".pydefect.yaml for potcar setup.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['dvs'])

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
    parser_vasp_set.add_argument(
        "-w", "--make_wavecar", dest="wavecar", default=True, type=bool,
        help="Whether to make WAVECAR file or not.")
    parser_vasp_set.set_defaults(func=defect_vasp_set)

    # -- defect_entry ----------------------------------------------------------
    parser_defect_entry = subparsers.add_parser(
        name="defect_entry",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations. By default, print defect_entry.",
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

    parser_defect_entry.set_defaults(func=defect_entry)

    # -- supercell_calc_results ------------------------------------------------
    parser_supercell_results = subparsers.add_parser(
        name="supercell_results",
        description="Tools for analyzing vasp supercell results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sr'])

    defaults = get_default_args(SupercellCalcResults.from_vasp_files)

    parser_supercell_results.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Directory names.")
    parser_supercell_results.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and " "perfect directory.")
    parser_supercell_results.add_argument(
        "-v", dest="vasprun", default=defaults["vasprun"], type=str)
    parser_supercell_results.add_argument(
        "-c", dest="contcar", default=defaults["contcar"], type=str)
    parser_supercell_results.add_argument(
        "-o", dest="outcar", default=defaults["outcar"], type=str)
    parser_supercell_results.add_argument(
        "-s", dest="not_check_shallow", action="store_false",
        help="Not check whether the defects are shallow.")
    parser_supercell_results.add_argument(
        "-pe", dest="perfect_results", type=str, default="perfect")
    parser_supercell_results.add_argument(
        "-de", dest="defect_entry_name", type=str, default="defect_entry.json")
    parser_supercell_results.add_argument(
        "-be", dest="band_edge", type=str, nargs="+", default=None)
    parser_supercell_results.add_argument(
        "--json", dest="json", type=str, default="dft_results.json",
        help="dft_results.json type file name.")
    parser_supercell_results.add_argument(
        "--symprec", dest="symprec", type=float,
        default=defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_supercell_results.add_argument(
        "--angle_tol", dest="angle_tolerance", type=float,
        default=defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
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

    # -- FNV correction --------------------------------------------------------
    parser_correction = subparsers.add_parser(
        name="extended_fnv_correction",
        description="Tools for extended FNV correction for error of defect "
                    "formation energy due to finite supercell size.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['efc'])

    # when no corrections are required for some reasons
    parser_correction.add_argument(
        "--nocorr", "-nc", dest="nocorr", action="store_true",
        help="Set when no corrections are required for some reasons.")
    parser_correction.add_argument(
        "-m", "--manual", dest="manual", type=float, default=0.0,
        help="Set manual correction energy.")
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
    parser_correction.add_argument(
        "--json_file", dest="json_file", default="correction.json", type=str,
        help="Json file for the correction.")
    parser_correction.add_argument(
        "--print", dest="print", action="store_true",
        help="Print FNV correction information.")

    parser_correction.set_defaults(func=efnv_correction)

    # -- vasp_set ----------------------------------------------------------
    parser_vasp_oba_set = subparsers.add_parser(
        name="vasp_set",
        description="Tools for constructing vasp input set with oba_set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vos'])

    parser_vasp_oba_set.add_argument(
        "-p", "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR-type file name for the supercell.")
    parser_vasp_oba_set.add_argument(
        "-x", "--xc", dest="xc", default="pbesol", type=str,
        help="XC interaction treatment.")
    parser_vasp_oba_set.add_argument(
        "-t", "--task", dest="task", default="structure_opt", type=str,
        help="The task name.")
    parser_vasp_oba_set.add_argument(
        "-k", "--kpt_density", dest="kpt_density", default=3, type=float,
        help="K-points density.")
    parser_vasp_oba_set.add_argument(
        "-s", "--standardize", dest="standardize", action="store_false",
        help="Store if one doesn't want the cell to be transformed to a "
             "primitive one.")
    parser_vasp_oba_set.add_argument(
        "-d", "--prev_dir", dest="prev_dir", type=str,
        help=".")
    # parser_vasp_oba_set.add_argument(
    #     "-gga", "--prev_gga", dest="prev_gga",
    #     action="store_true", help=".")
    # parser_vasp_oba_set.add_argument(
    #     "-gpre", "--prev_dir_gw_pre2", dest="prev_dir_gw_pre2", type=str,
    #     help=".")
    # parser_vasp_oba_set.add_argument(
    #     "-gw", "--prev_dir_gw0", dest="prev_dir_gw0", type=str,
    #     help=".")
    parser_vasp_oba_set.add_argument(
        "-kw", "--kwargs", dest="kwargs", type=str, nargs="+",
        default=None, help="keyword arguments.")
    parser_vasp_oba_set.add_argument(
        "-is", "--incar_setting", dest="incar_setting", type=str, nargs="+",
        default=None, help="INCAR setting to be overwritten.")
    parser_vasp_oba_set.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Make vasp set for the directories.")
    parser_vasp_oba_set.add_argument(
        "-pi", "-prior_info", dest="prior_info", action="store_true",
        help="Set if prior_info.json is not read.")

    parser_vasp_oba_set.set_defaults(func=vasp_set)

    # -- diagnose --------------------------------------------------------------
    parser_diagnose = subparsers.add_parser(
        name="diagnose",
        description="Tools for diagnosing results related to defects.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['d'])

    parser_diagnose.add_argument(
        "--defect_dirs", dest="defect_dirs", nargs="+", type=str,
        help="Directory names for the defect supercell results. "
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in each directory.")
    parser_diagnose.add_argument(
        "--unitcell", dest="unitcell", type=str, default="unitcell.json",
        help="UnitcellCalcResults class object json file name.")
    parser_diagnose.add_argument(
        "--perfect", dest="perfect", type=str,
        default="perfect/dft_results.json",
        help="Json file name for the SupercellCalcResults class object of the "
             "perfect supercell result.")
    parser_diagnose.add_argument(
        "--json", dest="json", type=str, default="dft_results.json",
        help="dft_results.json type file name.")

    parser_diagnose.set_defaults(func=diagnose)

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
        "--de", dest="defect_entry", type=str, default="defect_entry.json")
    parser_plot_energy.add_argument(
        "--dft_results", dest="dft_results", type=str,
        default="dft_results.json")
    parser_plot_energy.add_argument(
        "--correction", dest="correction", type=str, default="correction.json")
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
        "-stl", "--show_transition_level", dest="show_tl", action="store_true",
        help="Show the transition levels.")
    parser_plot_energy.add_argument(
        "-sa", "--show_all", dest="show_all", action="store_true",
        help="Show all the energy lines.")
    parser_plot_energy.add_argument(
        "-t", "--temperature", dest="temperature", nargs="+", type=float,
        help="Temperature for calculating the Fermi level. When two "
             "temperatures are supplied, the first temperature is quenched to "
             "the second temperature.")
    # parser_plot_energy.add_argument(
    #     "-u", dest="u", action="store_true",
    #     help="Calculate the U values at given defect name and three charges")

    parser_plot_energy.set_defaults(func=plot_energy)

    # -- parse_eigenvalues -----------------------------------------------------
    parser_parse_eigenvalues = subparsers.add_parser(
        name="parse_eigenvalues",
        description="Tools for parsing defect eigenvalues",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['eig'])

    parser_parse_eigenvalues.add_argument(
        "--name", dest="name", type=str, default="",
        help="System name that is written in the title.")
    parser_parse_eigenvalues.add_argument(
        "-y", "--yrange", dest="y_range", type=float, nargs='+', default=None,
        help="Two float values for the y-range of the plot.")
    parser_parse_eigenvalues.add_argument(
        "-s", "--save_file", dest="save_file", type=str, default=None,
        help="PDF file name to save the plot.")
    parser_parse_eigenvalues.add_argument(
        "--unitcell", dest="unitcell", type=str, default="unitcell.json",
        help="UnitcellCalcResults class object json file name.")
    parser_parse_eigenvalues.add_argument(
        "--perfect", dest="perfect", type=str,
        default="perfect/dft_results.json",
        help="Json file name for the SupercellCalcResults class object of the "
             "perfect supercell result.")
    parser_parse_eigenvalues.add_argument(
        "--de", dest="defect_entry", type=str, default="defect_entry.json")
    parser_parse_eigenvalues.add_argument(
        "--dft_results", dest="dft_results", type=str,
        default="dft_results.json")
    parser_parse_eigenvalues.add_argument(
        "--correction", dest="correction", type=str, default="correction.json")
    parser_parse_eigenvalues.add_argument(
        "--defect_dir", dest="defect_dir", type=str,
        help="Directory name for the defect supercell result. "
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in the directory.")

    parser_parse_eigenvalues.set_defaults(func=parse_eigenvalues)

    # -- vasp_parchg_set -------------------------------------------------------
    parser_vasp_parchg_set = subparsers.add_parser(
        name="vasp_parchg_set",
        description="Tools for constructing vasp set for generating PARCHG",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vps'])

    parser_vasp_parchg_set.add_argument(
        "-rd", "--read_dir", dest="read_dir", type=str,
        help="Read directory name.")
    parser_vasp_parchg_set.add_argument(
        "-wd", "--write_dir", dest="write_dir", type=str, default=".",
        help="Write directory name.")
    parser_vasp_parchg_set.add_argument(
        "-b", "--bands", dest="band_indices", type=int, nargs='+',
        help="Band indices.")
    parser_vasp_parchg_set.add_argument(
        "-k", "--kpoints", dest="kpoint_indices", type=int, nargs='+',
        help="K-point indices.")

    parser_vasp_parchg_set.set_defaults(func=vasp_parchg_set)

    args = parser.parse_args()
    args.func(args)


def recommend_supercell(args):
    s = Supercells(structure=Structure.from_file(args.poscar),
                   is_conventional=not args.primitive,
                   max_num_atoms=args.max_num_atoms,
                   min_num_atoms=args.min_num_atoms,
                   isotropy_criterion=args.isotropy_criterion)

    if s.supercells:
        if args.set:
            logger.info("The number of supercells:" + str(len(s.supercells)))

            for supercell in s.supercells:
                # Suffix "c" means conventional cell, while "p" primitive cell.
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


def initial_setting(args):
    if args.print_dopant:
        print(dopant_info(args.print_dopant))

    else:
        structure = Structure.from_file(args.poscar)
        with open(args.poscar) as f:
            first_line = f.readline()
            transformation_matrix = [int(x) for x in
                                     first_line.split(":")[1].split(",")[
                                         0].split()]
            cell_multiplicity = int(
                first_line.split(":")[2].split(",")[0].split()[0])

        defect_setting = \
            DefectInitialSetting.from_basic_settings(
                structure=structure,
                transformation_matrix=transformation_matrix,
                cell_multiplicity=cell_multiplicity,
                dopants=args.dopants,
                is_antisite=args.is_antisite,
                en_diff=args.en_diff,
                included=args.included,
                excluded=args.excluded,
                displacement_distance=args.displacement_distance,
                cutoff=args.cutoff,
                symprec=args.symprec,
                angle_tolerance=args.angle_tolerance,
                interstitial_site_names=args.interstitials)

        defect_setting.to()


def interstitial(args):
    try:
        interstitial_set = \
            InterstitialSiteSet.from_files(args.dposcar, args.yaml)
    except FileNotFoundError:
        structure = Structure.from_file(args.dposcar)
        interstitial_set = InterstitialSiteSet(structure=structure)

    coords = args.interstitial_coords
    interstitial_set.add_site(coord=coords,
                              site_name=args.site_name,
                              check_neighbor_radius=args.radius,
                              force_add=args.force_add,
                              symprec=args.symprec,
                              angle_tolerance=args.angle_tolerance)

    interstitial_set.site_set_to_yaml_file(filename=args.yaml)


def defect_vasp_set(args):

    def make_dir(name, obrs):
        """Helper function"""
        if args.force_overwrite and os.path.exists(name):
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
    defect_entry_set_maker = \
        DefectEntrySetMaker(dis, args.keywords, args.particular_defects)

    oba_set = ObaSet.make_input(
        structure=defect_entry_set_maker.perfect_structure,
        standardize_structure=False,
        task="defect",
        xc=args.xc,
        sort_structure=False,
        kpt_mode="manual",
        kpt_density=args.kpt_density,
        user_incar_settings={"ISPIN": 1})

    make_dir("perfect", oba_set)

    for de in defect_entry_set_maker.defect_entries:
        defect_name = "_".join([de.name, str(de.charge)])
        json_file_name = os.path.join(defect_name, "defect_entry.json")

        oba_set = ObaSet.make_input(structure=de.perturbed_initial_structure,
                                    charge=de.charge,
                                    standardize_structure=False,
                                    task="defect",
                                    xc=args.xc,
                                    sort_structure=False,
                                    weak_incar_settings={"LWAVE": args.wavecar},
                                    kpt_mode="manual",
                                    kpt_density=args.kpt_density)

        make_dir(defect_name, oba_set)
        de.to_json_file(json_file_name)

        if de.neighboring_sites:
            poscar_name = os.path.join(defect_name, "POSCAR")
            with open(poscar_name, "r") as f:
                lines = f.readlines()
                for index, line in enumerate(lines.copy()):
                    if index - 8 in de.neighboring_sites:
                        lines[index] = line.strip() + "  Disp\n"

            with open(poscar_name, "w") as f:
                for line in lines:
                    f.write(line)


def defect_entry(args):
    if args.make_defect_entry:
        defect_entry_from_yaml = DefectEntry.from_yaml(args.yaml)
        defect_entry_from_yaml.to_json_file(args.json)
    else:
        print(DefectEntry.load_json(args.json))


def supercell_calc_results(args):
    if args.print:
        print(SupercellCalcResults.load_json(args.json))
        return

    if args.dir_all:
        dirs = glob('*[0-9]/')
        dirs.insert(0, "perfect/")
    else:
        dirs = args.dirs

    for d in dirs:
        if os.path.isdir(d):
            logger.info("Parsing data in {}...".format(d))
            if args.band_edge:
                if args.band_edge[0] == "up":
                    spin = Spin.up
                elif args.band_edge[0] == "down":
                    spin = Spin.down
                else:
                    raise ValueError("band edge flag is inadequate. "
                                     "Ex. -be up no_in_gap")
                state = args.band_edge[1]
                dft_results = SupercellCalcResults.load_json(args.json)
                dft_results.set_band_edges(spin=spin, state=state)
            elif d in ["perfect", "perfect/"]:
                try:
                    dft_results = SupercellCalcResults.from_vasp_files(
                        directory_path=d,
                        vasprun=args.vasprun,
                        contcar=args.contcar,
                        procar=True,
                        outcar=args.outcar)
                except:
                    raise IOError("Parsing data in perfect failed")
            else:
                try:
                    filename = join(args.perfect_results, "dft_results.json")
                    perfect_results = SupercellCalcResults.load_json(filename)
                    de = DefectEntry.load_json(join(d, args.defect_entry_name))

                    dft_results = \
                        SupercellCalcResults.from_vasp_files(
                            directory_path=d,
                            vasprun=args.vasprun,
                            contcar=args.contcar,
                            outcar=args.outcar,
                            procar=True,
                            referenced_dft_results=perfect_results,
                            defect_entry=de,
                            symprec=args.symprec,
                            angle_tolerance=args.angle_tolerance)
                except ParseError:
                    logger.warning("Parsing data in {} failed.".format(d))
                    continue

            dft_results.to_json_file(
                filename=join(d, "dft_results.json"))
        else:
            logger.warning("{} does not exist, so nothing is done.".format(d))


def unitcell_calc_results(args):
    if args.print:
        print(UnitcellCalcResults.load_json(filename=args.json_file))
        return

    try:
        dft_results = UnitcellCalcResults.load_json(filename=args.json_file)
    except IOError:
        dft_results = UnitcellCalcResults()

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


def efnv_correction(args):
    if args.print:
        print(ExtendedFnvCorrection.load_json(args.json_file))
        return

    dirs = glob('*[0-9]/') if args.dir_all else args.dirs

    if args.nocorr:
        for directory in dirs:
            c = NoCorrection(manual_correction_energy=args.manual)
            c.to_json_file(join(directory, "correction.json"))

        return

    try:
        unitcell_dft_data = UnitcellCalcResults.load_json(args.unitcell_json)
    except IOError:
        raise FileNotFoundError("JSON for the unitcell info is not found.")

    try:
        perfect_dft_data = SupercellCalcResults.load_json(args.perfect_json)
    except IOError:
        raise FileNotFoundError("JSON for the perfect supercell is not found.")

    # Ewald parameter related
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
            compute_alignment_by_extended_fnv(defect_entry=entry,
                                              defect_dft=defect_dft_data,
                                              unitcell_dft=unitcell_dft_data,
                                              ewald_json=ewald_filename)

        c.plot_distance_vs_potential(join(directory, "potential.pdf"))
        c.to_json_file(join(directory, "correction.json"))


def vasp_set(args):

    #TODO: When writing GW part, refer oba_set_main.py in obadb

    base_kwargs = {"task":                  args.task,
                   "xc":                    args.xc,
                   "kpt_density":           args.kpt_density,
                   "standardize_structure": args.standardize}

    flags = list(chain.from_iterable(incar_flags.values()))
    base_user_incar_settings = list2dict(args.incar_setting, flags)

    flags = list(signature(ObaSet.make_input).parameters.keys())
    base_kwargs.update(list2dict(args.kwargs, flags))

    original_dir = os.getcwd()
    dirs = args.dirs if args.dirs else ["."]

    for d in dirs:
        os.chdir(os.path.join(original_dir, d))
        logger.info("Constructing vasp set in {}".format(d))
        user_incar_settings = deepcopy(base_user_incar_settings)
        kwargs = deepcopy(base_kwargs)

        if args.prior_info:
            if os.path.exists("prior_info.json"):
                prior_info = PriorInfo.load_json("prior_info.json")
                kwargs["band_gap"] = prior_info.band_gap
                kwargs["is_magnetic"] = prior_info.is_magnetic
                kwargs["is_cluster"] = prior_info.is_cluster

        if args.prev_dir:
            files = {"CHGCAR": "C", "WAVECAR": "L", "WAVEDER": "M"}
            oba_set = ObaSet.from_prev_calc(args.prev_dir,
                                            copied_file_names=files, **kwargs)
        else:
            s = Structure.from_file(args.poscar)
            oba_set = ObaSet.make_input(structure=s,
                                        user_incar_settings=user_incar_settings,
                                        **kwargs)

        oba_set.write_input(".")

    os.chdir(original_dir)


def diagnose(args):
    # try:
    #     unitcell = UnitcellCalcResults.load_json(args.unitcell)
    # except FileNotFoundError:
    #     print("{} not found".format(args.unitcell))

    # try:
    #     perfect = SupercellCalcResults.load_json(args.perfect)
    # except FileNotFoundError:
    #     print("{} not found".format(args.perfect))

    defects_dirs = args.defect_dirs if args.defect_dirs else glob('*[0-9]/')
    for d in defects_dirs:
        print(d.rjust(12), end="  ")
        dft_results = SupercellCalcResults.load_json(join(d, args.json))
        print(dft_results.diagnose)


def plot_energy(args):
    unitcell = UnitcellCalcResults.load_json(args.unitcell)
    perfect = SupercellCalcResults.load_json(args.perfect)

    defects_dirs = args.defect_dirs if args.defect_dirs else glob('*[0-9]/')

    defects = []
    for d in defects_dirs:
        logger.info("parsing directory {}...".format(d))
        input_objects = []
        files = [args.defect_entry, args.dft_results, args.correction]
        classes = [DefectEntry, SupercellCalcResults, ExtendedFnvCorrection]
        for f, c in zip(files, classes):
            try:
                input_objects.append(c.load_json(join(d, f)))
            except IOError:
                logger.warning("Parsing {} in {} failed.".format(f, d))
                continue

        defects.append(Defect(defect_entry=input_objects[0],
                              dft_results=input_objects[1],
                              correction=input_objects[2]))

    chem_pot = ChemPotDiag.load_vertices_yaml(args.chem_pot_yaml)

    # First construct DefectEnergies class object.
    defect_energies = \
        DefectEnergies.from_files(unitcell=unitcell,
                                  perfect=perfect,
                                  defects=defects,
                                  chem_pot=chem_pot,
                                  chem_pot_label=args.chem_pot_label,
                                  system=args.name)

    # if args.concentration:
    #     defect_concentration = \
    #         DefectConcentration.from_defect_energies(
    #             defect_energies=defect_energies,
    #             temperature=args.temperature[0],
    #             unitcell=unitcell,
    #             num_sites_filename=args.num_site_file)

    #     if len(args.temperature) == 2:
    #         defect_concentration = \
    #             DefectConcentration.from_defect_energies(
    #                 defect_energies=defect_energies,
    #                 temperature=args.temperature[1],
    #                 unitcell=unitcell,
    #                 num_sites_filename=args.num_site_file,
    #                 previous_concentration=defect_concentration)
    # else:
#        defect_concentration = None

    plotter = DefectEnergyPlotter(defect_energies)
    #plotter = DefectEnergyPlotter(defect_energies, defect_concentration)
    plt = plotter.plot_energy(filtering_words=args.filtering,
                              x_range=args.x_range,
                              y_range=args.y_range,
                              show_fermi_level=args.concentration,
                              show_transition_levels=args.show_tl,
                              show_all_energies=args.show_all)

    plt.savefig(args.save_file, format="pdf") if args.save_file else plt.show()


def parse_eigenvalues(args):
    unitcell = UnitcellCalcResults.load_json(args.unitcell)
    perfect = SupercellCalcResults.load_json(args.perfect)

    d = args.defect_dir
    logger.info("parsing directory {}...".format(d))
    input_objects = []
    files = [args.defect_entry, args.dft_results, args.correction]
    classes = [DefectEntry, SupercellCalcResults, ExtendedFnvCorrection]
    for f, c in zip(files, classes):
        try:
            input_objects.append(c.load_json(join(d, f)))
        except IOError:
            logger.warning("Parsing {} in {} failed.".format(f, d))
            continue

    defect = Defect(defect_entry=input_objects[0],
                    dft_results=input_objects[1],
                    correction=input_objects[2])

    defect_eigenvalues = DefectEigenvalue.from_files(unitcell=unitcell,
                                                     perfect=perfect,
                                                     defect=defect)

    defect_eigenvalues.plot()

#    plt.savefig(args.save_file, format="pdf") if args.save_file else plt.show()


def vasp_parchg_set(args):
    user_incar_settings = {"LPARD": True, "LSEPB": True, "KPAR": 1, "IBAND": args.band_indices}
    if args.kpoint_indices:
       user_incar_settings["KPUSE"] = args.kpoint_indices
    oba_set = ObaSet.from_prev_calc(dirname=args.read_dir,
                                    parse_calc_results=False,
                                    parse_magnetization=False,
                                    standardize_structure=False,
                                    sort_structure=False,
                                    parse_potcar=True,
                                    parse_incar=True,
                                    parse_kpoints=True,
                                    copied_file_names={"WAVECAR": "L"},
                                    user_incar_settings=user_incar_settings)
    oba_set.write_input(args.write_dir)


if __name__ == "__main__":
    main()
