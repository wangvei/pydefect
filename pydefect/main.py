#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.complex_defects import ComplexDefects
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.input_maker.supercell_maker import Supercells
from pydefect.main_functions import (
    recommend_supercell, initial_setting, interstitial, complex_defects,
    defect_vasp_oba_set, defect_entry, supercell_calc_results,
    unitcell_calc_results, efnv_correction, vasp_oba_set, defects, plot_energy,
    parse_eigenvalues, vasp_parchg_set, local_structure, concentration)
from pydefect.util.logger import get_logger
from pydefect.util.main_tools import get_default_args

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

    # -- recommend_supercell --------------------------------------------------
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
        default=defaults["criterion"],
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

    # -- initial_setting ------------------------------------------------------
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
        default=defaults["interstitial_sites"],
        help="Interstitial site names.")
    parser_initial.add_argument(
        "--print_dopant", dest="print_dopant", type=str,
        help="Print dopant information that can be added a posteriori.")

    parser_initial.set_defaults(func=initial_setting)

    # -- interstitial ---------------------------------------------------------
    parser_interstitial = subparsers.add_parser(
        name="interstitial",
        description="Tools for handling the interstitial sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['i'])

    defaults = get_default_args(InterstitialSiteSet.add_sites)
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
        default=defaults["vicinage_radius"],
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
    parser_interstitial.add_argument(
        "--chgcar", dest="chgcar",  type=str, default=None,
        help="CHGCAR type filename to determine the local charge minimum.")
    parser_interstitial.set_defaults(func=interstitial)

    # -- complex_defects -------------------------------------------------------
    parser_complex_defects = subparsers.add_parser(
        name="complex_defects",
        description="Tools for handling the complex defects.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cd'])

    defaults = get_default_args(ComplexDefects.add_defect)
    defaults.update(get_default_args(ComplexDefects.from_files))

    parser_complex_defects.add_argument(
        "--yaml", dest="yaml", type=str, default=defaults["filename"],
        help="Yaml-type file name to be read/written.")
    parser_complex_defects.add_argument(
        "--dposcar", dest="dposcar", type=str, default=defaults["structure"],
        help="DPOSCAR-type file name.")
    parser_complex_defects.add_argument(
        "-r", dest="removed_atom_indices", nargs="+", type=int,
        help="")
    parser_complex_defects.add_argument(
        "-i", dest="inserted_elements", nargs="+", type=str,
        help="")
    parser_complex_defects.add_argument(
        "-c", dest="inserted_coords", nargs="+", type=float,
        help="")
    parser_complex_defects.add_argument(
        "--name", dest="name", type=str,
        help="Set the complex defect name.")
    parser_complex_defects.add_argument(
        "--oxidation_state", dest="oxidation_state", type=str, default=None,
        help=".")
    parser_complex_defects.add_argument(
        "--annotation", dest="annotation", type=str, default=None,
        help=".")
    parser_complex_defects.add_argument(
        "--symprec", dest="symprec", type=float,
        default=defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_complex_defects.add_argument(
        "--angle_tol", dest="angle_tolerance", type=float,
        default=defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
    parser_complex_defects.set_defaults(func=complex_defects)

    # -- defect_vasp_set_maker ------------------------------------------------
    parser_defect_vasp_set = subparsers.add_parser(
        name="defect_vasp_oba_set",
        description="Tools for configuring vasp defect_set files for a set of "
                    "defect calculations. One needs to set "
                    ".pydefect.yaml for potcar setup.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['dvos'])

    parser_defect_vasp_set.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_defect_vasp_set.add_argument(
        "--dposcar", dest="dposcar", default="DPOSCAR", type=str,
        help="DPOSCAR-type file name.")
    parser_defect_vasp_set.add_argument(
        "-py", "--potcar_yaml", dest="potcar", default=None, type=str,
        help="User specifying Potcar yaml file.")
    parser_defect_vasp_set.add_argument(
        "-kw", "--keywords", dest="keywords", type=str, default=None,
        nargs="+", help="Filtering keywords.")
    parser_defect_vasp_set.add_argument(
        "-vos_kw", "--vos_kwargs", dest="vos_kwargs", type=str,
        default={}, nargs="+",
        help="Keywords for vasp_oba_set.")
    parser_defect_vasp_set.add_argument(
        "-d", dest="particular_defects", type=str, default=None, nargs="+",
        help="Particular defect names to be added.")
    parser_defect_vasp_set.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Set if the existing folders are overwritten.")
    parser_defect_vasp_set.add_argument(
        "-x", "--xc", dest="xc", default="pbe", type=str,
        help="XC interaction treatment.")
    parser_defect_vasp_set.add_argument(
        "-k", "--kpt_density", dest="kpt_density", default=2.3, type=float,
        help="K-points density.")
    parser_defect_vasp_set.add_argument(
        "-nw", "--no_wavecar", dest="wavecar", action="store_false",
        help="Do not make WAVECAR file or not.")
    parser_defect_vasp_set.add_argument(
        "-ldauu", dest="ldauu", type=str, default=None, nargs="+",
        help=".")
    parser_defect_vasp_set.add_argument(
        "-ldaul", dest="ldaul", type=str, default=None, nargs="+",
        help=".")
#    parser_recommend_supercell.add_argument(
#        "-set", dest="set", action="store_true",
#        help="Output all the supercells satisfying the criterion.")
    parser_defect_vasp_set.set_defaults(func=defect_vasp_oba_set)

    # -- defect_entry ---------------------------------------------------------
    parser_defect_entry = subparsers.add_parser(
        name="defect_entry",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations. By default, print defect_entry.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    defaults = get_default_args(DefectEntry.from_yaml)

    parser_defect_entry.add_argument(
        "--print", dest="print", action="store_true",
        help="Print DefectEntry class object information.")
    parser_defect_entry.add_argument(
        "--make_defect_entry", dest="make_defect_entry", action="store_true",
        help="Make defect_entry.json from yaml file.")
    parser_defect_entry.add_argument(
        "--yaml", dest="yaml", type=str, default=None,
        help="defect_entry.yaml-type file name to be read.")
    parser_defect_entry.add_argument(
        "--json", dest="json", type=str, default="defect_entry.json",
        help="defect_entry.json type file name.")
    parser_defect_entry.add_argument(
        "--cutoff", "-c", dest="cutoff", type=float,
        default=defaults["cutoff"],
        help="Cutoff radius to determine the neighboring atoms.")
    parser_defect_entry.add_argument(
        "--no_calc_sites", dest="calc_sites", action="store_false",
        help="Do not calculate number of equivalent sites (clusters).")

    parser_defect_entry.set_defaults(func=defect_entry)

    # -- supercell_calc_results -----------------------------------------------
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

    # -- unitcell_calc_results ------------------------------------------------
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

    # -- FNV correction -------------------------------------------------------
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

    # -- vasp_oba_set ---------------------------------------------------------
    parser_vasp_oba_set = subparsers.add_parser(
        name="vasp_oba_set",
        description="Tools for constructing vasp input set with oba_set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['vos'])

    # write use potcar setting
    parser_vasp_oba_set.add_argument(
        "-p", "--poscar", dest="poscar", default="POSCAR", type=str,
        help="POSCAR-type file name for the supercell.")
    parser_vasp_oba_set.add_argument(
        "-py", "--potcar_yaml", dest="potcar", default=None, type=str,
        help="User specifying Potcar yaml file.")
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
    parser_vasp_oba_set.add_argument(
        "-c", "--charge", dest="charge", type=int,
        help="Charge state.")
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
        "-vos_kw", "--vos_kwargs", dest="vos_kwargs", type=str, nargs="+",
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
    parser_vasp_oba_set.add_argument(
        "-nw", "--no_wavecar", dest="wavecar", action="store_false",
        help="Do not make WAVECAR file or not.")
    parser_vasp_oba_set.add_argument(
        "-ldauu", dest="ldauu", type=str, default=None, nargs="+",
        help=".")
    parser_vasp_oba_set.add_argument(
        "-ldaul", dest="ldaul", type=str, default=None, nargs="+",
        help=".")

    parser_vasp_oba_set.set_defaults(func=vasp_oba_set)

    # -- defects -------------------------------------------------------------
    parser_defects = subparsers.add_parser(
        name="defects",
        description="Tools for generating and analyzing Defect objects that "
                    "are used for analyzing results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['d'])

    parser_defects.add_argument(
        "--defect_dirs", dest="defect_dirs", nargs="+", type=str,
        help="Directory names for the defect supercell results. "
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in each directory.")
    parser_defects.add_argument(
        "-d", dest="diagnose", action="store_true",
        help="Show diagnosed defect results.")
    parser_defects.add_argument(
        "--json", dest="json", type=str, default="defect.json",
        help="defect.json type file name.")
    parser_defects.add_argument(
        "--perfect", dest="perfect", type=str,
        default="perfect/dft_results.json",
        help="Json file name for the SupercellCalcResults class object of the "
             "perfect supercell result.")
    parser_defects.add_argument(
        "--de", dest="defect_entry", type=str, default="defect_entry.json")
    parser_defects.add_argument(
        "--correction", dest="correction", type=str, default="correction.json")
    parser_defects.add_argument(
        "--dft_results", dest="dft_results", type=str,
        default="dft_results.json")

    parser_defects.set_defaults(func=defects)

    # -- plot_energy ----------------------------------------------------------
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
        help="File name to save the plot.")
    parser_plot_energy.add_argument(
        "--energies", dest="energies", type=str,
        default="defect_energies.json",
        help="DefectEnergies class object json file name.")
    parser_plot_energy.add_argument(
        "--unitcell", dest="unitcell", type=str, default="unitcell.json",
        help="UnitcellCalcResults class object json file name.")
    parser_plot_energy.add_argument(
        "--perfect", dest="perfect", type=str,
        default="perfect/dft_results.json",
        help="Json file name for the SupercellCalcResults class object of the "
             "perfect supercell result.")
    parser_plot_energy.add_argument(
        "--d", dest="defect", type=str, default="defect.json")
    parser_plot_energy.add_argument(
        "--defect_dirs", dest="dirs", nargs="+", type=str,
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
    # parser_plot_energy.add_argument(
    #     "-u", dest="u", action="store_true",
    #     help="Calculate the U values at given defect name and three charges")

    parser_plot_energy.set_defaults(func=plot_energy)

    # -- parse_eigenvalues ----------------------------------------------------
    parser_parse_eigenvalues = subparsers.add_parser(
        name="parse_eigenvalues",
        description="Tools for parsing defect eigenvalues",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['eig'])

    parser_parse_eigenvalues.add_argument(
        "--title", dest="title", type=str, default="",
        help="Title of the plot.")
    parser_parse_eigenvalues.add_argument(
        "-y", "--yrange", dest="y_range", type=float, nargs='+', default=None,
        help="Two float values for the y-range of the plot.")
    parser_parse_eigenvalues.add_argument(
        "-s", "--save_file", dest="save_file", type=str, default=None,
        help="File name to save the plot.")
    parser_parse_eigenvalues.add_argument(
        "--unitcell", dest="unitcell", type=str, default="unitcell.json",
        help="UnitcellCalcResults class object json file name.")
    # parser_parse_eigenvalues.add_argument(
    #     "--perfect", dest="perfect", type=str,
    #     default="perfect/dft_results.json",
    #     help="Json file name for the SupercellCalcResults class object of the "
    #          "perfect supercell result.")
    parser_parse_eigenvalues.add_argument(
        "--d", dest="defect", type=str, default="defect.json")
    parser_parse_eigenvalues.add_argument(
        "--defect_dir", dest="defect_dir", type=str,
        help="Directory name for the defect supercell result. "
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in the directory.")

    parser_parse_eigenvalues.set_defaults(func=parse_eigenvalues)

    # -- vasp_parchg_set ------------------------------------------------------
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

    # -- local_structure ------------------------------------------------------
    parser_local_structure = subparsers.add_parser(
        name="local_structure",
        description="Tools for analyzing local structure",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ls'])

    parser_local_structure.add_argument(
        "--d", dest="defect", type=str, default="defect.json")
    parser_local_structure.add_argument(
        "--show_all", dest="show_all", action="store_true")
    parser_local_structure.add_argument(
        "--defect_dirs", dest="defect_dirs", type=str, default=None, nargs="+",
        help="Directory names for the defect supercell result."
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in the directory.")

    parser_local_structure.set_defaults(func=local_structure)

    # -- concentrations -------------------------------------------------------
    parser_concentration = subparsers.add_parser(
        name="concentrations",
        description="Tools for estimating carrier and defect concentrations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['c'])

    parser_concentration.add_argument(
        "--energies", dest="energies", type=str,
        default="defect_energies.json",
        help="DefectEnergies class object json file name.")
    parser_concentration.add_argument(
        "--unitcell", dest="unitcell", type=str, default="unitcell.json",
        help="UnitcellCalcResults class object json file name.")
    parser_concentration.add_argument(
        "--filtering", dest="filtering", type=str, help="Filtering word.")
    parser_concentration.add_argument(
        "-fmto", dest="fmto", action="store_true",
        help="Set fractional magnetization to 1.")
    parser_concentration.add_argument(
        "-v", "-verbose", dest="verbose", action="store_true",
        help="Show information on estimation of concentrations.")
    parser_concentration.add_argument(
        "-t", "--temperature", dest="temperature", type=float,
        help="Temperature for calculating the Fermi level. When two "
             "temperatures are supplied, the first temperature is quenched to "
             "the second temperature.")

    parser_concentration.set_defaults(func=concentration)

    # # -- refined_structure_vasp_oba_set ---------------------------------------
    # parser_rs_vasp_oba_set = subparsers.add_parser(
    #     name="refined_structure_vasp_oba_sst",
    #     description="Tools for constructing vasp set for refined structures",
    #     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    #     aliases=['rsvos'])

    # parser_rs_vasp_oba_set.add_argument(
    #     "--defect_dirs", dest="defect_dirs", nargs="+", type=str,
    #     help="Directory names for the defect supercell results. "
    #          "defect_entry.json, dft_results.json, and correction.json files "
    #          "are required in each directory.")
    # parser_rs_vasp_oba_set.add_argument(
    #     "--dft_results", dest="dft_results", type=str,
    #     default="dft_results.json")

    # parser_rs_vasp_oba_set.set_defaults(func=make_refined_structure)

    args = parser.parse_args()
    args.func(args)

#    filtering_words=args.filtering_words,


if __name__ == "__main__":
    main()
