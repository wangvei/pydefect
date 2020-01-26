#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from typing import Union

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.cluster_defects import ClusterDefects
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.input_maker.supercell_maker import Supercells
from pydefect.core.config import DEFECT_KPT_DENSITY
from pydefect.main_functions import (
    initial_setting, interstitial, cluster_defects,
    defect_vasp_set, defect_entry, supercell_calc_results,
    unitcell_calc_results, efnv_correction, defects, plot_energy,
    parse_eigenvalues, vasp_parchg_set, local_structure, concentration)
from pydefect.util.logger import get_logger
from pydefect.util.main_tools import (
    get_default_args)
from vise.util.main_tools import get_user_settings, dict2list
from pydefect.corrections.efnv_corrections import Ewald

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

__version__ = '0.0.1'
__date__ = 'will be inserted'


def main():
    # The following keys can be set by pydefect.yaml
    setting_keys = ["symprec",
                    "defect_symprec",
                    "angle_tolerance",
                    "kpt_density",
                    "defect_kpt_density",
                    "defect_incar_setting",
                    "defect_vise_kwargs",
                    "ldauu",
                    "ldaul",
                    "xc",
                    "no_wavecar",
                    "potcar_set",
                    "outcar",
                    "contcar",
                    "vasprun",
                    "procar",
                    "vicinage_radius",
                    "displacement_distance",
                    "volume_dir",
                    "static_diele_dir",
                    "ionic_diele_dir",
                    "band_edge_dir",
                    "dos_dir",
                    "unitcell_json",
                    "perfect_json",
                    "chem_pot_yaml"]

    user_settings = get_user_settings(yaml_filename="pydefect.yaml",
                                      setting_keys=setting_keys)

    def simple_override(d: dict, overridden_keys: Union[list, str]) -> None:
        """Override dict if keys exist in user_settings.

        When the value in the user_settings is a dict, it will be changed to
        list using dict2list.
        """
        if isinstance(overridden_keys, str):
            overridden_keys = [overridden_keys]
        for key in overridden_keys:
            if key in user_settings:
                v = user_settings[key]
                if isinstance(v, dict):
                    v = dict2list(v)
                d[key] = v

    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.""",
        epilog=f"""                                 
    Author: Yu Kumagai
    Version: {__version__}                                                                 
    Last updated: {__date__}""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    # -- unitcell_calc_results ------------------------------------------------
    parser_unitcell_results = subparsers.add_parser(
        name="unitcell_results",
        description="Tools for analyzing vasp unitcell results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ur'])

    ur_defaults = {"volume_dir": None,
                   "static_diele_dir": None,
                   "ionic_diele_dir": None,
                   "band_edge_dir": None,
                   "dos_dir": None,
                   "unitcell_json": "unitcell.json",
                   "contcar": "CONTCAR",
                   "outcar": "OUTCAR",
                   "vasprun": "vasprun.xml"}

    simple_override(ur_defaults,
                    ["volume_dir",
                     "static_diele_dir",
                     "ionic_diele_dir",
                     "band_edge_dir",
                     "dos_dir",
                     "unitcell_json",
                     "contcar",
                     "outcar",
                     "vasprun"])

    parser_unitcell_results.add_argument(
        "--json_file", dest="json_file", default="unitcell.json", type=str,
        help="Json file for the unitcell info.")
    parser_unitcell_results.add_argument(
        "--static_diele", dest="static_diele", default=None, type=float,
        nargs="+",
        help="Set static (electronic) dielectric constant.")
    parser_unitcell_results.add_argument(
        "--ionic_diele", dest="ionic_diele", default=None, type=float,
        nargs="+",
        help="Set ionic dielectric constant.")

    parser_unitcell_results.add_argument(
        "--band_edge_dir", dest="band_edge_dir",
        default=ur_defaults["band_edge_dir"], type=str,
        help="Set band edge from a vasprun.xml file")

    parser_unitcell_results.add_argument(
        "--static_diele_dir", dest="static_diele_dir",
        default=ur_defaults["static_diele_dir"], type=str,
        help="Set static dielectric constant from an OUTCAR file")
    parser_unitcell_results.add_argument(
        "--ionic_diele_dir", dest="ionic_diele_dir",
        default=ur_defaults["ionic_diele_dir"], type=str,
        help="Set ionic dielectric constant from an OUTCAR file")

    parser_unitcell_results.add_argument(
        "--volume_dir", dest="volume_dir",
        default=ur_defaults["volume_dir"], type=str,
        help="Set volume from a POSCAR file")

    parser_unitcell_results.add_argument(
        "--total_dos_dir", dest="total_dos_dir",
        default=ur_defaults["dos_dir"], type=str,
        help="Set total density of states from a vasprun.xml file")
    parser_unitcell_results.add_argument(
        "-c", dest="contcar", type=str, default=ur_defaults["contcar"])
    parser_unitcell_results.add_argument(
        "-o", dest="outcar", type=str, default=ur_defaults["outcar"])
    parser_unitcell_results.add_argument(
        "-v", dest="vasprun", type=str, default=ur_defaults["vasprun"])
    parser_unitcell_results.add_argument(
        "--print", dest="print", action="store_true",
        help="Print Unitcell class object information.")

    del ur_defaults

    parser_unitcell_results.set_defaults(func=unitcell_calc_results)

    # -- initial_setting ------------------------------------------------------
    parser_initial = subparsers.add_parser(
        name="initial_setting",
        description="Tools for recommending an optimal supercell(s) and "
                    "configuring initial settings for a standard set "
                    "of point-defect calculations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['is'])

    is_defaults = get_default_args(DefectInitialSetting.from_basic_settings)
    is_defaults.update(get_default_args(Supercells))
    simple_override(is_defaults,
                    ["symprec", "angle_tolerance", "displacement_distance"])

    parser_initial.add_argument(
        "--symprec", dest="symprec", type=float,
        default=is_defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")

    # supercell
    parser_initial.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR",
        help="Input poscar file name.")
    parser_initial.add_argument(
        "--matrix", dest="matrix", type=int, nargs="+", default=None,
        help="Generate the supercell based on the Transformation matrix."
             "1, 3, or 9 components are accepted..")
    parser_initial.add_argument(
        "-c", "--criterion", dest="isotropy_criterion", type=float,
        default=is_defaults["criterion"],
        help="Criterion used for screening candidate supercells.")
    parser_initial.add_argument(
        "--min_num_atoms", dest="min_num_atoms", type=int,
        default=is_defaults["min_num_atoms"],
        help="Minimum number of atoms in the candidate supercells")
    parser_initial.add_argument(
        "--max_num_atoms", dest="max_num_atoms", type=int,
        default=is_defaults["max_num_atoms"],
        help="Maximum number of atoms in the candidate supercells")
    parser_initial.add_argument(
        "-pr", "--primitive", dest="primitive", action="store_true",
        help="Set when the supercell is expanded based on the primitive cell."
             "When the conventional and primitive unit cells are the same,"
             "this flag has no meaning.")
    parser_initial.add_argument(
        "-i", "--most_isotropic", dest="most_isotropic", action="store_true",
        help="Generate the supercell with the smallest isotropy instead of the "
             "smallest supercell.")

    # defect.in
    parser_initial.add_argument(
        "--rhombohedral_angle", dest="rhombohedral_angle", type=float,
        default=is_defaults["rhombohedral_angle"],
        help="Only the supercells with rhombohedral_angle <= lattice angle <= "
             "180 - rhombohedral_angle are returned. ")
    parser_initial.add_argument(
        "-set", dest="set", action="store_true",
        help="Output all the supercells satisfying the criterion.")
    parser_initial.add_argument(
        "-d", "--dopants", dest="dopants", nargs="+", type=str,
        default=is_defaults["dopants"],
        help="Dopant elements, e.g., Ga In")
    parser_initial.add_argument(
        "-a", "--antisite", dest="is_antisite", action="store_false",
        help="Set if antisite defects are *not* considered.")
    parser_initial.add_argument(
        "--en_diff", dest="en_diff", type=float, default=is_defaults["en_diff"],
        help="Criterion of the electronegativity difference that determines "
             "substituted impurity types.")
    parser_initial.add_argument(
        "--included", dest="included", type=str, nargs="+",
        default=is_defaults["included"],
        help="Exceptionally included defects with full names. E.g., Va_O2_-1.")
    parser_initial.add_argument(
        "--excluded", dest="excluded", type=str, nargs="+",
        default=is_defaults["excluded"],
        help="Exceptionally excluded defects with full names. E.g., Va_O2_0.")
    parser_initial.add_argument(
        "--displacement_distance", dest="displacement_distance", type=float,
        default=is_defaults["displacement_distance"],
        help="Displacement distance. 0 means that random displacement is not "
             "considered.")
    parser_initial.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=is_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
    parser_initial.add_argument(
        "--interstitial_sites", dest="interstitials", type=str, nargs="+",
        default=is_defaults["interstitial_sites"],
        help="Interstitial site names.")
    parser_initial.add_argument(
        "--cluster_defect_names", dest="cluster_defect_names", type=str,
        nargs="+", default=is_defaults["cluster_defect_names"],
        help="Complex defect names.")
    parser_initial.add_argument(
        "--print_dopant", dest="print_dopant", type=str,
        help="Print dopant information that can be added a posteriori.")

    del is_defaults

    parser_initial.set_defaults(func=initial_setting)

    # -- interstitial ---------------------------------------------------------
    parser_interstitial = subparsers.add_parser(
        name="interstitial",
        description="Tools for handling the interstitial sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['i'])

    i_defaults = get_default_args(InterstitialSiteSet.add_sites)
    i_defaults.update(get_default_args(InterstitialSiteSet.from_files))

    simple_override(i_defaults,
                    ["vicinage_radius", "symprec", "angle_tolerance"])

    parser_interstitial.add_argument(
        "--yaml", dest="yaml", type=str, default=i_defaults["yaml_filename"],
        help="interstitial.yaml file name.")
    parser_interstitial.add_argument(
        "--dposcar", dest="dposcar", type=str, default=i_defaults["dposcar"],
        help="DPOSCAR-type file name.")
    parser_interstitial.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_interstitial.add_argument(
        "-c", dest="interstitial_coords", nargs="+", type=float,
        help="Interstitial coordinates in the UPOSCAR cell. Eg., 0.5 0.5 0.5.")
    parser_interstitial.add_argument(
        "--name", dest="site_name", type=str, default=None,
        help="Interstitial site name.")
    parser_interstitial.add_argument(
        "--radius", dest="radius", type=str,
        default=i_defaults["vicinage_radius"],
        help="Radius to judge whether too close atoms exist.")
    parser_interstitial.add_argument(
        "--force_add", dest="force_add", action="store_true",
        help="Set if the interstitial site is forcibly added, although it is "
             "too close to other atoms.")
    parser_interstitial.add_argument(
        "--symprec", dest="symprec", type=float,
        default=i_defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_interstitial.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=i_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
    parser_interstitial.add_argument(
        "--method", dest="method", type=str, default=i_defaults["method"],
        help="Name of method determining the interstitial site.")
    parser_interstitial.add_argument(
        "--chgcar", dest="chgcar",  type=str, default=None,
        help="CHGCAR type filename to determine the local charge minima.")

    del i_defaults

    parser_interstitial.set_defaults(func=interstitial)

    # -- cluster_defects -------------------------------------------------------
    parser_cluster_defects = subparsers.add_parser(
        name="cluster_defects",
        description="Tools for handling the cluster defects.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cd'])

    cd_defaults = get_default_args(ClusterDefects.add_defect)
    cd_defaults.update(get_default_args(ClusterDefects.from_files))
    simple_override(cd_defaults, ["symprec", "angle_tolerance"])

    parser_cluster_defects.add_argument(
        "--yaml", dest="yaml", type=str, default=cd_defaults["yaml_filename"],
        help="cluster_defect.yaml file name.")
    parser_cluster_defects.add_argument(
        "--dposcar", dest="dposcar", type=str, default=cd_defaults["structure"],
        help="DPOSCAR-type file name.")
    parser_cluster_defects.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_cluster_defects.add_argument(
        "-r", dest="removed_atom_indices", nargs="+", type=int,
        help="Removed atom indices (beginning at 0) from the pristine supercell"
             " when constructing a cluster defect.")
    parser_cluster_defects.add_argument(
        "-i", dest="inserted_elements", nargs="+", type=str,
        help="Inserted atom elements when constructing a cluster defect."
             "E.g., Mg Al")
    parser_cluster_defects.add_argument(
        "-c", dest="inserted_coords", nargs="+", type=float,
        help="Inserted atom coordinates when constructing a cluster defect."
             "E.g., 0 0 0  0.25 0.25 0.25")
    parser_cluster_defects.add_argument(
        "--name", dest="name", type=str, help="Complex defect name.")
    parser_cluster_defects.add_argument(
        "--extreme_charge_state", dest="extreme_charge_state", type=int,
        help="Extreme charge state of the cluster defect. For example, in case"
             "of V_Na + V_O, the extreme charge state would be "
             "(-1) + (+2) = 1.")
    parser_cluster_defects.add_argument(
        "--annotation", dest="annotation", type=str, default=None,
        help="Annotation of the cluster defect.")
    parser_cluster_defects.add_argument(
        "--symprec", dest="symprec", type=float,
        default=cd_defaults["symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_cluster_defects.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=cd_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")

    del cd_defaults

    parser_cluster_defects.set_defaults(func=cluster_defects)

    # -- defect_vasp_set_maker ------------------------------------------------
    parser_defect_vasp_set = subparsers.add_parser(
        name="defect_vasp_set",
        description="Tools for configuring vasp defect_set files for a set of "
                    "defect calculations. One needs to set "
                    ".pydefect.yaml for potcar setup.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['dvs'])

    # all the defaults must be declared here.
    dvs_defaults = {"dvs_kwargs": dict(),
                    "xc":         "pbesol",
                    "defect_kpt_density": DEFECT_KPT_DENSITY,
                    "defect_incar_setting": None,
                    "potcar_set": None,
                    "ldauu":      None,
                    "ldaul":      None}

    simple_override(dvs_defaults,
                    ["xc",
                     "defect_kpt_density",
                     "defect_incar_setting",
                     "potcar_set",
                     "ldauu",
                     "ldaul"])

    dvs_defaults["dvs_kwargs"].update(
        user_settings.get("defect_vise_kwargs", {}))
    simple_override(
        dvs_defaults["dvs_kwargs"], ["symprec", "angle_tolerance"])
    dvs_defaults["dvs_kwargs"] = dict2list(dvs_defaults["dvs_kwargs"])

    parser_defect_vasp_set.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_defect_vasp_set.add_argument(
        "--dposcar", dest="dposcar", default="DPOSCAR", type=str,
        help="DPOSCAR-type file name.")
    parser_defect_vasp_set.add_argument(
        "--potcar", dest="potcar_set", default=dvs_defaults["potcar_set"],
        type=str, nargs="+",
        help="User specifying POTCAR set. E.g., Mg_pv O_h")
    parser_defect_vasp_set.add_argument(
        "-x", "--xc", dest="xc", default=dvs_defaults["xc"], type=str,
        help="XC interaction treatment.")
    parser_defect_vasp_set.add_argument(
        "-kw", "--keywords", dest="keywords", type=str, default=None,
        nargs="+", help="Filtering keywords.")
    parser_defect_vasp_set.add_argument(
        "-vos_kw", "--vos_kwargs", dest="vos_kwargs", type=str,
        default=dvs_defaults["dvs_kwargs"], nargs="+",
        help="Keywords for vasp_oba_set.")
    parser_defect_vasp_set.add_argument(
        "-is", "--incar_setting", dest="incar_setting", type=str, nargs="+",
        default=dvs_defaults["defect_incar_setting"],
        help="user_incar_setting in make_input classmethod of ObaSet in vise. "
             "See document in vise for details.")
    parser_defect_vasp_set.add_argument(
        "-d", dest="specified_defects", type=str, default=None, nargs="+",
        help="Particularly specified defect names to be added.")
    parser_defect_vasp_set.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Set if the existing folders are overwritten.")
    parser_defect_vasp_set.add_argument(
        "-k", "--kpt_density", dest="kpt_density",
        default=dvs_defaults["defect_kpt_density"], type=float,
        help="K-point density in Angstrom along each direction .")
    parser_defect_vasp_set.add_argument(
        "-nw", "--no_wavecar", dest="no_wavecar", action="store_true",
        help="Do not make WAVECAR file or not.")
    parser_defect_vasp_set.add_argument(
        "-ldauu", dest="ldauu", type=dict, default=dvs_defaults["ldauu"],
        nargs="+", help=".")
    parser_defect_vasp_set.add_argument(
        "-ldaul", dest="ldaul", type=str, default=dvs_defaults["ldaul"],
        nargs="+", help=".")

    del dvs_defaults

    parser_defect_vasp_set.set_defaults(func=defect_vasp_set)

    # -- defect_entry ---------------------------------------------------------
    parser_defect_entry = subparsers.add_parser(
        name="defect_entry",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    de_defaults = get_default_args(DefectEntry.from_defect_structure)
    simple_override(de_defaults, ["displacement_distance"])

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
        "--perfect_poscar", dest="perfect_poscar", type=str,
        default="../perfect/POSCAR",
        help="Perfect supercell POSCAR filename with a path from the current"
             "working directory.")
    parser_defect_entry.add_argument(
        "--defect_poscar", dest="defect_poscar", type=str,
        default="POSCAR",
        help="Defect supercell POSCAR filename.")
    parser_defect_entry.add_argument(
        "--displacement_distance", dest="displacement_distance", type=float,
        default=de_defaults["displacement_distance"],
        help="Displacement distance.")
    parser_defect_entry.add_argument(
        "--defect_name", dest="defect_name", type=str,
        default=de_defaults["defect_name"],
        help="Defect name.")

    del de_defaults

    parser_defect_entry.set_defaults(func=defect_entry)

    # -- supercell_calc_results -----------------------------------------------
    parser_supercell_results = subparsers.add_parser(
        name="supercell_results",
        description="Tools for analyzing vasp supercell results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sr'])

    sr_defaults = get_default_args(SupercellCalcResults.from_vasp_files)
    simple_override(sr_defaults,
                    ["vasprun",
                     "contcar",
                     "outcar",
                     "procar",
                     "cutoff",
                     "defect_symprec",
                     "angle_tolerance"])

    parser_supercell_results.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Directory names.")
    parser_supercell_results.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and " "perfect directory.")
    parser_supercell_results.add_argument(
        "-v", dest="vasprun", default=sr_defaults["vasprun"], type=str)
    parser_supercell_results.add_argument(
        "-c", dest="contcar", default=sr_defaults["contcar"], type=str)
    parser_supercell_results.add_argument(
        "-o", dest="outcar", default=sr_defaults["outcar"], type=str)
    parser_supercell_results.add_argument(
        "-p", dest="procar", default=sr_defaults["procar"], type=str)
    parser_supercell_results.add_argument(
        "-s", dest="not_check_shallow", action="store_false",
        help="Not check whether the defects are shallow.")
    parser_supercell_results.add_argument(
        "--center", dest="defect_center", nargs="+", type=float,
        help="Set defect center in fractional coordinates or atomic index.")
    parser_supercell_results.add_argument(
        "-de", dest="defect_entry_name", type=str, default="defect_entry.json")
    parser_supercell_results.add_argument(
        "--json", dest="json", type=str, default="dft_results.json",
        help="dft_results.json type file name.")
    parser_supercell_results.add_argument(
        "--cutoff", dest="cutoff", type=float, default=sr_defaults["cutoff"],
        help="Cutoff radius to determine the neighboring atoms.")
    parser_supercell_results.add_argument(
        "--defect_symprec", dest="symprec", type=float,
        default=sr_defaults["defect_symprec"],
        help="Set length precision used for symmetry analysis [A].")
    parser_supercell_results.add_argument(
        "--angle_toleranec", dest="angle_tolerance", type=float,
        default=sr_defaults["angle_tolerance"],
        help="Set angle precision used for symmetry analysis.")
    parser_supercell_results.add_argument(
        "--print", dest="print", action="store_true",
        help="Print SupercellCalcResults class object information.")

    del sr_defaults

    parser_supercell_results.set_defaults(func=supercell_calc_results)

    # -- extended FNV correction -----------------------------------------------
    parser_correction = subparsers.add_parser(
        name="extended_fnv_correction",
        description="Tools for extended FNV correction for error of defect "
                    "formation energy due to finite supercell size.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['efc'])

    efc_defaults = get_default_args(Ewald.from_optimization)
    efc_defaults.update({"unitcell_json": "../unitcell/unitcell.json",
                         "perfect_json": "perfect/dft_results.json"})
    simple_override(efc_defaults, ["unitcell_json", "perfect_json"])

    # when no corrections are required for some reasons
    parser_correction.add_argument(
        "--nocorr", "-nc", dest="nocorr", action="store_true",
        help="Set when no corrections are required for some reasons.")
    parser_correction.add_argument(
        "-m", "--manual", dest="manual", type=float, default=0.0,
        help="Set manual correction energy.")

    # needed files
    parser_correction.add_argument(
        "--unitcell_json", dest="unitcell_json",
        default=efc_defaults["unitcell_json"],
        type=str)
    parser_correction.add_argument(
        "--perfect_json", dest="perfect_json",
        default=efc_defaults["perfect_json"], type=str)

    # ewald
    parser_correction.add_argument(
        "--ewald_json", dest="ewald_json", default="ewald.json",
        type=str)
    parser_correction.add_argument(
        "--ewald_initial_param", dest="ewald_initial_param",
        default=efc_defaults["initial_ewald_param"])
    parser_correction.add_argument(
        "--ewald_convergence", dest="ewald_convergence",
        default=efc_defaults["convergence"])
    parser_correction.add_argument(
        "--ewald_accuracy", dest="ewald_accuracy",
        default=efc_defaults["prod_cutoff_fwhm"])

    # correction
    parser_correction.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str, help="Directory names.")
    parser_correction.add_argument(
        "--dir_all", dest="dir_all", action="store_true",
        help="Make dft_results.json for *[0-9] and perfect directory.")
    parser_correction.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Overwrite already existing correction.json.")
    parser_correction.add_argument(
        "--center", dest="defect_center", nargs="+", type=float, default=None,
        help="Set defect center in fractional coordinates.")
    parser_correction.add_argument(
        "--json_file", dest="json_file", default="correction.json", type=str,
        help="Json file for the correction.")
    parser_correction.add_argument(
        "--print", dest="print", action="store_true",
        help="Print FNV correction information.")

    # plot potential
    parser_correction.add_argument(
        "--plot_potential", dest="plot_potential", action="store_true",
        help="Only replot the potential profile.")
    parser_correction.add_argument(
        "-y", "--y_range", dest="y_range", type=float, nargs='+', default=None,
        help="Two float values for the y-range of the plot.")

    del efc_defaults

    parser_correction.set_defaults(func=efnv_correction)

    # -- defects -------------------------------------------------------------
    parser_defects = subparsers.add_parser(
        name="defects",
        description="Tools for generating and analyzing Defect objects that "
                    "are used for analyzing results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['d'])

    d_defaults = {"perfect_json": "perfect/dft_results.json"}
    simple_override(d_defaults, "perfect_json")

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
        default=d_defaults["perfect_json"],
        help="Json file name for the SupercellCalcResults class object of the "
             "perfect supercell result.")
    parser_defects.add_argument(
        "--de", dest="defect_entry", type=str, default="defect_entry.json")
    parser_defects.add_argument(
        "--correction", dest="correction", type=str, default="correction.json")
    parser_defects.add_argument(
        "--dft_results", dest="dft_results", type=str,
        default="dft_results.json")
    parser_defects.add_argument(
        "-be", dest="band_edge", type=str, nargs="+", default=None)

    del d_defaults

    parser_defects.set_defaults(func=defects)

    # -- plot_energy ----------------------------------------------------------
    parser_plot_energy = subparsers.add_parser(
        name="plot_energy",
        description="Tools for plotting defect formation energies as a "
                    "function of Fermi level",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pe'])

    pe_defaults = {"unitcell_json": "../unitcell/unitcell.json",
                   "perfect_json": "perfect/dft_results.json",
                   "chem_pot_yaml": "../competing_phases/vertices.yaml"}
    simple_override(pe_defaults,
                    ["unitcell_json", "perfect_json", "chem_pot_yaml"])

    parser_plot_energy.add_argument(
        "--name", dest="name", type=str, default="",
        help="System name that is written in the title.")
    parser_plot_energy.add_argument(
        "-x", "--x_range", dest="x_range", type=float, nargs='+', default=None,
        help="Two float values for the x-range of the plot wrt the VBM.")
    parser_plot_energy.add_argument(
        "-y", "--y_range", dest="y_range", type=float, nargs='+', default=None,
        help="Two float values for the y-range of the plot.")
    parser_plot_energy.add_argument(
        "-s", "--save_file", dest="save_file", type=str, default=None,
        help="File name to save the plot.")
    parser_plot_energy.add_argument(
        "--energies", dest="energies", type=str,
        default="defect_energies.json",
        help="DefectEnergies class object json file name.")
    parser_plot_energy.add_argument(
        "--unitcell", dest="unitcell", type=str,
        default=pe_defaults["unitcell_json"],
        help="UnitcellCalcResults class object json file name.")
    parser_plot_energy.add_argument(
        "--perfect", dest="perfect", type=str,
        default=pe_defaults["perfect_json"],
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
        default=pe_defaults["chem_pot_yaml"],
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

    del pe_defaults

    parser_plot_energy.set_defaults(func=plot_energy)

    # -- parse_eigenvalues ----------------------------------------------------
    parser_parse_eigenvalues = subparsers.add_parser(
        name="parse_eigenvalues",
        description="Tools for parsing defect eigenvalues",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['eig'])

    eig_defaults = {"unitcell_json": "../unitcell/unitcell.json"}
    simple_override(eig_defaults, "unitcell_json")

    parser_parse_eigenvalues.add_argument(
        "--title", dest="title", type=str, default="",
        help="Title of the plot.")
    parser_parse_eigenvalues.add_argument(
        "-y", "--y_range", dest="y_range", type=float, nargs='+', default=None,
        help="Two float values for the y-range of the plot.")
    parser_parse_eigenvalues.add_argument(
        "-s", "--save_file", dest="save_file", type=str, default=None,
        help="File name to save the plot.")
    parser_parse_eigenvalues.add_argument(
        "--unitcell", dest="unitcell", type=str,
        default=eig_defaults["unitcell_json"],
        help="UnitcellCalcResults class object json file name.")
    parser_parse_eigenvalues.add_argument(
        "--d", dest="defect", type=str, default="defect.json")
    parser_parse_eigenvalues.add_argument(
        "--defect_dir", dest="defect_dir", type=str,
        help="Directory name for the defect supercell result. "
             "defect_entry.json, dft_results.json, and correction.json files "
             "are required in the directory.")

    del eig_defaults

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
        "--cs", dest="compare_structure", action="store_true",
        help="Compare the structures between different charge states.")
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

    c_defaults = {"unitcell_json": "../unitcell/unitcell.json"}
    simple_override(c_defaults, "unitcell_json")

    parser_concentration.add_argument(
        "--energies", dest="energies", type=str,
        default="defect_energies.json",
        help="DefectEnergies class object json file name.")
    parser_concentration.add_argument(
        "--unitcell", dest="unitcell", type=str,
        default=c_defaults["unitcell_json"],
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

    del c_defaults

    parser_concentration.set_defaults(func=concentration)


    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

