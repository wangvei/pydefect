#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from glob import glob
import os
import warnings

from pymatgen.core.structure import Structure

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.input_maker.defect_initial_setting \
    import print_dopant_info, DefectInitialSetting
from pydefect.input_maker.vasp_input_maker \
    import make_hpkot_primitive_poscar, make_supercell_poscar, make_incar, \
    make_kpoints, make_potcar
from pydefect.input_maker.vasp_defect_set_maker import VaspDefectInputSetMaker
from pydefect.input_maker.recommend_supercell_ase_cythonized.\
    ase_generation_supercell \
    import recommend_supercell_ase
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.util.vasp_process_analyzer \
    import check_vasp_output, vasp_convergence_ionic, \
    vasp_convergence_electronic
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.core.correction import Ewald, Correction
from pydefect.analysis.chempotdiag.chem_pot_diag \
    import ChemPotDiag
from pydefect.analysis.chempotdiag.make_inputs import make_vasp_inputs_from_mp

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
_ANGLE_TOLERANCE = -1.0


def main():
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package for first-principles point defect calculations. It 
    allows us to construct input files, parse first-principles calculation 
    results, and analyze data.""",
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
        aliases=['is'])

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

    # -- vasp_defect_set_maker -------------------------------------------------
    parser_vasp_defect_set = subparsers.add_parser(
        name="vasp_defect_set",
        description="Tools for configuring vasp defect_set files for a set of "
                    "defect calculations. One needs to set "
                    ".pydefect.yaml for potcar setup.",
        aliases=['vds'])

    parser_vasp_defect_set.add_argument(
        "--defect_in", dest="defect_in", default="defect.in", type=str,
        help="defect.in-type file name.")
    parser_vasp_defect_set.add_argument(
        "--dposcar", dest="dposcar", default="DPOSCAR", type=str,
        help="DPOSCAR-type file name.")
    parser_vasp_defect_set.add_argument(
        "--incar", dest="incar", default="INCAR", type=str,
        help="INCAR-type file name.")
    parser_vasp_defect_set.add_argument(
        "--kpoints", dest="kpoints", default="KPOINTS", type=str,
        help="KPOINTS-type file name.")
    parser_vasp_defect_set.add_argument(
        "--filtering", dest="filtering", type=str, default=None, nargs="+",
        help="Filtering kwargs.")
    parser_vasp_defect_set.add_argument(
        "--add", dest="add", type=str, default=None, nargs="+",
        help="Particular defect names to be added.")
    parser_vasp_defect_set.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Set if the existing folders are overwritten.")
    parser_vasp_defect_set.add_argument(
        "--make_incar", dest="make_incar", action="store_true",
        help="Make INCAR file using several default setting.")
    parser_vasp_defect_set.add_argument(
        "--make_kpoints", dest="make_kpoints", action="store_true",
        help="Make KPOINTS file based on the lattice constants.")

    parser_vasp_defect_set.set_defaults(func=vasp_defect_set)

    # -- vasp_kpoints_maker ----------------------------------------------------
    parser_vasp_kpoints_maker = subparsers.add_parser(
        name="vasp_kpoints_maker",
        description="Tools for configuring vasp KPOINTS file depending on the "
                    "task",
        aliases=['vkm'])

    parser_vasp_kpoints_maker.add_argument(
        "--task", "-t", dest="task", type=str, required=True,
        choices=["structure_opt", "band", "dos", "dielectric",
                 "dielectric_function", "competing_phase",
                 "competing_phase_molecule", "defect"],
        help="Task.")
    parser_vasp_kpoints_maker.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_vasp_kpoints_maker.add_argument(
        "-pp", dest="pposcar", type=str, default="PPOSCAR")
    parser_vasp_kpoints_maker.add_argument(
        "--num_split_kpoints", dest="num_split_kpoints", type=int, default=1,
        help="Number of the divided KPOINTS.")
    parser_vasp_kpoints_maker.add_argument(
        "--is_metal", dest="is_metal", action="store_true", default=False,
        help="Set if the system metal is metal, for which k-point density is "
             "increased.")
    parser_vasp_kpoints_maker.add_argument(
        "--kpts_shift", dest="kpts_shift", default=None, nargs="+", type=int,
        help="Origin of the k-points.")
    parser_vasp_kpoints_maker.add_argument(
        "--kpts_density_opt", dest="kpts_density_opt", type=float, default=3,
        help="K-point density used for the structure optimization of systems "
             "with band gaps ")
    parser_vasp_kpoints_maker.add_argument(
        "--kpts_density_defect", dest="kpts_density_defect", type=float,
        default=1.5,
        help="K-point density used for the calculations of point defects.")
    parser_vasp_kpoints_maker.add_argument(
        "--multiplier_factor", dest="multiplier_factor", type=float, default=2,
        help="Multiplier_factor for the calculations of density of states, "
             "dielectric constants, and dielectric function.")
    parser_vasp_kpoints_maker.add_argument(
        "--multiplier_factor_metal", dest="multiplier_factor_metal", type=float,
        default=2,
        help="Multiplier factor the structure optimization of metallic systems")

    parser_vasp_kpoints_maker.set_defaults(func=vasp_kpoints_maker)

    # -- vasp_incar_maker ----------------------------------------------------
    parser_vasp_incar_maker = subparsers.add_parser(
        name="vasp_incar_maker",
        description="Tools for configuring vasp INCAR file depending on the "
                    "task",
        aliases=['vim'])

    parser_vasp_incar_maker.add_argument(
        "--task", "-t", dest="task", type=str, required=True,
        choices=["structure_opt", "band", "dos", "dielectric",
                 "dielectric_function", "competing_phase",
                 "competing_phase_molecule", "defect"],
        help="Task.")
    parser_vasp_incar_maker.add_argument(
        "--functional", "-f", dest="functional", type=str, required=True,
        choices=["pbe", "hse06", "pbesol", "pbe_d3"],
        help="Functional.")
    parser_vasp_incar_maker.add_argument(
        "--hfscreen", dest="hfscreen", type=float,
        help="Screening distance for exchange interaction.")
    parser_vasp_incar_maker.add_argument(
        "--aexx", dest="aexx", type=float,
        help="Mixing parameter for exchange interaction.")
    parser_vasp_incar_maker.add_argument(
        "--is_magnetization", dest="is_magnetization", action="store_true",
        help="Set if the system metal is spin polarized.")
    parser_vasp_incar_maker.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")

    parser_vasp_incar_maker.set_defaults(func=vasp_incar_maker)

    # -- vasp_poscar_maker ----------------------------------------------------
    parser_vasp_poscar_maker = subparsers.add_parser(
        name="vasp_poscar_maker",
        description="Tools for configuring vasp POSCAR file. By default, "
                    "standardized primitive cell is generated.",
        aliases=['vpsm'])

    parser_vasp_poscar_maker.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR")
    parser_vasp_poscar_maker.add_argument(
        "-pp", dest="pposcar", type=str, default="PPOSCAR")
    parser_vasp_poscar_maker.add_argument(
        "--supercell", "-s", dest="supercell", type=float, nargs="+",
        help="Construct a supercell.")
    parser_vasp_poscar_maker.add_argument(
        "--symprec", dest="symprec", type=float, default=_SYMPREC,
        help="Set precision used for symmetry analysis [A].")
    parser_vasp_poscar_maker.add_argument(
        "--angle_tolerance", dest="angle_tolerance", type=float,
        default=_ANGLE_TOLERANCE,
        help="Set the angle_tolerance used for symmetry analysis.")

    parser_vasp_poscar_maker.set_defaults(func=vasp_poscar_maker)

    # -- vasp_potcar_maker ----------------------------------------------------
    parser_vasp_potcar_maker = subparsers.add_parser(
        name="vasp_potcar_maker",
        description="Tools for configuring vasp POTCAR file.",
        aliases=['vptm'])

    parser_vasp_potcar_maker.add_argument(
        "--path", dest="path", type=str, default=".")

    parser_vasp_potcar_maker.add_argument(
        "--elements", "-e", dest="elements", type=str, nargs="+",
        help="Element names.")
    parser_vasp_potcar_maker.add_argument(
        "-p", dest="poscar", type=str, default="POSCAR",
        help="Element names are also obtained from a POSCAR file.")

    parser_vasp_potcar_maker.set_defaults(func=vasp_potcar_maker)

    # -- recommend_supercell ---------------------------------------------------
    parser_recommend_supercell = subparsers.add_parser(
        name="recommend_supercell",
        description="Tools for recommendation of optimal supercell",
        aliases=['rs'])
    parser_recommend_supercell.add_argument(
        "--poscar_path", dest="poscar_path", type=str,
        help="Path of poscar to make supercell")
    parser_recommend_supercell.add_argument(
        "--criterion", dest="criterion", type=float,
        help="Criterion of ASE optimality measure (default is 0.5)")
    parser_recommend_supercell.add_argument(
        "--min_natom", dest="min_natom", type=int,
        help="Minimum number of atom (default is 50)")
    parser_recommend_supercell.add_argument(
        "--max_natom", dest="max_natom", type=int,
        help="Maximum number of atom (default is 400)")

    parser_recommend_supercell.set_defaults(func=recommend_supercell)

    # -- defect_entry ----------------------------------------------------------
    parser_defect_entry = subparsers.add_parser(
        name="defect_entry",
        description="Tools for configuring defect_entry files for post process"
                    "of defect calculations.",
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

    # -- supercell_dft_results -------------------------------------------------
    parser_supercell_results = subparsers.add_parser(
        name="supercell_results",
        description="Tools for analyzing vasp supercell results",
        aliases=['sr'])

    parser_supercell_results.add_argument(
        "-c", "--convergence", dest="convergence", action="store_true",
        help="Check convergence of vasp calc at electronic and ionic steps.")
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
    parser_supercell_results.add_argument(
        "--json", dest="json", type=str, default="dft_results.json",
        help="dft_results.json type file name.")
    parser_supercell_results.add_argument(
        "--print", dest="print", action="store_true",
        help="Print SupercellDftResults class object information.")

    parser_supercell_results.set_defaults(func=supercell_results)

    # -- unitcell_dft_results -------------------------------------------------
    parser_unitcell_results = subparsers.add_parser(
        name="unitcell_results",
        description="Tools for analyzing vasp unitcell results",
        aliases=['ur'])

    parser_unitcell_results.add_argument(
        "--json_file", dest="json_file", default="unitcell.json", type=str)

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

    parser_unitcell_results.set_defaults(func=unitcell_results)

    # -- correction ------------------------------------------------------------
    parser_correction = subparsers.add_parser(
        name="correction",
        description="Tools for correction of error of defect formation energy"
                    " due to finite cell size.",
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
    parser_correction.add_argument(
        "--force_overwrite", dest="force_overwrite", action="store_true",
        help="Overwrite already existing correction.json.")

    parser_correction.set_defaults(func=correction)

    # -- chempotdiag -----------------------------------------------------------
    parser_chempotdiag = subparsers.add_parser(
        name="chempotdiag",
        description="",
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
                                    help="")
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
        aliases=['pe'])

    parser_plot_energy.add_argument("--name", dest="name", type=str, default="")
#    parser_plot_energy.add_argument("--xrange", dest="xrange", type=str,
#                                    nargs='+', default=None,
#                                    help="X range for the plot.")
#    parser_plot_energy.add_argument("--yrange", dest="yrange", type=str,
#                                    nargs='+', default=None,
#                                    help="Y range for the plot.")
    parser_plot_energy.add_argument("-s", "--save_file", dest="save_file",
                                    type=str, default=None,
                                    help="File name to save the drawn plot.")
    parser_plot_energy.add_argument("--unitcell", dest="unitcell", type=str,
                                    default="unitcell.json")
    parser_plot_energy.add_argument("--perfect", dest="perfect", type=str,
                                    default="perfect/dft_results.json")
    parser_plot_energy.add_argument("--chem_pot_yaml", dest="chem_pot_yaml",
                                    type=str,
                                    default="chem_pot.yaml")
    parser_plot_energy.add_argument("--chem_pot_label", dest="chem_pot_label",
                                    type=str,
                                    default="A")
    parser_plot_energy.add_argument("--defect_dirs", dest="defect_dirs",
                                    default="", type=str)
    parser_plot_energy.add_argument("--show_tls", dest="--show_tls",
                                    action="store_true",
                                    help="Show the transition levels.")

    parser_plot_energy.set_defaults(func=plot_energy)

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


def vasp_defect_set(args):
    defect_initial_setting = DefectInitialSetting. \
        from_defect_in(poscar=args.dposcar, defect_in_file=args.defect_in)

    VaspDefectInputSetMaker(defect_initial_setting=defect_initial_setting,
                            filtering_words=args.filtering,
                            particular_defects=args.add,
                            incar=args.incar,
                            kpoints=args.kpoints,
                            force_overwrite=args.force_overwrite)


def vasp_kpoints_maker(args):

    make_kpoints(task=args.task,
                 poscar=args.poscar,
                 num_split_kpoints=args.num_split_kpoints,
                 is_metal=args.is_metal,
                 kpts_shift=args.kpts_shift,
                 kpts_density_opt=args.kpts_density_opt,
                 kpts_density_defect=args.kpts_density_defect,
                 multiplier_factor=args.multiplier_factor,
                 multiplier_factor_metal=args.multiplier_factor_metal)


def vasp_incar_maker(args):

    make_incar(task=args.task,
               functional=args.functional,
               hfscreen=args.hfscreen,
               aexx=args.aexx,
               is_magnetization=args.is_magnetization,
               poscar=args.poscar)


def vasp_poscar_maker(args):

    if args.supercell:
        make_supercell_poscar(args.supercell, args.poscar, args.sposcar)
    else:
        make_hpkot_primitive_poscar(poscar=args.poscar,
                                    pposcar=args.pposcar,
                                    symprec=args.symprec,
                                    angle_tolerance=args.angle_tolerance)


def vasp_potcar_maker(args):
    if args.elements:
        make_potcar(args.path, args.elements)
    elif args.poscar:
        elements = Structure.from_file(args.poscar).symbol_set
        make_potcar(args.path, elements)


def recommend_supercell(args):
    if not args.poscar_path:
        raise ValueError("Specify path of poscar!")
    kw_args = {}
    if args.criterion:
        kw_args["criterion"] = args.criterion
    if args.min_natom:
        kw_args["min_natom"] = args.min_natom
    if args.max_natom:
        kw_args["max_natom"] = args.max_natom
    recommend_supercell_ase(args.poscar_path, **kw_args)


def defect_entry(args):
    if args.make_defect_entry:
        defect_entry_from_yaml = DefectEntry.from_yaml(args.yaml)
        defect_entry_from_yaml.to_json_file("defect_entry.json")
    elif args.print:
        print(DefectEntry.load_json(args.json))


def supercell_results(args):
    if args.print:
        print(SupercellDftResults.load_json(args.json))
        return True

    if args.dir_all:
        dirs = glob('*[0-9]/')
        dirs.append("perfect/")
    else:
        dirs = args.dirs

    for d in dirs:
        if os.path.isdir(d):
            if args.convergence:
                check = check_vasp_output(d,
                                          contcar_name=args.poscar,
                                          outcar_name=args.outcar,
                                          vasprun_name=args.vasprun)
                if check["all"]:
                    if vasp_convergence_ionic(d, args.vasprun):
                        ionic = "Y"
                    else:
                        ionic = "N"
                    if vasp_convergence_electronic(d, args.vasprun):
                        electronic = "Y"
                    else:
                        electronic = "N"
                    print("{:>20}  ionic:{:>3}  electronic:{:>3}".
                          format(d, ionic, electronic))
                else:
                    contcar = "Y"
                    outcar = "Y"
                    vasprun = "Y"

                    if check["contcar"] is None:
                        contcar = "NA"
                    elif check["contcar"] is False:
                        contcar = "N"

                    if check["outcar"] is None:
                        outcar = "NA"
                    elif check["outcar"] is False:
                        outcar = "N"

                    if check["vasprun"] is None:
                        vasprun = "NA"
                    elif check["vasprun"] is False:
                        vasprun = "N"

                    print("{:>20}  CONTCAR:{:>3}  OUTCAR:{:>3}, vasprun:{:>3}".
                          format(d, contcar, outcar, vasprun))

            else:
                print(d)
                try:
                    dft_results = SupercellDftResults.from_vasp_files(
                        d,
                        contcar_name=args.poscar,
                        outcar_name=args.outcar,
                        vasprun_name=args.vasprun)

                    dft_results.to_json_file(
                        filename=os.path.join(d, "dft_results.json"))
                except:
                    warnings.warn(
                        message="Parsing data in " + d + " is failed.")
        else:
            warnings.warn(message=d + " does not exist, so nothing is done.")


def unitcell_results(args):
    try:
        dft_results = UnitcellDftResults.load_json(filename=args.json_file)
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
            dft_results. \
                set_static_dielectric_tensor_from_vasp(args.static_diele_dir,
                                                       outcar_name=args.outcar)
        except IOError:
            print(args.static_diele_dir, "is not appropriate.")

    if args.ionic_diele:
        dft_results.ionic_dielectric_tensor = args.ionic_diele
    elif args.ionic_diele_dir:
        try:
            dft_results. \
                set_ionic_dielectric_tensor_from_vasp(args.ionic_diele_dir,
                                                      outcar_name=args.outcar)
        except IOError:
            print(args.ionic_diele_dir, "is not appropriate.")

    if args.volume_dir:
        try:
            dft_results.set_volume_from_vasp(args.volume_dir,
                                             contcar_name=args.poscar)
        except IOError:
            print(args.volume_dir, "is not appropriate.")

    if args.total_dos_dir:
        try:
            dft_results.set_total_dos_from_vasp(args.total_dos_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            print(args.total_dos_dir, "is not appropriate.")

    dft_results.to_json_file(args.json_file)


def correction(args):
    try:
        unitcell_dft_data = UnitcellDftResults.load_json(args.unitcell_json)
    except IOError:
        raise FileNotFoundError("JSON of unitcell was not found.")

    try:
        perfect_dft_data = SupercellDftResults.load_json(args.perfect_json)
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
        json_to_make = os.path.join(directory, "correction.json")
        if os.path.exists(json_to_make) and not args.force_overwrite:
            print("{} exists. Correction was not done.".format(json_to_make))
            continue
        print("correcting {0} ...".format(directory))
        try:
            entry = DefectEntry.load_json(
                os.path.join(directory, "defect_entry.json"))
            defect_dft_data = \
                SupercellDftResults.load_json(
                    os.path.join(directory, "dft_results.json"))
            c = Correction.compute_alignment_by_extended_fnv(entry,
                                                             defect_dft_data,
                                                             perfect_dft_data,
                                                             unitcell_dft_data,
                                                             ewald_data)
            c.plot_distance_vs_potential(
                file_name=os.path.join(directory, "potential.eps"))
            c.to_json_file(os.path.join(directory, "correction.json"))
        except Exception as e:
            warnings.warn("Correction for {0} is failed. "
                          "The calculation for {0} is skipped."
                          "Exception type and message is {1}, {2}".
                          format(directory, type(e), e.args))


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
            cp = ChemPotDiag.from_file(args.energy_file)
        if args.vasp_dirs:
            poscar_paths = [d + args.poscar_name for d in args.vasp_dirs]
            outcar_paths = [d + args.outcar_name for d in args.vasp_dirs]
            cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                          outcar_paths)
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

    unitcell = UnitcellDftResults.load_json(args.unitcell)
    perfect = SupercellDftResults.load_json(args.perfect)

    if not args.defect_dirs:
        from glob import glob
        defects_dirs = glob('*[0-9]/')
    else:
        defects_dirs = args.defects

    defects = []
    for d in defects_dirs:
        try:

            de = DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
            dr = SupercellDftResults.\
                load_json(os.path.join(d, "dft_results.json"))
            co = Correction.load_json(os.path.join(d, "correction.json"))

            defects.append(Defect(defect_entry=de,
                                  dft_results=dr,
                                  correction=co))
        except:
            warnings.warn(message="Parsing data in " + d + " is failed.")

    chem_pot = ChemPotDiag.load_vertices_yaml(args.chem_pot_yaml)

    defect_energies = DefectEnergies(unitcell=unitcell,
                                     perfect=perfect,
                                     defects=defects,
                                     chem_pot=chem_pot,
                                     chem_pot_label=args.chem_pot_label,
                                     system_name=args.name)

    defect_energies.calc_transition_levels()
    defect_energies.plot_energy()


if __name__ == "__main__":
    main()
