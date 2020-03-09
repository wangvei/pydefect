# -*- coding: utf-8 -*-

import os
import shutil
from copy import deepcopy
from glob import glob
from inspect import signature
from itertools import chain
from pathlib import Path

import numpy as np

from pydefect.analysis.defect import Defect
from pydefect.analysis.defect_carrier_concentration import DefectConcentration
from pydefect.analysis.defect_eigenvalues import DefectEigenvalue
from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.analysis.defect_structure import (
    DefectStructure, defect_structure_matcher)
from pydefect.core.complex_defects import ComplexDefects
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.error_classes import StructureError
from pydefect.core.interstitial_site import (
    InterstitialSiteSet, interstitials_from_charge_density)
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.corrections.corrections import ManualCorrection
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection, Ewald
from pydefect.corrections.vertical_transition_energy_correction import (
    VerticalTransitionEnergyCorrection)
from pydefect.input_maker.defect_initial_setting import (
    dopant_info, DefectInitialSetting)
from pydefect.input_maker.supercell_maker import Supercell, Supercells
from pydefect.util.logger import get_logger
from pydefect.cli.main_tools import (
    generate_objects_from_json_files)

from pymatgen import Structure, Spin
from pymatgen.core.periodic_table import Element

from vise.input_set.incar import incar_flags
from vise.input_set.input_set import ViseInputSet
from vise.cli.main_function import vasp_settings_from_args
from vise.cli.main_tools import potcar_str2dict, list2dict
from vise.chempotdiag.chem_pot_diag import ChemPotDiag
logger = get_logger(__name__)


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
                                                vasprun_name=args.vasprun,
                                                outcar_name=args.outcar)
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
        except AttributeError as e:
            logger.error(str(e))

    if args.ionic_diele:
        dft_results.ionic_dielectric_tensor = args.ionic_diele
    elif args.ionic_diele_dir:
        try:
            dft_results.set_ionic_dielectric_tensor_from_vasp(
                args.ionic_diele_dir, outcar_name=args.outcar)
        except IOError:
            raise FileNotFoundError(args.ionic_diele_dir, "not appropriate.")

    if args.volume_dir:
        try:
            dft_results.set_volume_from_vasp(
                args.volume_dir, contcar_name=args.contcar)
        except IOError:
            raise FileNotFoundError(args.volume_dir, "not appropriate.")

    if args.total_dos_dir:
        try:
            dft_results.set_total_dos_from_vasp(args.total_dos_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            raise FileNotFoundError(args.total_dos_dir, "not appropriate.")

    dft_results.to_json_file(args.json_file)
    print(UnitcellCalcResults.load_json(filename=args.json_file))


def initial_setting(args):
    if args.print_dopant:
        print(dopant_info(args.print_dopant))
        return

    structure = Structure.from_file(args.poscar)

    kwargs = {"dopants": args.dopants,
              "is_antisite": args.antisite,
              "en_diff": args.en_diff,
              "included": args.included,
              "excluded": args.excluded,
              "displacement_distance": args.displacement_distance,
              "symprec": args.symprec,
              "angle_tolerance": args.angle_tolerance,
              "interstitial_sites": args.interstitials,
              "complex_defect_names": args.complex_defect_names}

    if args.matrix:
        supercell = Supercell(structure=structure,
                              trans_mat=args.matrix,
                              check_unitcell=True,
                              symprec=args.symprec,
                              angle_tolerance=args.angle_tolerance)
        supercell.to(poscar="DPOSCAR", uposcar="UPOSCAR")

        defect_setting = \
            DefectInitialSetting.from_basic_settings(
                structure=supercell.structure,
                transformation_matrix=supercell.trans_mat.tolist(),
                cell_multiplicity=supercell.multiplicity, **kwargs)
        defect_setting.to()
        return

    supercells = Supercells(structure=structure,
                            conventional_base=args.conventional_base,
                            max_num_atoms=args.max_num_atoms,
                            min_num_atoms=args.min_num_atoms,
                            criterion=args.isotropy_criterion,
                            rhombohedral_angle=args.rhombohedral_angle,
                            symprec=args.symprec,
                            angle_tolerance=args.angle_tolerance)

    if not supercells.supercells:
        logger.warning("No supercell satisfies the criterion.")
        return False

    unitcell = supercells.unitcell
    if unitcell != structure:
        logger.warning(
            "Unitcell is different from input structure, so generate UPOSCAR.")
        supercells.to_uposcar(uposcar="UPOSCAR")
    else:
        logger.info("Input structure is the primitive cell.")

    if not args.supercell_set:
        if args.most_isotropic:
            supercell = supercells.most_isotropic_supercell
        else:
            supercell = supercells.smallest_supercell

        supercell.to("DPOSCAR")
        defect_setting = \
            DefectInitialSetting.from_basic_settings(
                structure=supercell.structure,
                transformation_matrix=supercell.trans_mat.tolist(),
                cell_multiplicity=supercell.multiplicity, **kwargs)
        defect_setting.to()

    else:
        logger.info(f"Number of supercells: {len(supercells.supercells)}")

        for supercell in supercells.supercells:
            # Suffix "c" means conventional cell, while "p" primitive cell.
            prefix = "c" if supercells.conventional_base else "p"
            isotropy, _ = supercell.isotropy
            if np.count_nonzero(supercell.trans_mat) == 3:
                mat = "x".join(
                    [str(supercell.trans_mat[i][i]) for i in range(3)])
                name = f"{prefix + mat}_{supercell.num_atoms}_{isotropy}"
            else:
                name = f"{prefix}_{supercell.num_atoms}_{isotropy}"

            os.mkdir(name)
            p = Path(name)
            supercell.to(poscar=p / "DPOSCAR")
            defect_setting = \
                DefectInitialSetting.from_basic_settings(
                    structure=supercell.structure,
                    transformation_matrix=supercell.trans_mat.tolist(),
                    cell_multiplicity=supercell.multiplicity, **kwargs)
            defect_setting.to(p / "defect.in")


def interstitial(args):

    if args.chgcar:
        interstitials_from_charge_density(
            chgcar_filename=args.chgcar,
            interstitial_symprec=args.interstitial_symprec,
            angle_tolerance=args.angle_tolerance)

    else:
        try:
            interstitial_set = \
                InterstitialSiteSet.from_files(args.dposcar, args.yaml)
        except FileNotFoundError:
            structure = Structure.from_file(args.dposcar)
            interstitial_set = InterstitialSiteSet(structure=structure)

        coords = args.interstitial_coords
        if len(coords) == 3:
            coords = [coords]
        elif len(coords) % 3 == 0:
            n_coords = int(len(coords) / 3)
            coords = \
                [[coords[3 * i + j] for j in range(3)] for i in range(n_coords)]
        else:
            raise ValueError(f"Interstitial coordinates "
                             f"{args.interstitial_coords} are invalid")

        defect_initial_setting = \
            DefectInitialSetting.from_defect_in(poscar=args.dposcar,
                                                defect_in_file=args.defect_in)
        # To change coords from unitcell to supercell, multiply inverse of
        # trans_mat to coords.
        tm_list = defect_initial_setting.transformation_matrix
        trans_mat = [[tm_list[3 * i + j] for j in range(3)] for i in range(3)]
        inv_trans_mat = np.linalg.inv(trans_mat)
#        supercell_coords = [np.dot(inv_trans_mat, c).tolist() for c in coords]
#        print(supercell_coords)
        # TODO: understand the reason why inv_trans_mat must be multiplied from
        #  right.
        supercell_coords = [np.dot(c, inv_trans_mat).tolist() for c in coords]

        interstitial_set.add_sites(frac_coords=supercell_coords,
                                   vicinage_radius=args.radius,
                                   defect_symprec=args.symprec,
                                   angle_tolerance=args.angle_tolerance)

        interstitial_set.site_set_to_yaml_file(yaml_filename=args.yaml)


def complex_defects(args):
    try:
        complex_defects_obj = ComplexDefects.from_files(args.dposcar, args.yaml)
    except FileNotFoundError:
        structure = Structure.from_file(args.dposcar)
        complex_defects_obj = ComplexDefects(structure=structure)

    if (args.inserted_elements and not args.inserted_coords) or \
            (not args.inserted_elements and args.inserted_coords):
        raise ValueError(f"For interstitial sites, both elements and "
                         f"fractional coordinates need to be entered.")
    elif args.inserted_elements and args.inserted_coords and \
            len(args.inserted_elements) * 3 != len(args.inserted_coords):
        raise ValueError(f"The numbers of inserted elements "
                         f"{args.inserted_elements} and coords "
                         f"{args.inserted_coords} are invalid")

    inserted_elements = args.inserted_elements if args.inserted_elements else []
    inserted_coords = args.inserted_coords if args.inserted_coords else []

    inserted_atoms = []
    for i, e in enumerate(inserted_elements):
        coords = [inserted_coords[3 * i + j] for j in range(3)]
        inserted_atoms.append({"element": e, "coords": coords})

    # defect_initial_setting = \
    #     DefectInitialSetting.from_defect_in(poscar=args.dposcar,
    #                                         defect_in_file=args.defect_in)

    complex_defects_obj.add_defect(
        removed_atom_indices=args.removed_atom_indices,
        inserted_atoms=inserted_atoms,
        name=args.name,
        extreme_charge_state=args.extreme_charge_state,
        annotation=args.annotation,
        symprec=args.defect_symprec,
        angle_tolerance=args.angle_tolerance)

    complex_defects_obj.site_set_to_yaml_file(yaml_filename=args.yaml)


def make_dir(name: str, vis: ViseInputSet, force_overwrite: bool) -> None:
    """Helper function"""
    if force_overwrite and os.path.exists(name):
        logger.warning(f"{name:>10} is being removed.")
        shutil.rmtree(name)

    if os.path.exists(name):
        logger.warning(f"{name:>10} already exists, so nothing is done.")
    else:
        logger.warning(f"{name:>10} is being constructed.")
        os.makedirs(name)
        vis.write_input(name)
        vis.to_json_file("/".join([name, "vise.json"]))


def defect_vasp_set(args):

    user_incar_settings, vis_kwargs = vasp_settings_from_args(args)
    user_incar_settings["LWAVE"] = args.wavecar
    vis_kwargs.update(
        {"xc": args.xc,
         "task": "defect",
         "kpt_density": args.kpt_density,
         "kpt_mode": "manual_set",
         "only_even": False,
         "sort_structure": False,
         "standardize_structure": False,
         })

    defect_initial_setting = DefectInitialSetting.from_defect_in(
        poscar=args.dposcar, defect_in_file=args.defect_in)

    defect_initial_setting.make_defect_set(
        keywords=args.keywords, specified_defects=args.specified_defects)

    if not args.specified_defects:
        perfect_incar_setting = deepcopy(user_incar_settings)
        vise_set = ViseInputSet.make_input(
            structure=defect_initial_setting.structure,
            user_incar_settings=perfect_incar_setting,
            **vis_kwargs)

        make_dir("perfect", vise_set, args.force_overwrite)

    if args.spin_polarize:
        user_incar_settings["ISPIN"] = 2

    for de in defect_initial_setting.defect_entries:
        defect_name = "_".join([de.name, str(de.charge)])
        json_file_name = os.path.join(defect_name, "defect_entry.json")

        vise_set = ViseInputSet.make_input(
            structure=de.perturbed_initial_structure,
            charge=de.charge,
            user_incar_settings=user_incar_settings,
            **vis_kwargs)

        make_dir(defect_name, vise_set, args.force_overwrite)
        de.to_json_file(json_file_name)

        if de.neighboring_sites:
            poscar_name = os.path.join(defect_name, "POSCAR")
            with open(poscar_name, "r") as f:
                lines = f.readlines()
                for index, line in enumerate(lines.copy()):
                    if index - 8 in de.neighboring_sites:
                        lines[index] = line.strip() + "  Neighbor\n"

            with open(poscar_name, "w") as f:
                for line in lines:
                    f.write(line)


def vertical_transition_input_maker(args):
    if abs(args.additional_charge) != 1:
        raise ValueError(f"{args.additional_charge} is invalid.")

    initial_dirname = Path(args.initial_dir_name)
    de_filename = initial_dirname / "defect_entry.json"
    de = DefectEntry.load_json(de_filename)
    src_filename = initial_dirname / "dft_results.json"
    src = SupercellCalcResults.load_json(src_filename)

    new_charge = de.charge + args.additional_charge
    vis = ViseInputSet.from_prev_calc(
        initial_dirname,
        user_incar_settings={"NSW": 0},
        parse_calc_results=False,
        contcar_filename=args.contcar,
        charge=new_charge)

    de.charge += args.additional_charge
    de.initial_structure = src.final_structure
    de.perturbed_initial_structure = src.final_structure
    de.initial_site_symmetry = src.site_symmetry
    # FIX MAGNETIZATION?
    new_dirname = initial_dirname / f"add_charge_{args.additional_charge}"
    make_dir(str(new_dirname), vis, force_overwrite=False)

    new_de_filename = new_dirname / "defect_entry.json"
    de.to_json_file(new_de_filename)


def defect_entry(args):
    if args.print:
        print(DefectEntry.load_json(args.json))
    elif args.make_defect_entry:
        defect_structure = Structure.from_file(args.defect_poscar)
        perfect_structure = Structure.from_file(args.perfect_poscar)

        defect_entry_from_yaml = DefectEntry.from_defect_structure(
            defect_structure=defect_structure,
            perfect_structure=perfect_structure,
            displacement_distance=args.displacement_distance,
            defect_name=args.defect_name)

        defect_entry_from_yaml.to_json_file(args.json)
    else:
        logger.warning("Set make_defect_entry or print option.")


def supercell_calc_results(args):

    if args.print:
        print(SupercellCalcResults.load_json(args.json))
        return

    if args.defect_center:
        if len(args.defect_center) != 1 and len(args.defect_center) != 3:
            raise ValueError("Length of the defect center is neither 1 or 3")
        results = SupercellCalcResults.load_json(args.json)
        results.defect_center = args.defect_center
        results.to_json_file(filename=args.json)
        return

    if args.dir_all:
        dirs = glob('*[0-9]/')
        dirs.insert(0, "perfect/")
    else:
        dirs = args.dirs

    for d in dirs:
        if os.path.isdir(d):
            logger.info(f"Parsing data in {d} ...")

            if d in ["perfect", "perfect/"]:
                try:
                    dft_results = SupercellCalcResults.from_vasp_files(
                        directory_path=d,
                        vasprun=args.vasprun,
                        contcar=args.contcar,
                        procar=args.procar,
                        outcar=args.outcar)
                except IOError:
                    raise IOError("Parsing data in perfect failed.")
            else:
                try:
                    de = DefectEntry.load_json(
                        os.path.join(d, args.defect_entry_name))

                    dft_results = \
                        SupercellCalcResults.from_vasp_files(
                            directory_path=d,
                            vasprun=args.vasprun,
                            contcar=args.contcar,
                            outcar=args.outcar,
                            procar=args.procar,
                            cutoff=args.cutoff,
                            defect_entry=de,
                            defect_symprec=args.defect_symprec,
                            angle_tolerance=args.angle_tolerance)
                except IOError:
                    logger.warning(f"Parsing data in {d} failed.")
                    continue

            dft_results.to_json_file(
                filename=os.path.join(d, "dft_results.json"))
        else:
            logger.warning(f"{d} does not exist, so nothing is done.")


def efnv_correction(args):
    if args.print:
        print(ExtendedFnvCorrection.load_json(args.json_file))
        return

    dirs = glob('*[0-9]/') if args.dir_all else args.dirs

    if args.plot_potential:
        for directory in dirs:
            json_file = os.path.join(directory, "correction.json")
            c = ExtendedFnvCorrection.load_json(json_file)
            c.plot_potential(os.path.join(directory, "potential.pdf"),
                             args.y_range)
        return

    if args.nocorr:
        for directory in dirs:
            c = ManualCorrection(manual_correction_energy=args.manual)
            c.to_json_file(os.path.join(directory, "correction.json"))
        return

    try:
        ucr = UnitcellCalcResults.load_json(args.unitcell_json)
        dielectric_tensor = ucr.total_dielectric_tensor
    except IOError:
        raise FileNotFoundError("JSON for the unitcell info is not found.")

    try:
        perfect_dft_data = SupercellCalcResults.load_json(args.perfect_json)
    except IOError:
        raise FileNotFoundError("JSON for the perfect supercell is not found.")

    # Ewald parameter related
    if not Path(args.ewald_json).is_file():
        logger.info("optimizing ewald...")
        ewald = Ewald.from_optimization(
            structure=perfect_dft_data.final_structure,
            dielectric_tensor=dielectric_tensor,
            initial_ewald_param=args.ewald_initial_param,
            convergence=args.ewald_convergence,
            prod_cutoff_fwhm=args.ewald_accuracy)
        ewald.to_json_file(args.ewald_json)

    for directory in dirs:
        json_to_make = os.path.join(directory, "correction.json")

        if os.path.exists(json_to_make) and not args.force_overwrite:
            logger.warning(f"{json_to_make} already exists, so nothing done.")
            continue

        logger.info(f"correcting {directory} ...")
        entry = DefectEntry.load_json(os.path.join(directory,
                                                   "defect_entry.json"))
        try:
            defect_dft_data = SupercellCalcResults.load_json(
                os.path.join(directory, "dft_results.json"))
        except IOError:
            logger.warning(f"dft_results.json in {directory} does not exist.")
            continue

        c = ExtendedFnvCorrection. \
            compute_correction(defect_entry=entry,
                               defect_dft=defect_dft_data,
                               perfect_dft=perfect_dft_data,
                               dielectric_tensor=dielectric_tensor,
                               defect_center=args.defect_center,
                               ewald=args.ewald_json)

        c.plot_potential(os.path.join(directory, "potential.pdf"),
                         args.y_range)
        c.to_json_file(os.path.join(directory, "correction.json"))


def vertical_transition_energy(args):

    initial_dir = Path(args.initial_dir)
    initial_calc_results = \
        SupercellCalcResults.load_json(initial_dir / "dft_results.json")
    final_dir = Path(args.dir)
    final_calc_results = \
        SupercellCalcResults.load_json(final_dir / "dft_results.json")
    unitcell = UnitcellCalcResults.load_json(args.unitcell_json)

    if args.print:
        vtec = VerticalTransitionEnergyCorrection.load_json(args.json)
        vtec.plot_potential()
        print(vtec)
        if vtec.additional_charge == 1:
            cbm = unitcell.band_edge[1]
            print(f"CBM position (eV): {cbm}")
            band_edge_related_energy = cbm

        else:
            vbm = unitcell.band_edge[0]
            print(f"VBM position (eV): {vbm}")
            band_edge_related_energy = -vbm

        vte_wo_corr = (final_calc_results.total_energy
                       - initial_calc_results.total_energy
                       + band_edge_related_energy)
        vte = vte_wo_corr + vtec.correction_energy
        print(f"Vertical transition energy w/o correction (eV): {vte_wo_corr}")
        print(f"Vertical transition energy w/  correction (eV): {vte}")
        return

    dielectric_tensor = unitcell.total_dielectric_tensor
    static_dielectric_tensor = unitcell.static_dielectric_tensor

    initial_efnv = \
        ExtendedFnvCorrection.load_json(initial_dir / "correction.json")
    initial_calc_results = \
        SupercellCalcResults.load_json(initial_dir / "dft_results.json")

    final_defect_entry = DefectEntry.load_json(final_dir / "defect_entry.json")

    c = VerticalTransitionEnergyCorrection.from_files(dielectric_tensor,
                                                      static_dielectric_tensor,
                                                      initial_efnv,
                                                      initial_calc_results,
                                                      final_defect_entry,
                                                      final_calc_results)
    c.to_json_file(args.json)
    print(c)


def defects(args):

    if args.band_edge:
        if args.band_edge[0] == "up":
            spin = Spin.up
        elif args.band_edge[0] == "down":
            spin = Spin.down
        else:
            raise ValueError("band edge is inadequate (e.g. -be up no_in_gap).")
        state = args.band_edge[1]
        defect = Defect.load_json(args.json)
        defect.set_band_edge_state(spin=spin, state=state)
        defect.to_json_file(args.json)
        return True

    try:
        perfect = SupercellCalcResults.load_json(args.perfect)
    except FileNotFoundError:
        print(f"{args.perfect} not found.")
        raise

    defects_dirs = args.defect_dirs or glob('*[0-9]/')
    for d in defects_dirs:
        filename = os.path.join(d, args.json)
        if args.diagnose:
            print(d.rjust(12), end="  ")
            try:
                print(Defect.load_json(filename).diagnose)
            except FileNotFoundError:
                logger.warning("No supercell results file.")
            except Exception as e:
                logger.warning(f"An error {e} is caught.")
            continue

        logger.info(f"parsing directory {d}...")
        files = [args.defect_entry, args.dft_results, args.correction]
        classes = [DefectEntry, SupercellCalcResults, ExtendedFnvCorrection]
        input_objects = generate_objects_from_json_files(d, files, classes,
                                                         raise_error=False)
        if input_objects:
            try:
                defect = Defect.from_objects(defect_entry=input_objects[0],
                                             dft_results=input_objects[1],
                                             perfect_dft_results=perfect,
                                             correction=input_objects[2])
                defect.to_json_file(filename)
            except StructureError:
                logger.warning(f"defect.json is not generated in {d}.")

        else:
            logger.warning(f"Generating {filename} failed.")
            continue


def plot_energy(args):

    if args.reload_defects:
        os.remove(args.energies)

    try:
        defect_energies = DefectEnergies.load_json(args.energies)
    except FileNotFoundError:
        unitcell = UnitcellCalcResults.load_json(args.unitcell)
        perfect = SupercellCalcResults.load_json(args.perfect)

        defects_dirs = args.defect_dirs or glob('*[0-9]/')

        defect_list = []
        for d in defects_dirs:
            filename = os.path.join(d, args.defect)
            logger.info(f"parsing directory {d}...")
            try:
                defect_list.append(Defect.load_json(filename))
            except FileNotFoundError:
                logger.warning(f"Parsing {filename} failed.")
                continue

        chem_pot = ChemPotDiag.load_json(args.chem_pot_json)

        # First construct DefectEnergies class object.
        defect_energies = \
            DefectEnergies.from_objects(unitcell=unitcell,
                                        perfect=perfect,
                                        defects=defect_list,
                                        chem_pot=chem_pot,
                                        chem_pot_label=args.chem_pot_label,
                                        filtering_words=args.filtering,
                                        system=args.name)

    if args.print:
        print(defect_energies)
        return

    defect_energies.to_json_file(filename=args.energies)

    # if args.concentration:
    #     defect_concentration = \
    #         DefectConcentration.from_defect_energies(
    #             energies=energies,
    #             temperature=args.temperature[0],
    #             unitcell=unitcell,
    #             num_sites_filename=args.num_site_file)

 #  #     if len(args.temperature) == 2:
    #         defect_concentration = \
    #             DefectConcentration.from_defect_energies(
    #                 energies=energies,
    #                 temperature=args.temperature[1],
    #                 unitcell=unitcell,
    #                 num_sites_filename=args.num_site_file,
    #                 previous_concentration=defect_concentration)
    # else:
#         defect_concentration = None

    # plt = defect_energies.plot_energy(filtering_words=args.filtering,
    #                                   x_range=args.x_range,
    #                                   y_range=args.y_range,
   #                                   show_transition_levels=args.show_tl,
   #                                   show_all_energies=args.show_all)
    plt = defect_energies.plot_energy(
        x_range=args.x_range,
        y_range=args.y_range,
        show_transition_levels=args.show_transition_level,
        show_all_energies=args.show_all)

    if args.save_file:
        plt.savefig(args.save_file, format="pdf", transparent=True)
    else:
        plt.show()


def parse_eigenvalues(args):
    unitcell = UnitcellCalcResults.load_json(args.unitcell)

    # TODO: Modify here to run w/o correction
    logger.info(f"parsing directory {args.defect_dir}...")
    filename = os.path.join(args.defect_dir, args.defect)
    defect = Defect.load_json(filename)

    defect_eigenvalues = DefectEigenvalue.from_files(unitcell=unitcell,
                                                     defect=defect)

    defect_eigenvalues.plot(y_range=args.y_range,
                            title=args.title,
                            filename=args.save_file)


def vasp_parchg_set(args):
    user_incar_settings = {"LPARD": True,
                           "LSEPB": True,
                           "KPAR": 1,
                           "IBAND": args.band_indices}

    if args.kpoint_indices:
        user_incar_settings["KPUSE"] = args.kpoint_indices

    vasp_set = ViseInputSet.from_prev_calc(
        dirname=args.read_dir,
        parse_calc_results=False,
        parse_incar=True,
        sort_structure=False,
        standardize_structure=False,
        files_to_transfer={"WAVECAR": "L"},
        user_incar_settings=user_incar_settings,
        contcar_filename=args.contcar)

    vasp_set.write_input(args.write_dir)


def local_structure(args):
    defects_dirs = args.defect_dirs if args.defect_dirs else glob('*[0-9]/')

    d_list = []
    for d in defects_dirs:
        filename = os.path.join(d, args.defect)
        logger.info(f"parsing directory {d}...")
        try:
            defect = Defect.load_json(filename)
            d_list.append(defect)
        except FileNotFoundError:
            logger.warning(f"Parsing {filename} failed.")
            print("")
            print("")
            continue

        defect_structure = DefectStructure.from_defect(defect)
        print("-" * 82)
        print(defect_structure.show_displacements(all_atoms=args.show_all))
        print("")
        print("")

    if args.compare_structure:
        print("-" * 82)
        print(defect_structure_matcher(d_list, args.site_tolerance))
        print("")


def concentration(args):
    defect_energies = DefectEnergies.load_json(args.energies).defect_energies
    unitcell = UnitcellCalcResults.load_json(args.unitcell)

    defect_concentration = DefectConcentration.from_calc_results(
        defect_energies=defect_energies,
        unitcell=unitcell,
        round_magnetization=args.fmto)

    defect_concentration.calc_equilibrium_concentration(
        temperature=args.temperature, verbose=args.verbose)

    defect_concentration.calc_quenched_equilibrium_concentration(
        verbose=args.verbose)

    print(defect_concentration)
    defect_concentration.calc_concentrations(args.temperature)
    plt = defect_concentration.plot_carrier_concentrations()
    plt.show()

