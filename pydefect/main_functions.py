# -*- coding: utf-8 -*-

import os
import shutil
from copy import deepcopy
from glob import glob
from inspect import signature
from itertools import chain
from os.path import join

import numpy as np
from chempotdiag.chem_pot_diag import ChemPotDiag
from vise.input_set.incar import incar_flags
from vise.input_set.input_set import ObaSet
from pydefect.analysis.defect import Defect
from pydefect.analysis.defect_carrier_concentration import DefectConcentration
from pydefect.analysis.defect_eigenvalues import DefectEigenvalue
from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.analysis.defect_structure import DefectStructure
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.interstitial_site import (
    InterstitialSiteSet, interstitials_from_charge_density)
from pydefect.core.prior_info import PriorInfo
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection, Ewald
from pydefect.corrections.corrections import ManualCorrection
from pydefect.input_maker.defect_initial_setting import (
    dopant_info, DefectInitialSetting)
from pydefect.input_maker.supercell_maker import Supercells
from pydefect.util.logger import get_logger
from pydefect.util.main_tools import list2dict, generate_objects
from pymatgen import Structure, Spin
from pymatgen.core.periodic_table import Element

__author__ = "Yu Kumagai, Akira Takahashi"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def recommend_supercell(args):
    supercells = Supercells(structure=Structure.from_file(args.poscar),
                            conventional_base=not args.primitive,
                            max_num_atoms=args.max_num_atoms,
                            min_num_atoms=args.min_num_atoms,
                            criterion=args.isotropy_criterion)

    if supercells.supercells:
        if args.set:
            logger.info(f"Number of supercells: {len(supercells.supercells)}")

            for supercell in supercells.supercells:
                # Suffix "c" means conventional cell, while "p" primitive cell.
                prefix = "c" if supercells.conventional_base else "p"
                if np.count_nonzero(supercell.trans_mat) == 3:
                    mat = [str(supercell.trans_mat[i][i]) for i in range(3)]
                else:
                    mat = [str(j) for i in supercell.trans_mat for j in i]
                name = "_".join([args.sposcar,
                                 prefix + "x".join(mat),
                                 str(supercell.num_atoms),
                                 str(supercell.isotropy[0])])
                supercell.to_poscar(poscar_filename=name)

            supercells.supercells[0].to_uposcar(uposcar_filename=args.uposcar)

        else:
            if args.most_isotropic:
                supercell = supercells.most_isotropic_supercell
            else:
                supercell = supercells.smallest_supercell

            supercell.to_poscar(poscar_filename=args.sposcar)
            supercell.to_uposcar(uposcar_filename=args.uposcar)

    else:
        logger.warning("No supercell satisfies the criterion.")


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
                interstitial_sites=args.interstitials)

        defect_setting.to()


def interstitial(args):

    if args.chgcar:
        interstitials_from_charge_density(
            chgcar_filename=args.chgcar,
            symprec=args.symprec,
            angle_tol=args.angle_tolerance)

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
            length = int(len(coords) / 3)
            coords = [[coords[3 * i + j]
                       for j in range(3)] for i in range(length)]
        else:
            raise ValueError(
                f"Interstitial coordinates {args.interstitial_coords} invalid")

        dis = DefectInitialSetting.from_defect_in(poscar=args.dposcar,
                                                  defect_in_file="defect.in")
        trans_mat = [[dis.transformation_matrix[3 * i + j]
                      for j in range(3)] for i in range(3)]
        # Change coords from unitcell to supercell
        # multiply inverse of trans_mat to coords
        inv_trans_mat = np.linalg.inv(trans_mat)
        supercell_coords = [np.dot(inv_trans_mat, c).tolist() for c in coords]

        interstitial_set.add_sites(coords=supercell_coords,
                                   vicinage_radius=args.radius,
                                   symprec=args.symprec,
                                   angle_tol=args.angle_tolerance)

        interstitial_set.site_set_to_yaml_file(filename=args.yaml)


def defect_vasp_oba_set(args):

    flags = list(signature(ObaSet.make_input).parameters.keys())
    kwargs = list2dict(args.vos_kwargs, flags)

    def make_dir(name, obrs):
        """Helper function"""
        if args.force_overwrite and os.path.exists(name):
            logger.warning(f"{name:>10} is being removed.")
            shutil.rmtree(name)

        if os.path.exists(name):
            logger.warning(f"{name:>10} already exists, so nothing is done.")
        else:
            logger.warning(f"{name:>10} is being constructed.")
            os.makedirs(name)
            obrs.write_input(name)

    dis = DefectInitialSetting.from_defect_in(
        poscar=args.dposcar, defect_in_file=args.defect_in)
    dis.make_defect_set(keywords=args.keywords,
                        specified_defects=args.particular_defects)

    if not args.particular_defects:
        oba_set = ObaSet.make_input(
            structure=dis.structure,
            standardize_structure=False,
            task="defect",
            xc=args.xc,
            additional_user_potcar_yaml=args.potcar,
            sort_structure=False,
            weak_incar_settings={"LWAVE": args.wavecar},
            kpt_mode="manual",
            kpt_density=args.kpt_density,
            only_even=False,
            user_incar_settings={"ISPIN": 1},
            **kwargs)

        make_dir("perfect", oba_set)

    for de in dis.defect_entries:
        defect_name = "_".join([de.name, str(de.charge)])
        json_file_name = os.path.join(defect_name, "defect_entry.json")

        oba_set = ObaSet.make_input(structure=de.perturbed_initial_structure,
                                    charge=de.charge,
                                    standardize_structure=False,
                                    task="defect",
                                    xc=args.xc,
                                    additional_user_potcar_yaml=args.potcar,
                                    sort_structure=False,
                                    weak_incar_settings=
                                    {"LWAVE": args.wavecar},
                                    kpt_mode="manual",
                                    kpt_density=args.kpt_density,
                                    only_even=False,
                                    **kwargs)

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
    if args.print:
        print(DefectEntry.load_json(args.json))
    elif args.make_defect_entry:
        defect_entry_from_yaml = \
            DefectEntry.from_yaml(filename=args.yaml,
                                  cutoff=args.cutoff,
                                  calc_num_equiv_site=args.calc_sites)
        defect_entry_from_yaml.to_json_file(args.json)
    else:
        logger.warning("Set make_defect_entry or print option.")


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
                except IOError:
                    raise IOError("Parsing data in perfect failed")
            else:
                try:
    #                filename = join(args.perfect_results, "dft_results.json")
    #                perfect_results = SupercellCalcResults.load_json(filename)
                    de = DefectEntry.load_json(join(d, args.defect_entry_name))

                    dft_results = \
                        SupercellCalcResults.from_vasp_files(
                            directory_path=d,
                            vasprun=args.vasprun,
                            contcar=args.contcar,
                            outcar=args.outcar,
                            procar=True,
                            defect_entry=de,
                            symprec=args.symprec,
                            angle_tolerance=args.angle_tolerance)
                except IOError:
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
                args.volume_dir, contcar_name=args.poscar)
        except IOError:
            raise FileNotFoundError(args.volume_dir, "not appropriate.")

    if args.total_dos_dir:
        try:
            dft_results.set_total_dos_from_vasp(args.total_dos_dir,
                                                vasprun_name=args.vasprun)
        except IOError:
            raise FileNotFoundError(args.total_dos_dir, "not appropriate.")

    dft_results.to_json_file(args.json_file)


def efnv_correction(args):
    if args.print:
        print(ExtendedFnvCorrection.load_json(args.json_file))
        return

    dirs = glob('*[0-9]/') if args.dir_all else args.dirs

    if args.nocorr:
        for directory in dirs:
            c = ManualCorrection(manual_correction_energy=args.manual)
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
            compute_correction(defect_entry=entry,
                               defect_dft=defect_dft_data,
                               perfect_dft=perfect_dft_data,
                               unitcell_dft=unitcell_dft_data,
                               ewald_json=args.read_ewald_json,
                               to_filename=args.dump_ewald_json)

        c.plot_distance_vs_potential(join(directory, "potential.pdf"))
        c.to_json_file(join(directory, "correction.json"))


def vasp_oba_set(args):

    #TODO: When writing GW part, refer oba_set_main.py in obadb

    flags = [str(s) for s in list(Element)]
    ldauu = list2dict(args.ldauu, flags)
    ldaul = list2dict(args.ldaul, flags)

    base_kwargs = {"task":                  args.task,
                   "xc":                    args.xc,
                   "kpt_density":           args.kpt_density,
                   "standardize_structure": args.standardize,
                   "ldauu": ldauu,
                   "ldaul": ldaul}

    flags = list(chain.from_iterable(incar_flags.values()))
    base_user_incar_settings = list2dict(args.incar_setting, flags)

    flags = list(signature(ObaSet.make_input).parameters.keys())
    base_kwargs.update(list2dict(args.vos_kwargs, flags))

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
                kwargs["band_gap"] = prior_info["band_gap"]
                kwargs["is_magnetization"] = True \
                    if abs(prior_info["total_magnetization"]) > 0.1 else False

        if args.prev_dir:
            files = {"CHGCAR": "C", "WAVECAR": "M", "WAVEDER": "M"}
            oba_set = ObaSet.from_prev_calc(args.prev_dir,
                                            charge=args.charge,
                                            copied_file_names=files, **kwargs)
        else:
            s = Structure.from_file(args.poscar)
            oba_set = \
                ObaSet.make_input(structure=s,
                                  charge=args.charge,
                                  user_incar_settings=user_incar_settings,
                                  weak_incar_settings={"LWAVE": args.wavecar},
                                  additional_user_potcar_yaml=args.potcar,
                                  **kwargs)

        oba_set.write_input(".")

    os.chdir(original_dir)


def defects(args):
    # try:
    #     unitcell = UnitcellCalcResults.load_json(args.unitcell)
    # except FileNotFoundError:
    #     print("{} not found".format(args.unitcell))

    try:
        perfect = SupercellCalcResults.load_json(args.perfect)
    except FileNotFoundError:
        print(f"{args.perfect} not found.")
        raise

    defects_dirs = args.defect_dirs if args.defect_dirs else glob('*[0-9]/')
    for d in defects_dirs:
        filename = join(d, args.json)
        if args.diagnose:
            print(d.rjust(12), end="  ")
            try:
                print(Defect.load_json(filename).diagnose)
            except FileNotFoundError:
                logger.warning("No supercell results file.")
            except Exception as e:
                logger.warning(f"An error {e} is caught.")
            continue

        logger.info("parsing directory {}...".format(d))
        files = [args.defect_entry, args.dft_results, args.correction]
        classes = [DefectEntry, SupercellCalcResults, ExtendedFnvCorrection]
        input_objects = generate_objects(d, files, classes, raise_error=False)

        if input_objects:
            defect = Defect.from_objects(defect_entry=input_objects[0],
                                         dft_results=input_objects[1],
                                         perfect_dft_results=perfect,
                                         correction=input_objects[2])
            defect.to_json_file(filename)
        else:
            logger.warning(f"Generating {filename} failed.")
            continue


def plot_energy(args):

    try:
        defect_energies = DefectEnergies.load_json(args.energies)
    except FileNotFoundError:
        unitcell = UnitcellCalcResults.load_json(args.unitcell)
        perfect = SupercellCalcResults.load_json(args.perfect)

        defects_dirs = args.dirs if args.dirs else glob('*[0-9]/')

        defect_list = []
        for d in defects_dirs:
            filename = join(d, args.defect)
            logger.info("parsing directory {}...".format(d))
            try:
                defect_list.append(Defect.load_json(filename))
            except FileNotFoundError:
                logger.warning(f"Parsing {filename} failed.")
                continue

        chem_pot = ChemPotDiag.load_vertices_yaml(args.chem_pot_yaml)

        # First construct DefectEnergies class object.
        defect_energies = DefectEnergies.from_objects(
            unitcell=unitcell, perfect=perfect, defects=defect_list,
            chem_pot=chem_pot, chem_pot_label=args.chem_pot_label,
            system=args.name)

    defect_energies.to_json_file(filename=args.energies)

    # if args.concentration:
    #     defect_concentration = \
    #         DefectConcentration.from_defect_energies(
    #             energies=energies,
    #             temperature=args.temperature[0],
    #             unitcell=unitcell,
    #             num_sites_filename=args.num_site_file)

    #     if len(args.temperature) == 2:
    #         defect_concentration = \
    #             DefectConcentration.from_defect_energies(
    #                 energies=energies,
    #                 temperature=args.temperature[1],
    #                 unitcell=unitcell,
    #                 num_sites_filename=args.num_site_file,
    #                 previous_concentration=defect_concentration)
    # else:
#        defect_concentration = None

    plt = defect_energies.plot_energy(filtering_words=args.filtering,
                                      x_range=args.x_range,
                                      y_range=args.y_range,
                                      show_transition_levels=args.show_tl,
                                      show_all_energies=args.show_all)

    plt.savefig(args.save_file, format="pdf") if args.save_file else plt.show()


def parse_eigenvalues(args):
    unitcell = UnitcellCalcResults.load_json(args.unitcell)

    # TODO: Modify here to run w/o correction
    logger.info("parsing directory {}...".format(args.defect_dir))
    filename = join(args.defect_dir, args.defect)
    defect = Defect.load_json(filename)

    defect_eigenvalues = DefectEigenvalue.from_files(unitcell=unitcell,
                                                     defect=defect)

    defect_eigenvalues.plot(yrange=args.y_range, title=args.title,
                            filename=args.save_file)


def vasp_parchg_set(args):
    user_incar_settings = {"LPARD": True, "LSEPB": True, "KPAR": 1,
                           "IBAND": args.band_indices}

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


def local_structure(args):
    defects_dirs = args.defect_dirs if args.defect_dirs else glob('*[0-9]/')

    for d in defects_dirs:
        logger.info("parsing directory {}...".format(d))
        filename = join(d, args.defect)
        logger.info("parsing directory {}...".format(d))
        try:
            defect = Defect.load_json(filename)
        except FileNotFoundError:
            logger.warning(f"Parsing {filename} failed.")
            continue

        defect_structure = DefectStructure.from_defect(defect)
        print(defect_structure.show_displacements(all_atoms=args.show_all))


def concentration(args):
    defect_energies = DefectEnergies.load_json(args.energies)
    unitcell = UnitcellCalcResults.load_json(args.unitcell)
    defect_concentration = DefectConcentration.from_calc_results(
        defect_energies=defect_energies,
        unitcell=unitcell,
        fractional_magnetization_to_one=args.fmto)

    defect_concentration.calc_equilibrium_concentration(
        temperature=args.temperature, verbose=args.verbose)
    defect_concentration.calc_quenched_equilibrium_concentration(
        verbose=args.verbose)
    print(defect_concentration)
    defect_concentration.calc_concentrations()
    plt = defect_concentration.plot_carrier_concentrations()
    plt.show()


# def make_refined_structure(args):
#     defects_dirs = args.defect_dirs if args.defect_dirs else glob('*[0-9]/')
#     original_dir = os.getcwd()

    # for d in defects_dirs:
    #     print(d.rjust(12), end="  ")
    #     os.chdir(join(original_dir, d))
    #     calc_results = SupercellCalcResults.load_json(args.dft_results)
    #     if not calc_results.is_converged:
    #         raise NoConvergenceError("Vasp is not converged.")

        # if calc_results.symmetrized_structure:
        #     os.mkdir(join(original_dir, d, "refined"))
        #     os.chdir(join(original_dir, d, "refined"))
        #     oba_set = ObaSet.from_prev_calc(dirname="..",
        #                                     parse_calc_results=False,
        #                                     parse_magnetization=False,
        #                                     standardize_structure=False,
        #                                     sort_structure=False,
        #                                     parse_potcar=True,
        #                                     parse_incar=True,
        #                                     parse_kpoints=True,
        #                                     copied_file_names={"WAVECAR": "L"})
        #     oba_set.write_input(".")
        #     poscar = Poscar(calc_results.symmetrized_structure)
        #     poscar.write_file("POSCAR")
        # else:
        #     raise ValueError("No symmetrized structure in DFT results.")


