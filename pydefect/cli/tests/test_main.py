# -*- coding: utf-8 -*-

from argparse import Namespace
import os
from pathlib import Path
import shutil

from pydefect.analysis.defect_structure import defect_structure_matcher
from pydefect.core.complex_defects import ComplexDefects
from pydefect.core.config import (
    DEFECT_SYMMETRY_TOLERANCE, SYMMETRY_TOLERANCE, ANGLE_TOL, DEFECT_KPT_DENSITY)
from pydefect.core.interstitial_site import (
    InterstitialSiteSet, interstitials_from_charge_density)
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.cli.main import simple_override, parse_args
from pydefect.core.defect_entry import DefectEntry
from pydefect.corrections.efnv_corrections import Ewald
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.input_maker.supercell_maker import Supercells
from pydefect.util.testing import PydefectTest

from vise.cli.main_tools import get_default_args

parent_dir = Path(__file__).parent


class SimpleOverrideTest(PydefectTest):
    def setUp(self) -> None:
        pydefect_yml = parent_dir / ".." / ".." / "test_files" / "pydefect.yaml"
        shutil.copy(str(pydefect_yml), "pydefect.yaml")

    def test(self):
        d = {"xc": "pbesol", "defect_symprec": 1}
        simple_override(d, ["defect_symprec"])
        expected = {'xc': 'pbesol', 'defect_symprec': 0.01}

        self.assertEqual(expected, d)

    def tearDown(self) -> None:
        os.remove("pydefect.yaml")


symprec_args = {"symprec": SYMMETRY_TOLERANCE,  "angle_tolerance": ANGLE_TOL}
defect_prec_args = {"defect_symprec": DEFECT_SYMMETRY_TOLERANCE,
                    "angle_tolerance": ANGLE_TOL}


class MainUnitcellCalcResultsTest(PydefectTest):

    def test_unitcell_calc_results_wo_options(self):
        actual = parse_args(["ur"])
#        d = get_default_args(UnitcellCalcResults)
        # func is a pointer so need to point the same address.
        expected = Namespace(
            json_file="unitcell.json",
            static_diele=None,
            ionic_diele=None,
            band_edge_dir=None,
            static_diele_dir=None,
            ionic_diele_dir=None,
            volume_dir=None,
            total_dos_dir=None,
            contcar="CONTCAR",
            outcar="OUTCAR",
            vasprun="vasprun.xml",
            print=False,
            func=actual.func)
        self.assertEqual(expected, actual)

    def test_unitcell_calc_results_w_options(self):
        actual = parse_args(["ur",
                             "--json_file", "a",
                             "--static_diele", "0.0", "0.1", "0.2",
                             "--ionic_diele", "0.3", "0.4", "0.5",
                             "--band_edge_dir", "b-dir",
                             "--static_diele_dir", "s-dir",
                             "--ionic_diele_dir", "i-dir",
                             "--volume_dir", "v-dir",
                             "--total_dos_dir", "d-dir",
                             "--contcar", "b",
                             "--outcar", "c",
                             "--vasprun", "d",
                             "--print"])
        expected = Namespace(
            json_file="a",
            static_diele=[0.0, 0.1, 0.2],
            ionic_diele=[0.3, 0.4, 0.5],
            band_edge_dir="b-dir",
            static_diele_dir="s-dir",
            ionic_diele_dir="i-dir",
            volume_dir="v-dir",
            total_dos_dir="d-dir",
            contcar="b",
            outcar="c",
            vasprun="d",
            print=True,
            func=actual.func)
        self.assertEqual(expected, actual)


class MainInitialSettingTest(PydefectTest):

    def test_initial_setting_wo_options(self):
        actual = parse_args(["is"])
        d = get_default_args(DefectInitialSetting.from_basic_settings)
        d.update(get_default_args(Supercells))
        # func is a pointer so need to point the same address.
        expected = Namespace(
            poscar="POSCAR",
            matrix=None,
            isotropy_criterion=d["criterion"],
            min_num_atoms=d["min_num_atoms"],
            max_num_atoms=d["max_num_atoms"],
            conventional_base=True,
            most_isotropic=False,
            rhombohedral_angle=d["rhombohedral_angle"],
            supercell_set=False,
            dopants=d["dopants"],
            antisite=True,
            en_diff=d["en_diff"],
            included=d["included"],
            excluded=d["excluded"],
            displacement_distance=d["displacement_distance"],
            interstitials=d["interstitial_sites"],
            complex_defect_names=d["complex_defect_names"],
            print_dopant=None,
            func=actual.func,
            **symprec_args)
        self.assertEqual(expected, actual)

    def test_initial_setting_w_options(self):
        actual = parse_args(["is",
                             "--poscar", "a",
                             "--matrix", "10",
                             "--criterion", "0.1",
                             "--min_num_atoms", "1",
                             "--max_num_atoms", "2",
                             "-cb", "F",
                             "--most_isotropic",
                             "--rhombohedral_angle", "3.1",
                             "--supercell_set",
                             "--dopants", "Ga", "In",
                             "--antisite", "F",
                             "--en_diff", "4.1",
                             "--included", "b",
                             "--excluded", "c",
                             "--displacement_distance", "5.1",
                             "--interstitial_sites", "i1", "i2",
                             "--complex_defect_names", "d", "e",
                             "--print_dopant", "Cs",
                             "--symprec", "6.1",
                             "--angle_tolerance", "7.1"])
        expected = Namespace(
            poscar="a",
            matrix=[10],
            isotropy_criterion=0.1,
            min_num_atoms=1,
            max_num_atoms=2,
            conventional_base=False,
            most_isotropic=True,
            rhombohedral_angle=3.1,
            supercell_set=True,
            dopants=["Ga", "In"],
            antisite=False,
            en_diff=4.1,
            included=["b"],
            excluded=["c"],
            displacement_distance=5.1,
            interstitials=["i1", "i2"],
            complex_defect_names=["d", "e"],
            print_dopant="Cs",
            symprec=6.1,
            angle_tolerance=7.1,
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainInterstitialTest(PydefectTest):

    def test_interstitial_wo_options(self):
        actual = parse_args(["i"])
        d = get_default_args(InterstitialSiteSet.add_sites)
        d.update(get_default_args(InterstitialSiteSet.from_files))
        d.update(get_default_args(interstitials_from_charge_density))
        # func is a pointer so need to point the same address.
        expected = Namespace(
            yaml=d["yaml_filename"],
            dposcar=d["dposcar"],
            defect_in="defect.in",
            interstitial_coords=None,
            site_name=None,
            radius=d["vicinage_radius"],
            force_add=False,
            interstitial_symprec=d["interstitial_symprec"],
            method=d["method"],
            chgcar=None,
            func=actual.func,
            **defect_prec_args)
        self.assertEqual(expected, actual)

    def test_interstitial_w_options(self):
        actual = parse_args(["i",
                             "--yaml", "a",
                             "--dposcar", "b",
                             "--defect_in", "c",
                             "-c", "0.1", "0.2", "0.3",
                             "--name", "d",
                             "--radius", "1.1",
                             "--force_add",
                             "--interstitial_symprec", "2.1",
                             "--method", "e",
                             "--chgcar", "f",
                             "--defect_symprec", "3.1",
                             "--angle_tolerance", "4.1",
                             ])
        expected = Namespace(
            yaml="a",
            dposcar="b",
            defect_in="c",
            interstitial_coords=[0.1, 0.2, 0.3],
            site_name="d",
            radius=1.1,
            force_add=True,
            interstitial_symprec=2.1,
            method="e",
            chgcar="f",
            defect_symprec=3.1,
            angle_tolerance=4.1,
            func=actual.func,
            )
        self.assertEqual(expected, actual)


class MainComplexDefectsTest(PydefectTest):

    def test_complex_defects_wo_options(self):
        actual = parse_args(["cd"])
        d = get_default_args(ComplexDefects.add_defect)
        d.update(get_default_args(ComplexDefects.from_files))
        # func is a pointer so need to point the same address.
        expected = Namespace(
            yaml=d["yaml_filename"],
            dposcar=d["dposcar"],
            defect_in="defect.in",
            removed_atom_indices=None,
            inserted_elements=None,
            inserted_coords=None,
            name=None,
            extreme_charge_state=None,
            annotation=None,
            func=actual.func,
            **defect_prec_args,
        )
        self.assertEqual(expected, actual)

    def test_complex_defects_w_options(self):
        actual = parse_args(["cd",
                             "--yaml", "a",
                             "--dposcar", "b",
                             "--defect_in", "c",
                             "--removed_atom_indices", "1",
                             "--inserted_elements", "H",
                             "--inserted_coords", "0.1", "0.2",
                             "--name", "d",
                             "--extreme_charge_state", "100",
                             "--annotation", "e",
                             "--defect_symprec", "1.1",
                             "--angle_tolerance", "2.1",
                             ])
        d = get_default_args(ComplexDefects.add_defect)
        d.update(get_default_args(ComplexDefects.from_files))
        # func is a pointer so need to point the same address.
        expected = Namespace(
            yaml="a",
            dposcar="b",
            defect_in="c",
            removed_atom_indices=[1],
            inserted_elements=["H"],
            inserted_coords=[0.1, 0.2],
            name="d",
            extreme_charge_state=100,
            annotation="e",
            defect_symprec=1.1,
            angle_tolerance=2.1,
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainDefectVaspSetMakerTest(PydefectTest):
    def test_defect_vasp_set_maker_wo_options(self):
        actual = parse_args(["dvs"])
        # default set used in the vise unittest
        from vise.cli.tests.test_main import default_vasp_args
        default_vasp_args.pop("charge")
        default_vasp_args.pop("task")
        expected = Namespace(
            spin_polarize=True,
            defect_in="defect.in",
            dposcar="DPOSCAR",
            keywords=None,
            specified_defects=None,
            force_overwrite=False,
            kpt_density=DEFECT_KPT_DENSITY,
            wavecar=True,
            func=actual.func,
            **default_vasp_args,
        )
        self.assertEqual(expected, actual)

    def test_defect_vasp_set_maker_w_options(self):
        actual = parse_args(["dvs",
                             "-s", "F",
                             "--potcar", "Mg_pv", "O_h",
                             "--potcar_set_name", "gw",
                             "-x", "pbesol",
                             "--vise_opts", "encut", "800",
                             "--user_incar_settings", "LREAD", "F",
                             "-auis", "ALGO", "D",
                             "--ldauu", "Mg", "5",
                             "--ldaul", "Mg", "1",

                             "--defect_in", "a",
                             "--dposcar", "b",
                             "-kw", "c", "d",
                             "-d", "e", "f",
                             "--force_overwrite",
                             "-k", "0.1",
                             "-w", "F",
                             ])
        expected = Namespace(
            spin_polarize=False,
            potcar_set=["Mg_pv", "O_h"],
            potcar_set_name="gw",
            xc="pbesol",
            vise_opts=["encut", "800"],
            user_incar_settings=["LREAD", "F"],
            additional_user_incar_settings=["ALGO", "D"],
            ldauu=["Mg", "5"],
            ldaul=["Mg", "1"],
            defect_in="a",
            dposcar="b",
            keywords=["c", "d"],
            specified_defects=["e", "f"],
            force_overwrite=True,
            kpt_density=0.1,
            wavecar=False,
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainDefectEntryTest(PydefectTest):
    def test_defect_entry_wo_options(self):
        actual = parse_args(["de"])
        d = get_default_args(DefectEntry.from_defect_structure)
        expected = Namespace(
            print=False,
            make_defect_entry=False,
            yaml=None,
            json="defect_entry.json",
            perfect_poscar="../perfect/POSCAR",
            defect_poscar="POSCAR",
            displacement_distance=d["displacement_distance"],
            defect_name=d["defect_name"],
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_defect_entry_w_options(self):
        actual = parse_args(["de",
                             "--print",
                             "--make_defect_entry",
                             "--yaml", "a",
                             "--json", "b",
                             "--perfect_poscar", "c",
                             "--defect_poscar", "d",
                             "--displacement_distance", "0.1",
                             "--defect_name", "e",
                             ])
        expected = Namespace(
            print=True,
            make_defect_entry=True,
            yaml="a",
            json="b",
            perfect_poscar="c",
            defect_poscar="d",
            displacement_distance=0.1,
            defect_name="e",
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainSupercellCalcResultsTest(PydefectTest):
    def test_supercell_calc_results_wo_options(self):
        actual = parse_args(["sr"])
        d = get_default_args(SupercellCalcResults.from_vasp_files)
        expected = Namespace(
            dirs=None,
            dir_all=False,
            vasprun=d["vasprun"],
            contcar=d["contcar"],
            outcar=d["outcar"],
            procar=d["procar"],
            defect_center=None,
            defect_entry_name="defect_entry.json",
            json="dft_results.json",
            cutoff=d["cutoff"],
            print=False,
            func=actual.func,
            **defect_prec_args
        )
        self.assertEqual(expected, actual)

    def test_supercell_calc_results_w_options(self):
        actual = parse_args(["sr",
                             "--dirs", "a", "b",
                             "--dir_all",
                             "-v", "c",
                             "-c", "d",
                             "-o", "e",
                             "-p", "f",
                             "--center", "10",
                             "-de", "g",
                             "--json", "h",
                             "--cutoff", "0.1",
                             "--print",
                             "--defect_symprec", "1.1",
                             "--angle_tolerance", "2.1",
                             ])
        expected = Namespace(
            dirs=["a", "b"],
            dir_all=True,
            vasprun="c",
            contcar="d",
            outcar="e",
            procar="f",
            defect_center=[10.0], # be careful
            defect_entry_name="g",
            json="h",
            cutoff=0.1,
            print=True,
            defect_symprec=1.1,
            angle_tolerance=2.1,
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainExtendedFNVCorrectionTest(PydefectTest):
    def test_supercell_calc_results_wo_options(self):
        actual = parse_args(["efc"])
        d = get_default_args(Ewald.from_optimization)
        expected = Namespace(
            nocorr=False,
            manual=0.0,
            unitcell_json="../unitcell/unitcell.json",
            perfect_json="perfect/dft_results.json",
            ewald_json="ewald.json",
            ewald_initial_param=d["initial_ewald_param"],
            ewald_convergence=d["convergence"],
            ewald_accuracy=d["prod_cutoff_fwhm"],
            dirs=None,
            dir_all=False,
            force_overwrite=False,
            defect_center=None,
            json_file="correction.json",
            print=False,
            plot_potential=False,
            y_range=None,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_supercell_calc_results_w_options(self):
        actual = parse_args(["efc",
                             "--nocorr",
                             "--manual", "0.1",
                             "--unitcell_json", "a",
                             "--perfect_json", "b",
                             "--ewald_json", "c",
                             "--ewald_initial_param", "0.2",
                             "--ewald_convergence", "0.3",
                             "--ewald_accuracy", "0.4",
                             "--dirs", "d", "e",
                             "--dir_all",
                             "--force_overwrite",
                             "--center", "0.5", "0.6", "0.7",
                             "--json_file", "f",
                             "--print",
                             "--plot_potential",
                             "--y_range", "0.8", "0.9",
                             ])
        expected = Namespace(
            nocorr=True,
            manual=0.1,
            unitcell_json="a",
            perfect_json="b",
            ewald_json="c",
            ewald_initial_param=0.2,
            ewald_convergence=0.3,
            ewald_accuracy=0.4,
            dirs=["d", "e"],
            dir_all=True,
            force_overwrite=True,
            defect_center=[0.5, 0.6, 0.7],
            json_file="f",
            print=True,
            plot_potential=True,
            y_range=[0.8, 0.9],
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainVerticalTransitionInputMakerTest(PydefectTest):
    def test_vertical_transition_input_maker_wo_options(self):
        actual = parse_args(["vtim"])
        expected = Namespace(
            additional_charge=None,
            initial_dir_name=None,
            contcar="CONTCAR",
            func=actual.func,
            )
        self.assertEqual(expected, actual)

    def test_vertical_transition_input_maker_w_options(self):
        actual = parse_args(["vtim",
                             "--additional_charge", "1",
                             "--initial_dir_name", "a",
                             "--contcar", "b",
                             ])
        expected = Namespace(
            additional_charge=1,
            initial_dir_name="a",
            contcar="b",
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainVerticalTransitionEnergyTest(PydefectTest):
    def test_vertical_transition_energy_wo_options(self):
        actual = parse_args(["vte"])
        expected = Namespace(
            unitcell_json="../unitcell/unitcell.json",
            initial_dir=None,
            dir=None,
            print=False,
            y_range=None,
            json="vte_correction.json",
            show=False,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_vertical_transition_energy_w_options(self):
        actual = parse_args(["vte",
                             "--unitcell_json", "a",
                             "--initial_dir", "b",
                             "--dir", "c",
                             "--print",
                             "--y_range", "0.1", "0.2",
                             "--json", "d",
                             "--show",
                             ])
        expected = Namespace(
            unitcell_json="a",
            initial_dir="b",
            dir="c",
            print=True,
            y_range=[0.1, 0.2],
            json="d",
            show=True,
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainDefectTest(PydefectTest):
    def test_defect_wo_options(self):
        actual = parse_args(["d"])
        expected = Namespace(
            defect_dirs=None,
            diagnose=False,
            json="defect.json",
            perfect="perfect/dft_results.json",
            defect_entry="defect_entry.json",
            correction="correction.json",
            dft_results="dft_results.json",
            band_edge=None,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_defect_w_options(self):
        actual = parse_args(["d",
                             "--defect_dirs", "a", "b",
                             "--diagnose",
                             "--json", "c",
                             "--perfect", "d",
                             "--defect_entry", "e",
                             "--correction", "f",
                             "--dft_results", "g",
                             "--band_edge", "h", "i",
                             ])
        expected = Namespace(
            defect_dirs=["a", "b"],
            diagnose=True,
            json="c",
            perfect="d",
            defect_entry="e",
            correction="f",
            dft_results="g",
            band_edge=["h", "i"],
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainPlotEnergyTest(PydefectTest):
    def test_plot_energy_wo_options(self):
        actual = parse_args(["pe"])
        expected = Namespace(
            name="",
            x_range=None,
            y_range=None,
            save_file=None,
            energies="defect_energies.json",
            unitcell="../unitcell/unitcell.json",
            perfect="perfect/dft_results.json",
            defect="defect.json",
            defect_dirs=None,
            chem_pot_json="../competing_phases/cpd.json",
            chem_pot_label="A",
            filtering=None,
            concentration=False,
            show_transition_level=False,
            show_all=False,
            reload_defects=False,
            print=False,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_plot_energy_w_options(self):
        actual = parse_args(["pe",
                             "--name", "a",
                             "--x_range", "0.1", "0.2",
                             "--y_range", "0.3", "0.4",
                             "--save_file", "b",
                             "--energies", "c",
                             "--unitcell", "d",
                             "--perfect", "e",
                             "--defect", "f",
                             "--defect_dirs", "g", "h",
                             "--chem_pot_json", "i",
                             "--chem_pot_label", "j",
                             "--filtering", "k",
                             "--concentration",
                             "--show_transition_level",
                             "--show_all",
                             "--reload_defects",
                             "--print",
                             ])
        expected = Namespace(
            name="a",
            x_range=[0.1, 0.2],
            y_range=[0.3, 0.4],
            save_file="b",
            energies="c",
            unitcell="d",
            perfect="e",
            defect="f",
            defect_dirs=["g", "h"],
            chem_pot_json="i",
            chem_pot_label="j",
            filtering="k",
            concentration=True,
            show_transition_level=True,
            show_all=True,
            reload_defects=True,
            print=True,
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainParseEigenvaluesTest(PydefectTest):
    def test_parse_eigenvalues_wo_options(self):
        actual = parse_args(["eig"])
        expected = Namespace(
            title="",
            y_range=None,
            save_file=None,
            unitcell="../unitcell/unitcell.json",
            defect="defect.json",
            defect_dir=None,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_parse_eigenvalues_w_options(self):
        actual = parse_args(["eig",
                             "--title", "a",
                             "--y_range", "0.1", "0.2",
                             "--save_file", "b",
                             "--unitcell", "c",
                             "--defect", "d",
                             "--defect_dir", "e",
                             ])
        expected = Namespace(
            title="a",
            y_range=[0.1, 0.2],
            save_file="b",
            unitcell="c",
            defect="d",
            defect_dir="e",
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainVaspParchgSetTest(PydefectTest):
    def test_vasp_parchg_set_wo_options(self):
        actual = parse_args(["vps"])
        expected = Namespace(
            read_dir=None,
            write_dir=".",
            band_indices=None,
            kpoint_indices=None,
            contcar="CONTCAR",
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_vasp_parchg_set_w_options(self):
        actual = parse_args(["vps",
                             "--read_dir", "a",
                             "--write_dir", "b",
                             "--band_indices", "1", "2", "3",
                             "--kpoint_indices", "4", "5", "6",
                             "--contcar", "c",
                             ])
        expected = Namespace(
            read_dir="a",
            write_dir="b",
            band_indices=[1, 2, 3],
            kpoint_indices=[4, 5, 6],
            contcar="c",
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainLocalStructureTest(PydefectTest):
    def test_local_structure_wo_options(self):
        actual = parse_args(["ls"])
        d = get_default_args(defect_structure_matcher)

        expected = Namespace(
            defect="defect.json",
            show_all=False,
            compare_structure=False,
            site_tolerance=d["site_tolerance"],
            defect_dirs=None,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_local_structure_w_options(self):
        actual = parse_args(["ls",
                             "--defect", "a",
                             "--show_all",
                             "--compare_structure",
                             "--site_tolerance", "0.1",
                             "--defect_dirs", "b", "c"
                             ])
        expected = Namespace(
            defect="a",
            show_all=True,
            compare_structure=True,
            site_tolerance=0.1,
            defect_dirs=["b", "c"],
            func=actual.func,
        )
        self.assertEqual(expected, actual)


class MainConcentrationsTest(PydefectTest):
    def test_concentrations_wo_options(self):
        actual = parse_args(["c"])
        expected = Namespace(
            energies="defect_energies.json",
            unitcell="../unitcell/unitcell.json",
            filtering=None,
            frac_mag_to_one=False,
            verbose=False,
            temperature=None,
            func=actual.func,
        )
        self.assertEqual(expected, actual)

    def test_concentrations_w_options(self):
        actual = parse_args(["c",
                             "--energies", "a",
                             "--unitcell", "b",
                             "--filtering", "c",
                             "--frac_mag_to_one",
                             "--verbose",
                             "--temperature", "1.1", "2.1",
                             ])
        expected = Namespace(
            energies="a",
            unitcell="b",
            filtering="c",
            frac_mag_to_one=True,
            verbose=True,
            temperature=[1.1, 2.1],
            func=actual.func,
        )
        self.assertEqual(expected, actual)
