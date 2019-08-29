# -*- coding: utf-8 -*-

import os
import tempfile

from pydefect.analysis.defect_carrier_concentration import *
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.testing import PydefectTest
from scipy.constants import physical_constants

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")

EV = physical_constants['Boltzmann constant in eV/K'][0]


class AllCombinationTest(PydefectTest):
    def test(self):
        # a = {'Va_O1': {0: {None: {'energy': 1.362400134999998, 'convergence': True, 'is_shallow': False}}, 1: {None: {'energy': -4.2058905866258645, 'convergence': True, 'is_shallow': False}}, 2: {None: {'energy': -7.983747657760918, 'convergence': True, 'is_shallow': False}, 'inward': {'energy': -6.016447897493341, 'convergence': True, 'is_shallow': False}}}}
        # print(all_combination(a))
        b = {0: {None: {'energy': 1.362400134999998, 'convergence': True, 'is_shallow': False}}, 1: {None: {'energy': -4.2058905866258645, 'convergence': True, 'is_shallow': False}}, 2: {None: {'energy': -7.983747657760918, 'convergence': True, 'is_shallow': False}, 'inward': {'energy': -6.016447897493341, 'convergence': True, 'is_shallow': False}}}
        # print(all_combination(a))
        print(flatten_dict(b))


class CalcCarrierConcentrationTest(PydefectTest):
    def setUp(self):
        # unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        # unitcell = UnitcellCalcResults.load_json(unitcell_file)

        self.vbm = 0
        self.cbm = 10
        energies = np.linspace(-1, 11, num=121)
        # Dos is halved at the boundaries.
        # dos -> 0.5 at -0.2, 0, 10, and 10.2, and 1 at -0.1 and 10.1
        dos = [0.0] * 8 + [0.5] + [1.0] + [0.5] + [0.0] * 99 + [0.5] + [1.0] + \
              [0.5] + [0.0] * 8
        self.total_dos = np.array([dos, energies])
        self.volume = 100
        self.temperature = 10000
        # [3.0036, 7.4934]
        # self.vbm, self.cbm = unitcell.band_edge
        # self.volume = unitcell.volume

    def test_hole_concentration(self):
        actual = hole_concentration(temperature=self.temperature, e_f=4,
                                    total_dos=self.total_dos,
                                    vbm=self.vbm, volume=self.volume)
        fd_part = [0.5 / (np.exp((4 - 0) / (10000 * 8.6173303e-05)) + 1),
                   1.0 / (np.exp((4 - (-0.1)) / (10000 * 8.6173303e-05)) + 1),
                   0.5 / (np.exp((4 - (-0.2)) / (10000 * 8.6173303e-05)) + 1)]
        expected = sum(fd_part)
        expected *= 0.1 / (100 / 10 ** 24)
        # expected = 1.707781233296676e+19
        self.assertAlmostEqual(actual / expected, 1, places=6)

    def test_electron_concentration(self):
        actual = electron_concentration(temperature=self.temperature, e_f=4,
                                        total_dos=self.total_dos,
                                        cbm=self.cbm, volume=self.volume)
        fd_part = [0.5 / (np.exp((10.0 - 4) / (10000 * 8.6173303e-05)) + 1),
                   1.0 / (np.exp((10.1 - 4) / (10000 * 8.6173303e-05)) + 1),
                   0.5 / (np.exp((10.2 - 4) / (10000 * 8.6173303e-05)) + 1)]
        expected = sum(fd_part)
        expected *= 0.1 / (100 / 10 ** 24)
        # expected = 1.6898803452706716e+18
        self.assertAlmostEqual(actual / expected, 1, places=6)


class CalcConcentrationTest(PydefectTest):
    def setUp(self):
        self.defect_energies = defaultdict(dict)
        self.defect_energies["Va_O1"][0] = {None: 4.0}
        self.defect_energies["Va_O1"][1] = {None: 3.0}
        self.defect_energies["Va_O1"][2] = {None: 1.0}
#        self.energies["Va_Mg1"][0] = {None: 4.0}
#        self.energies["Va_Mg1"][-1] = {None: 5.0}
        self.defect_energies["Va_Mg1"][-2] = {None: 6.0}

        self.vbm = 0
        self.cbm = 10
        energies = np.linspace(-1, 11, num=121)
        # Dos is halved at the boundaries.
        # dos -> 0.5 at -0.2, 0, 10, and 10.2, and 1 at -0.1 and 10.1
        dos = [0.0] * 8 + [0.5] + [1.0] + [0.5] + [0.0] * 99 + [0.5] + [1.0] + \
              [0.5] + [0.0] * 8
        self.total_dos = np.array([dos, energies])

        self.multiplicity = defaultdict(dict)
        self.multiplicity["Va_O1"][0] = {None: 1}
        self.multiplicity["Va_O1"][1] = {None: 1}
        self.multiplicity["Va_O1"][2] = {None: 2}
#        self.multiplicity["Va_Mg1"][0] = {None: 1}
#        self.multiplicity["Va_Mg1"][-1] = {None: 1}
        self.multiplicity["Va_Mg1"][-2] = {None: 1}

        self.magnetization = defaultdict(dict)
        self.magnetization["Va_O1"][0] = {None: 0}
        self.magnetization["Va_O1"][1] = {None: 1}
        self.magnetization["Va_O1"][2] = {None: 1}
#        self.magnetization["Va_Mg1"][0] = {None: 2}
#        self.magnetization["Va_Mg1"][-1] = {None: 1}
        self.magnetization["Va_Mg1"][-2] = {None: 0}

        self.temperature = 10000
        self.volume = 100

    def test(self):
        concentration = \
            calc_concentration(defect_energies=self.defect_energies,
                               temperature=self.temperature,
                               e_f=4,
                               vbm=self.vbm,
                               cbm=self.cbm,
                               total_dos=self.total_dos,
                               volume=self.volume)

        va_o1_2_expected = np.exp(-9 / (10000 * 8.6173303e-05)) * 2 * 2
        va_o1_2_expected /= (100 / 10 ** 24)
        va_mg1_m2_expected = np.exp(-(6.0 - 2*4) / (10000 * 8.6173303e-05)) * 1 * 1
        va_mg1_m2_expected /= (100 / 10 ** 24)

        self.assertAlmostEqual(
            concentration["Va_O1"][2][None] / va_o1_2_expected, 1, places=5)
        self.assertAlmostEqual(
            concentration["Va_Mg1"][-2][None] / va_mg1_m2_expected, 1, places=5)

    def test_ref_concentration(self):
        ref = defaultdict(dict)
        #        self.energies["Va_O1"][0] = {None: 4.0}
        #        self.energies["Va_O1"][1] = {None: 3.0}
        ref["Va_O1"][2] = {None: 1e+20}
        ref["Va_O1"][1] = {None: 1.5e+20}
        ref["Va_O1"][0] = {None: 0.5e+20}
        ref["Va_Mg1"][2] = {None: 1e+20}
        #        self.energies["Va_Mg1"][0] = {None: 4.0}

        concentration = \
            calc_concentration(defect_energies=self.defect_energies,
                               temperature=self.temperature,
                               e_f=4,
                               vbm=self.vbm,
                               cbm=self.cbm,
                               total_dos=self.total_dos,
                               volume=self.volume,
                               ref_concentration=ref)

        print(concentration)


    def test2(self):
        equiv_concentration = \
            calc_equilibrium_concentration(defect_energies=self.defect_energies,
                                           temperature=self.temperature,
                                           vbm=self.vbm,
                                           cbm=self.cbm,
                                           total_dos=self.total_dos,
                                           volume=self.volume,
                                           verbose=True)

        print(equiv_concentration[1])


class DefectConcentrationTest(PydefectTest):

    def setUp(self):
        """ """
        file = os.path.join(test_dir, "MgO/defects/defect_energies.json")
        defect_energies = DefectEnergies.load_json(file)
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        unitcell = UnitcellCalcResults.load_json(unitcell_file)

        self.defect_concentration = \
            DefectConcentration.from_calc_results(
                defect_energies=defect_energies, unitcell=unitcell)

    def test(self):
        # print(self.defect_concentration.energies)
        # print(self.defect_concentration.multiplicity)
        # print(self.defect_concentration.magnetization)
        # print(self.defect_concentration.volume)
        # print(self.defect_concentration.vbm)
        # print(self.defect_concentration.cbm)
        # print(self.defect_concentration.total_dos)

        # print(self.defect_concentration.temperature)
        # print(self.defect_concentration.concentrations)

        self.defect_concentration.calc_equilibrium_concentration(temperature=1000, verbose=False)
        print(self.defect_concentration.equilibrium_concentration)
        self.defect_concentration.calc_quenched_equilibrium_concentration(temperature=298, verbose=False)
        print(self.defect_concentration.quenched_equilibrium_concentration)

        self.defect_concentration.calc_concentrations(temperature=1000)
        print(self.defect_concentration)
#        self.defect_concentration.plot_carrier_concentrations()
#        print(self.defect_concentration.concentrations)
#         print(self.defect_concentration.equilibrium_concentration)
#         print("-----------------------------------------")
#        print(self.defect_concentration.concentrations)

    # def test_from_defect_energies(self):
    #     temperature = 10000
    #     num_sites_filename = os.path.join(test_dir,
    #                                       "MgO/defects/num_sites.yaml")

        # dc1 = DefectConcentration.from_defect_energies(
        #     energies=self.energies,
        #     temperature=temperature,
        #     unitcell=self.unitcell,
        #     num_sites_filename=num_sites_filename)

        # print(dc1.energies)
        # print(dc1.temperature)
        # print(dc1.e_f)
        # print(p)
        # print(n)
        # print(dc1.concentration)

        # temperature2 = 1000

#         # dc2 = DefectConcentration.from_defect_energies(
#         #     energies=self.energies,
#         #     temperature=temperature2,
#         #     unitcell=self.unitcell,
#         #     num_sites_filename=num_sites_filename,
#         #     previous_concentration=dc1,
#         #     verbose=True)

#         # print("-------------------------------")
#         # print(dc2.energies)
#         # print(dc2.temperature)
#         # print(dc2.e_f)
#         # print(dc2.p)
#         # print(dc2.n)
#         # print(dc2.concentration)

    def test_dict(self):
        """ round trip test of to_json and from_json """
        self.defect_concentration.calc_equilibrium_concentration(temperature=1000, verbose=False)
        self.defect_concentration.calc_quenched_equilibrium_concentration(temperature=298, verbose=False)
        self.defect_concentration.calc_concentrations(temperature=1000)

        d = self.defect_concentration.as_dict()
        defect_concentration_from_dict = DefectConcentration.from_dict(d)
        self.assertEqual(defect_concentration_from_dict.as_dict(), d)

    def test_json(self):
        """ round trip test of to_json and from_json """

        tmp_file = tempfile.NamedTemporaryFile().name
        self.defect_concentration.to_json_file(tmp_file)
        defect_concentration_from_json = DefectConcentration.load_json(tmp_file)
        print("")
        print(defect_concentration_from_json.as_dict())
        print(self.defect_concentration.as_dict())
        self.assertEqual(defect_concentration_from_json.as_dict(),
                         self.defect_concentration.as_dict())

