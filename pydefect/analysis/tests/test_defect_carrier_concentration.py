# -*- coding: utf-8 -*-

import os
import numpy as np

from pydefect.analysis.defect_carrier_concentration import (
    hole_concentration, electron_concentration, calc_concentration,
    calc_equilibrium_concentration, DefectConcentration)
from pydefect.analysis.defect_energies import DefectEnergy
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.testing import PydefectTest


class CalcCarrierConcentrationTest(PydefectTest):
    def setUp(self):
        self.vbm = 0
        self.cbm = 10
        energies = np.linspace(-1, 11, num=121).tolist()
        # Dos is halved at the boundaries.
        # dos -> 0.5 at -0.2, 0, 10, and 10.2, and 1 at -0.1 and 10.1
        dos = [0.0] * 8 + [0.5] + [1.0] + [0.5] + [0.0] * 99 + [0.5] + [1.0] + \
              [0.5] + [0.0] * 8
        self.total_dos = [dos, energies]
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
        self.assertAlmostEqual(actual / expected, 1, places=4)

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
        self.assertAlmostEqual(actual / expected, 1, places=4)


class CalcConcentrationTest(PydefectTest):
    def setUp(self):
        self.defect_energies = {"Va_O1": {}, "Va_Mg1": {}}

        self.defect_energies["Va_O1"][0] = \
            DefectEnergy(defect_energy=4.0,
                         annotation=None,
                         multiplicity=1,
                         magnetization=0.0,
                         convergence=True,
                         shallow=False)
        self.defect_energies["Va_O1"][1] = \
            DefectEnergy(3.0, None, 1, 1.0, True, False)
        self.defect_energies["Va_O1"][2] = \
            DefectEnergy(1.0, None, 2, 1.0, True, False)
        self.defect_energies["Va_Mg1"][-2] = \
            DefectEnergy(4.0, None, 1, 0.0, True, False)

        self.vbm = 0
        self.cbm = 10
        energies = np.linspace(-1, 11, num=121).tolist()
        # DOS is halved at the boundaries.
        # dos -> 0.5 at -0.2, 0, 10, and 10.2, and 1 at -0.1 and 10.1
        dos = [0.0] * 8 + [0.5] + [1.0] + [0.5] + [0.0] * 99 + [0.5] + [1.0] + \
              [0.5] + [0.0] * 8
        self.total_dos = [dos, energies]

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

        # Boltzmann distribution
        energy = 1 + 4 * 2
        t = 10000
        t_to_ev = 8.6173303e-05
        multiplicity = 2
        spin_deg = 2
        va_o1_2_expected = \
            np.exp(-energy / (t * t_to_ev)) * multiplicity * spin_deg
        va_o1_2_expected /= (100 / 10 ** 24)

        energy = 4 - 4 * 2
        va_mg1_m2_expected = \
            np.exp(-energy / (t * t_to_ev)) * 1 * 1
        va_mg1_m2_expected /= (100 / 10 ** 24)

        print(concentration["Va_O1"])

        self.assertAlmostEqual(
            concentration["Va_O1"][2] / va_o1_2_expected, 1, places=5)
        self.assertAlmostEqual(
            concentration["Va_Mg1"][-2] / va_mg1_m2_expected, 1, places=5)

    def test_ref_concentration(self):
        ref = {"Va_O1": {2: 1e+20, 1: 1.5e+20, 0: 0.5e+20},
               "Va_Mg1": {2: 1e+20}}

        concentration = \
            calc_concentration(defect_energies=self.defect_energies,
                               temperature=self.temperature,
                               e_f=4,
                               vbm=self.vbm,
                               cbm=self.cbm,
                               total_dos=self.total_dos,
                               volume=self.volume,
                               ref_concentration=ref)

        va_o_expected = 3.0e+20
        actual = (concentration["Va_O1"][0]
                  + concentration["Va_O1"][1]
                  + concentration["Va_O1"][2])
        self.assertAlmostEqual(actual / va_o_expected, 1, places=5)

        va_o_0_expected = (va_o_expected * 96.402441 /
                           (96.402441 + 5.931767 + 1.164818))
        actual = concentration["Va_O1"][0]
        self.assertAlmostEqual(actual / va_o_0_expected, 1, places=5)


class CalcEquilibriumConcentrationTest(PydefectTest):
    def setUp(self):
        self.defect_energies = {"Va_O1": {}, "Va_Mg1": {}}

        self.defect_energies["Va_O1"][0] = \
            DefectEnergy(defect_energy=4.0,
                         annotation=None,
                         multiplicity=1,
                         magnetization=0.0,
                         convergence=True,
                         shallow=False)
        self.defect_energies["Va_O1"][1] = \
            DefectEnergy(3.0, None, 1, 1.0, True, False)
        self.defect_energies["Va_O1"][2] = \
            DefectEnergy(1.0, None, 2, 1.0, True, False)
        self.defect_energies["Va_Mg1"][-2] = \
            DefectEnergy(4.0, None, 1, 0.0, True, False)

        self.vbm = 0
        self.cbm = 10
        energies = np.linspace(-1, 11, num=121).tolist()
        # DOS is halved at the boundaries.
        # dos -> 0.5 at -0.2, 0, 10, and 10.2, and 1 at -0.1 and 10.1
        dos = [0.0] * 8 + [0.5] + [1.0] + [0.5] + [0.0] * 99 + [0.5] + [1.0] + \
              [0.5] + [0.0] * 8
        self.total_dos = [dos, energies]

        self.temperature = 10000
        self.volume = 100

    def test(self):
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
        self.defect_energies = {"Va_O1": {}, "Va_Mg1": {}}

        self.defect_energies["Va_O1"][0] = \
            DefectEnergy(defect_energy=4.0,
                         annotation=None,
                         multiplicity=1,
                         magnetization=0.0,
                         convergence=True,
                         shallow=False)
        self.defect_energies["Va_O1"][1] = \
            DefectEnergy(3.0, None, 1, 1.0, True, False)
        self.defect_energies["Va_O1"][2] = \
            DefectEnergy(1.0, None, 2, 1.0, True, False)
        self.defect_energies["Va_Mg1"][-2] = \
            DefectEnergy(4.0, None, 1, 0.0, True, False)
        energies = np.linspace(-1, 11, num=121).tolist()
        # Dos is halved at the boundaries.
        # dos -> 0.5 at -0.2, 0, 10, and 10.2, and 1 at -0.1 and 10.1
        dos = [0.0] * 8 + [0.5] + [1.0] + [0.5] + [0.0] * 99 + [0.5] + [1.0] + \
              [0.5] + [0.0] * 8
        total_dos = [dos, energies]

        self.unitcell = UnitcellCalcResults(band_edge=[0, 10],
                                            static_dielectric_tensor=None,
                                            ionic_dielectric_tensor=None,
                                            total_dos=total_dos,
                                            volume=100)

        self.unitcell_no_band_edge = \
            UnitcellCalcResults(band_edge=None,
                                static_dielectric_tensor=None,
                                ionic_dielectric_tensor=None,
                                total_dos=total_dos,
                                volume=100)

        self.defect_concentration = \
            DefectConcentration.from_calc_results(self.defect_energies,
                                                  self.unitcell)

        self.defect_concentration.calc_equilibrium_concentration(
            temperature=10000, verbose=False)
        self.defect_concentration.calc_quenched_equilibrium_concentration(
            temperature=298, verbose=False)

    def test_unset_band_edge(self):
        with self.assertRaises(ValueError):
            DefectConcentration.from_calc_results(self.defect_energies,
                                                  self.unitcell_no_band_edge)

    def test_concentration(self):
        expected = 9.663745575119449e+20
        actual = self.defect_concentration.equilibrium_concentration["Va_O1"][2]
        self.assertAlmostEqual(expected, actual)

    def test_dict(self):
        expected = self.defect_concentration.as_dict()
        actual = DefectConcentration.from_dict(expected).as_dict()
        for k, v in expected.items():
            print(k)
            print(v)
            print("++++++++++++++++++++++++++++")
            print(actual[k])
            print("----------------------------")

            self.assertEqual(v, actual[k])

    def test_json(self):
        """ round trip test of to_json and from_json """

        self.defect_concentration.to_json_file()
        defect_concentration = \
            DefectConcentration.load_json("defect_concentrations.json")
        os.remove("defect_concentrations.json")
        expected = self.defect_concentration.as_dict()
        actual = defect_concentration.as_dict()
        for k, v in expected.items():
            print(k)
            print(v)
            print("++++++++++++++++++++++++++++")
            print(actual[k])
            print("----------------------------")

            self.assertEqual(v, actual[k])

