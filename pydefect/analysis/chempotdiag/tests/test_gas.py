import unittest

from pydefect.analysis.chempotdiag.gas \
    import FundamentalFrequencies, ShomateThermodynamicsFunction, Gas


class TestFundamentalFrequencies(unittest.TestCase):

    def test_from_yaml(self):
        file_path_nh3 = "../molecules/NH3/fundamental_frequencies.yaml"
        # FundamentalFrequencies.from_yaml(self.file_path_o2)
        nh3 = FundamentalFrequencies.from_yaml(file_path_nh3)
        expected = [{'Frequency': 3506, 'Degeneration': 1},
                    {'Frequency': 1022, 'Degeneration': 1},
                    {'Frequency': 3577, 'Degeneration': 2},
                    {'Frequency': 1691, 'Degeneration': 2}]
        actual = [f.as_dict() for f in nh3]
        self.assertListEqual(expected, actual)


class TestShomateThermodynamicsFunction(unittest.TestCase):

    def setUp(self):
        file_path_o2 = "../molecules/O2/shomate_nist.dat"
        file_path_no2 = "../molecules/NO2/shomate_nist.dat"
        self._o2 = \
            ShomateThermodynamicsFunction.from_nist_table(file_path_o2)
        self._no2 =\
            ShomateThermodynamicsFunction.from_nist_table(file_path_no2)

    def test_heat_capacity(self):
        self.assertAlmostEqual(self._o2.heat_capacity(400), 30.10, 1)
        self.assertAlmostEqual(self._o2.heat_capacity(1500), 36.55, 1)
        self.assertAlmostEqual(self._o2.heat_capacity(5500), 43.43, 1)

    def test_standard_enthalpy(self):
        self.assertAlmostEqual(self._o2.standard_enthalpy(298.15), 0, 1)
        self.assertAlmostEqual(self._o2.standard_enthalpy(400), 3.03/1000, 1)
        self.assertAlmostEqual(self._o2.standard_enthalpy(1500), 40.60/1000, 1)
        self.assertAlmostEqual(self._o2.standard_enthalpy(5500), 202.3/1000, 1)

    def test_standard_entropy(self):
        self.assertAlmostEqual(self._o2.standard_entropy(400), 213.9, 1)
        self.assertAlmostEqual(self._o2.standard_entropy(1500), 258.1, 1)
        self.assertAlmostEqual(self._o2.standard_entropy(5500), 309.6, 1)

    def test_pressure_term(self):
        self.assertAlmostEqual(self._o2.r_ln_p_p0(1e+5), 0)

    def test_zero_point(self):
        self.assertAlmostEqual(self._o2.enthalpy_zero, 8.683)


class TestGas(unittest.TestCase):

    def setUp(self):
        self._o2 = Gas.O2

    def test_read_properties(self):
        print(self._o2._properties)

    def test_str(self):
        self.assertEqual(str(self._o2), "O2")

    def test_thermodynamics(self):
        file_path_o2 = "../molecules/O2/shomate_nist.dat"
        o2_shomate = ShomateThermodynamicsFunction.from_nist_table(file_path_o2)
        self.assertEqual(self._o2.heat_capacity(400),
                         o2_shomate.heat_capacity(400))
        self.assertEqual(self._o2.standard_enthalpy(400),
                         o2_shomate.standard_enthalpy(400))
        self.assertEqual(self._o2.standard_entropy(400),
                         o2_shomate.standard_entropy(400))

    def test_temperature_range(self):
        self.assertEqual(self._o2.min_temperature, 100)
        self.assertEqual(self._o2.max_temperature, 6000)
        self.assertEqual(self._o2.temperature_range, (100, 6000))

    def test_zero_point(self):
        # Expected value is script of Prof. Kumagai
        self.assertAlmostEqual(Gas.O2.zero_point_vibrational_energy, 0.098/2, 2)
        self.assertAlmostEqual(Gas.N2.zero_point_vibrational_energy, 0.146/2, 2)

    def test_energy_shift(self):
        kumagai_script_o2_chemical_potential = \
            {100: -0.052,
             200: -0.243,
             300: -0.450,
             400: -0.667,
             500: -0.892,
             600: -1.124,
             700: -1.361,
             800: -1.604,
             900: -1.850,
             1000: -2.100}
        for t in range(100, 1001, 100):
            this = self._o2.energy_shift(pressure=1e+5, temperature=t)
            pre = kumagai_script_o2_chemical_potential[int(t)] / 2
            print(t, this, pre, (this-pre)/pre)


if __name__ == "__main__":
    unittest.main()

