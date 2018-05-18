import unittest

from pydefect.analysis.chempotdiag.gas \
    import Gas, ShomateThermodynamicsFunction


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
        self.assertAlmostEqual(self._o2.standard_enthalpy(400), 3.03, 1)
        self.assertAlmostEqual(self._o2.standard_enthalpy(1500), 40.60, 1)
        self.assertAlmostEqual(self._o2.standard_enthalpy(5500), 202.3, 1)

    def test_standard_entropy(self):
        self.assertAlmostEqual(self._o2.standard_entropy(400), 213.9, 1)
        self.assertAlmostEqual(self._o2.standard_entropy(1500), 258.1, 1)
        self.assertAlmostEqual(self._o2.standard_entropy(5500), 309.6, 1)


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


if __name__ == "__main__":
    unittest.main()

