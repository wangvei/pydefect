from abc import ABCMeta, abstractmethod
from collections import namedtuple
from enum import Enum
import json
import math
import os
import re

from scipy import constants
from ruamel import yaml

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar

"""
Data is from NIST Chemistry WebBook, SRD 69
"""


class FundamentalFrequency:

    def __init__(self, frequency, degeneration):
        """

        Args:
            frequency (float): frequency (cm^-1)
            degeneration (int): number of degeneration
        """
        self._frequency = frequency
        self._degeneration = degeneration

    def __str__(self):
        return "Frequency: {0} /cm, Degeneration: {1}".\
            format(self._frequency, self._degeneration)

    @property
    def frequency(self):
        """

        Returns (float): cm^(-1)

        """
        return self._frequency

    @property
    def degeneration(self):
        """

        Returns (int): number of degeneration

        """
        return self._degeneration

    def as_dict(self):
        return {"Frequency": self._frequency,
                "Degeneration": self._degeneration}

    @property
    def zero_point_vibrational_energy(self):
        """

        Returns (float): (eV)

        """
        h = constants.physical_constants["Planck constant in eV s"][0]
        return 1/2 * h * self._degeneration * self._frequency * 1e+2 *\
            constants.speed_of_light


class FundamentalFrequencies:

    def __init__(self, frequencies, reference):
        """
        Args:
            frequencies (list of FundamentalFrequency):
            reference (str):
        """
        self._frequencies = frequencies
        self._reference = reference

    def __iter__(self):
        return self._frequencies.__iter__()

    @classmethod
    def from_yaml(cls, file_path):
        """

        Args:
            file_path (str):

        Returns (FundamentalFrequencies):

        """

        with open(file_path, 'r') as fr:
            read_data = yaml.safe_load(fr)

        frequencies = \
            [FundamentalFrequency(f["Frequency"], f["Degeneration"])
             for f in read_data["Frequencies"]]
        try:
            reference = read_data["Reference"]
        except KeyError:
            reference = None

        return cls(frequencies, reference)

    @property
    def zero_point_vibrational_energy(self):
        """

        Returns (float): E_ZPVE (eV)

        """
        return sum(f.zero_point_vibrational_energy for f in self._frequencies)


class AbstractThermodynamicsFunction(metaclass=ABCMeta):
    # TODO: implement abstract factory

    @abstractmethod
    def heat_capacity(self, temperature):
        pass

    @abstractmethod
    def standard_enthalpy(self, temperature):
        pass

    @abstractmethod
    def standard_entropy(self, temperature):
        pass

    @abstractmethod
    def r_ln_p_p0(self, pressure):
        pass

    @property
    @abstractmethod
    def enthalpy_zero(self):
        """

        Returns: Energy at pressure = self._reference_pressure, T=0K

        """


class _ShomateParameters:
    """
        Parameters of Shomate Equations.
        Contains range of temperature like (100, 700).
    """

    def __init__(self, temperature_range,
                 a, b, c, d, e, f, g, h):
        """

        Args:
            temperature_range (tuple): Unit is (K).
            a (float): parameter A.
            b (float): parameter B.
            c (float): parameter C.
            d (float): parameter D.
            e (float): parameter E.
            f (float): parameter F.
            g (float): parameter G.
            h (float): parameter H.
        """
        self._temperature_range = temperature_range
        self._a = a
        self._b = b
        self._c = c
        self._d = d
        self._e = e
        self._f = f
        self._g = g
        self._h = h

    @property
    def temperature_range(self):
        return self._temperature_range

    @property
    def min_temperature(self):
        return self._temperature_range[0]

    @property
    def max_temperature(self):
        return self._temperature_range[1]

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def d(self):
        return self._d

    @property
    def e(self):
        return self._e

    @property
    def f(self):
        return self._f

    @property
    def g(self):
        return self._g

    @property
    def h(self):
        return self._h

    def can_apply(self, temp):
        if self.min_temperature <= temp <= self.max_temperature:
            return True
        else:
            return False

    def __str__(self):
        return "Temperature_range: {}\n"\
               "A: {}\n"\
               "B: {}\n"\
               "C: {}\n"\
               "D: {}\n"\
               "E: {}\n"\
               "F: {}\n"\
               "G: {}\n"\
               "H: {}\n".\
            format(self._temperature_range,
                   self._a, self._b, self._c, self._d,
                   self._e, self._f, self._g, self._h)


class ShomateThermodynamicsFunction(AbstractThermodynamicsFunction):
    """
    Shomate thermodynamics Function.
    """

    def __init__(self, func_param_list, enthalpy_zero,
                 reference_pressure=1e+5):
        """
        Args:
            func_param_list (list): parameters of functions.
                                    Must contain _Shomate_parameters.
            enthalpy_zero (float): enthalpy (kJ/mol) at T=0K,
                                   p=reference_pressure
            reference_pressure (float): reference_pressure (Pa)

        """
        self._params = func_param_list
        self._enthalpy_zero = enthalpy_zero
        self._reference_pressure = reference_pressure

    @classmethod
    def from_nist_table(cls, file_path):
        """
        shomate_nist.dat is written e.g.
        Enthalpy_zero   -8.825
        Temperature (K)	298. - 6000.
        A	31.44510
        B	8.413831
        C	-2.778850
        D	0.218104
        E	-0.211175
        F	-10.43260
        G	237.2770
        H	0.000000
        Reference	Chase, 1998
        Comment	Data last reviewed in June, 1982

        where Enthalpy_zero is kJ/mol

        Args:
            file_path (str):

        Returns:

        """
        # TODO: unit of file is confusing
        # Is there better expression with regular expression?
        temperature_ranges = []
        param_dicts = []
        enthalpy_zero = None
        with open(file_path, "r") as fr:
            for line in fr:
                if enthalpy_zero is None:
                    is_tmp_line = \
                        bool(re.match(r"Enthalpy_zero", line))
                    if is_tmp_line:
                        matched = re.findall(r"[\d\.]+", line)
                        enthalpy_zero = float(matched[0])
                    continue

                if not temperature_ranges:
                    is_tmp_line = \
                        bool(re.match(r"Temperature\s+\(K\)\s+", line))
                    if is_tmp_line:
                        matched = re.findall(r"[\d\.]+", line)
                        for i in range(int(len(matched)/2)):
                            temperature_ranges.append(
                                (float(matched[2*i]), float(matched[2*i+1])))
                        param_dicts = \
                            [dict() for _ in range(len(temperature_ranges))]
                else:
                    symbol_match = re.search(r"\s*([A-H])\s+", line)
                    if symbol_match:
                        symbol = symbol_match.groups()[0].lower()
                        values = re.findall(r"[\-\d\.]+", line)
                        if len(values) != len(temperature_ranges):
                            raise ValueError("Number of temperature_ranges"
                                             "({})"
                                             "and parameters ({})"
                                             "is not consistent.".
                                             format(len(temperature_ranges),
                                                    line))
                        for i in range(len(values)):
                            param_dicts[i][symbol] = float(values[i])
        return cls([_ShomateParameters(tr, **pd)
                    for tr, pd in zip(temperature_ranges, param_dicts)],
                   enthalpy_zero)

    def params(self, temperature, exception_out_of_apply_range=True):
        """
        Return parameters of function at temperature like {a:1.1, b:2.2, ...}
        Args:
            temperature:
            exception_out_of_apply_range(bool):
        Returns:

        """
        params = None
        for p in self._params:
            if p.can_apply(temperature):
                params = p
        if params is None:
            if exception_out_of_apply_range:
                raise ValueError("Temperature {} is out of temperature "
                                 "range to apply".format(temperature))
            else:
                max_index, max_data = \
                    max(enumerate(self._params),
                        key=lambda x: x[1].max_temperature)
                min_index, min_data = \
                    min(enumerate(self._params),
                        key=lambda x: x[1].min_temperature)
                if temperature < max_data.max_temperature:
                    params = self._params[max_index]
                elif temperature > min_data.min_temperature:
                    params = self._params[min_index]
                else:
                    raise ValueError("Temperature {} is strange. "
                                     "Parameters to apply were not found.".
                                     format(temperature))
        return params

    @property
    def temperature_range(self):
        return self.min_temperature, self.max_temperature

    @property
    def max_temperature(self):
        max_data = max(self._params, key=lambda x: x.max_temperature)
        return max_data.max_temperature

    @property
    def min_temperature(self):
        min_data = min(self._params, key=lambda x: x.min_temperature)
        return min_data.min_temperature

    def heat_capacity(self, temperature):
        """

        Args:
            temperature (float): (K)

        Returns (float): J/(K mol)

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e = \
            params.a, params.b, params.c, params.d, params.e
        return a + b * t + c * t**2 + d * t**3 + e / (t**2)

    def standard_enthalpy(self, temperature):
        """

        Args:
            temperature (float): (K)

        Returns (float): J/mol

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e, f, h = \
            params.a, params.b, params.c, params.d, params.e, \
            params.f, params.h
        enthalpy_kj = a * t + b * t**2 / 2 + c * t**3 / 3 + \
            d * t**4 / 4 - e / t + f - h
        return enthalpy_kj / 1000

    def standard_entropy(self, temperature):
        """

        Args:
            temperature (float): (K)

        Returns (float): J/(mol K)

        """
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e, g = \
            params.a, params.b, params.c, params.d, \
            params.e, params.g
        return a * math.log(t) + b * t + c * t**2 / 2 + d * t**3 / 3 - \
            e / (2 * t**2) + g

    def r_ln_p_p0(self, pressure):
        """
        R ln(p/p0)
        Args:
            pressure (float): MPa

        Returns (float): (J/(mol K))

        """
        return constants.R * math.log(pressure / self._reference_pressure)

    @property
    def enthalpy_zero(self):
        """

        Returns: Energy at pressure = self._reference_pressure, T=0K (J/mol)

        """
        return self._enthalpy_zero / 1000


class InvalidFileError(Exception):
    pass


class Gas(Enum):
    H2 = "H2"
    N2 = "N2"
    O2 = "O2"
    F2 = "F2"
    P2 = "P2"
    P4 = "P4"
    H2O = "H2O"
    NH3 = "NH3"
    NO2 = "NO2"

    def __init__(self, formula):
        self._formula = "%s" % formula

        dirname = os.path.dirname(__file__) + "/molecules/" + formula

        # Read prior_info.json
        try:
            path = dirname + "/prior_info.json"
            with open(path, "r") as fr:
                self._properties = json.load(fr)
        except FileNotFoundError:
            raise ValueError(
                "Data of {} molecule is not found.".format(formula))

        # Read data of thermodynamics function
        try:
            path = dirname + "/shomate_nist.dat"
            # TODO: Implement flexibly with abstract class
            self._thermodynamics_function = \
                ShomateThermodynamicsFunction.from_nist_table(path)
        except FileNotFoundError:
            raise ValueError(
                "shomate_nist.dat of {} molecule"
                " is not found.".format(formula))
        except:
            raise InvalidFileError(
                "shomate_nist.dat of {} molecule".
                format(formula))

        # Read POSCAR
        try:
            path = dirname + "/POSCAR"
            self._structure = Poscar.from_file(path).structure
        except FileNotFoundError:
            raise ValueError(
                "POSCAR of {} molecule"
                " is not found.".format(formula))
        except:
            raise InvalidFileError(
                "Failed to read POSCAR of {} molecule".
                format(formula))

        # Read Zero point vibrational frequency
        try:
            path = dirname + "/fundamental_frequencies.yaml"
            self._fundamental_frequencies = \
                FundamentalFrequencies.from_yaml(path)
        except FileNotFoundError:
            raise ValueError(
                "fundamental_frequencies.yaml of {} molecule"
                " is not found.".format(formula))
        except:
            raise InvalidFileError(
                "Failed to read fundamental_frequencies.yaml of {} molecule".
                format(formula))

    def __str__(self):
        return self.value

    @property
    def composition(self):
        return Composition(str(self))

    @property
    def properties(self):
        return self._properties

    @property
    def structure(self):
        return self._structure

    @property
    def temperature_range(self):
        return self._thermodynamics_function.temperature_range

    @property
    def min_temperature(self):
        return self._thermodynamics_function.min_temperature

    @property
    def max_temperature(self):
        return self._thermodynamics_function.max_temperature

    def r_ln_p0_p(self, pressure):
        return self._thermodynamics_function.r_ln_p_p0(pressure)

    @property
    def formula(self):
        return self._formula

    @property
    def n_atom(self):
        return Composition(self._formula).num_atoms

    @property
    def fundamental_frequency(self):
        """

        Returns: Frequency of vibration at zero point (cm^-1)

        """
        return self._fundamental_frequencies

    @property
    def zero_point_vibrational_energy(self):
        """

        Returns: Energy of vibrational energy at zero point (eV/atom)

        """
        return self._fundamental_frequencies.zero_point_vibrational_energy /\
            self.n_atom

    def heat_capacity(self, temperature):
        """
        Args:
            temperature(float): Unit is (K)

        Returns(float): Heat capacity C_p (J/mol*K).

        """
        return self._thermodynamics_function.heat_capacity(temperature)

    def standard_enthalpy(self, temperature):
        """
        Standard enthalpy from room temperature.
        Args:
            temperature(float): Unit is (K)

        Returns(float): Standard enthalpy (kJ/mol)
                        from room temperature (T=298.15(K))

        """
        return self._thermodynamics_function.standard_enthalpy(temperature)

    def standard_entropy(self, temperature):
        """
        Args:
            temperature(float): Unit is (K)

        Returns(float): Standard entropy (J/mol*K).

        """
        return self._thermodynamics_function.standard_entropy(temperature)

    def free_energy_shift(self, temperature):
        """
        G(T, P0) - H(T=0K, P0) (eV/mol)

        Returns:

        """
        val_joule =\
            self.standard_enthalpy(temperature) - \
            temperature * self.standard_entropy(temperature) - \
            self._thermodynamics_function.enthalpy_zero
        joule_to_ev = \
            constants.physical_constants["joule-electron volt relationship"][0]
        return val_joule * joule_to_ev / (self.n_atom * constants.Avogadro)

    def energy_shift(self, temperature, pressure):
        """

        Args:
            temperature (float): (K)
            pressure (float): (Pa)

        Returns (float): Free energy shift (eV/atom)

        """
        joule_to_ev =\
            constants.physical_constants["joule-electron volt relationship"][0]

        r_ln_p0_p_ev_per_atom = \
            joule_to_ev * self.r_ln_p0_p(pressure) / \
            (self.n_atom * constants.Avogadro)

        return self.zero_point_vibrational_energy + \
            self.free_energy_shift(temperature) + \
            temperature * r_ln_p0_p_ev_per_atom
