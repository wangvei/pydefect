from abc import ABCMeta, abstractmethod
from enum import Enum
import json
import math
import os
import re

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar

"""
Data is from NIST Chemistry WebBook, SRD 69
"""


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


class _ShomateParameters:
    """
        Parameters of Shomate Equations.
        Contains range of temperature like (100, 700).
    """

    def __init__(self, temperature_range,
                 a=None, b=None, c=None, d=None,
                 e=None, f=None, g=None, h=None):
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
               "H: {}\n".format(self._temperature_range,
                                self._a, self._b, self._c, self._d,
                                self._e, self._f, self._g, self._h)


class ShomateThermodynamicsFunction(AbstractThermodynamicsFunction):
    """
    Shomate thermodynamics Function.
    """

    def __init__(self, func_param_list):
        """
        Args:
            func_param_list (list):parameters of functions.
            Must contain _Shomate_parameters.
        """
        self._params = func_param_list

    @classmethod
    def from_nist_table(cls, file_path):
        # Is there better expression with regular expression?
        temperature_ranges = []
        param_dicts = []
        with open(file_path, "r") as fr:
            for line in fr:
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
                    for tr, pd in zip(temperature_ranges, param_dicts)])

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
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e = \
            params.a, params.b, params.c, params.d, params.e
        return a + b * t + c * t**2 + d * t**3 + e / (t**2)

    def standard_enthalpy(self, temperature):
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e, f, h = \
            params.a, params.b, params.c, params.d, params.e, \
            params.f, params.h
        return a * t + b * t**2 / 2 + c * t**3 / 3 + \
            d * t**4 / 4 - e / t + f - h

    def standard_entropy(self, temperature):
        params = self.params(temperature)
        t = temperature / 1000
        a, b, c, d, e, g = \
            params.a, params.b, params.c, params.d, \
            params.e, params.g
        return a * math.log(t) + b * t + c * t**2 / 2 + d * t**3 / 3 - \
            e / (2 * t**2) + g


class Gas(Enum):
    H2 = "H2"
    N2 = "N2"
    O2 = "O2"
    F2 = "F2"
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
                "Thermodynamic data of {} molecule"
                " is not found.".format(formula))

        # Read POSCAR
        try:
            path = dirname + "/POSCAR"
            self._structure = Poscar.from_file(path).structure
        except FileNotFoundError:
            raise ValueError(
                "POSCAR of {} molecule"
                " is not found.".format(formula))

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
