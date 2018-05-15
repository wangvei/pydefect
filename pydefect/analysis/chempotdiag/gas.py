from abc import ABCMeta, abstractmethod
from enum import Enum
import json
import math
import os
import re
import sys

from monty.re import regrep
from pymatgen.core.composition import Composition

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

    def __init__(self, temperature_range,
                 a=None, b=None, c=None, d=None,
                 e=None, f=None, g=None, h=None):
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
        if self._temperature_range[0] <= temp <= self._temperature_range[1]:
            return True
        else:
            return False

    def __str__(self):
        return f"Temperature_range: {self._temperature_range}\n"\
               f"A: {self._a}\n"\
               f"B: {self._b}\n"\
               f"C: {self._c}\n"\
               f"D: {self._d}\n"\
               f"E: {self._e}\n"\
               f"F: {self._f}\n"\
               f"G: {self._g}\n"\
               f"H: {self._h}\n"


class ShomateThermodynamicsFunction(AbstractThermodynamicsFunction):
    """
    Shomate thermodynamics Function.
    """

    def __init__(self, func_param_list):
        """
        Args:
            **func_param_list (list):parameters of functions.
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
                            raise ValueError(f"Number of temperature_ranges"
                                             f"({len(temperature_ranges)})"
                                             f"and parameters ({line})"
                                             f"is not consistent.")
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
                max_index, max_t = \
                    max(enumerate(self._params),
                        key=lambda x: x[1]["temperature_range"][1])
                min_index, min_t = \
                    min(enumerate(self._params),
                        key=lambda x: x[1]["temperature_range"][0])
                if temperature < max_t:
                    params = self._params[max_index]
                elif temperature > min_t:
                    params = self._params[min_index]
                else:
                    raise ValueError("Temperature {} is strange. "
                                     "Parameters to apply were not found.".
                                     format(temperature))
        return params

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
    NH3 = "NH3"
    NO2 = "NO2"

    def __init__(self, formula):
        self._formula = "%s" % formula

        # Read property.json
        try:
            path = os.path.dirname(__file__) + "/molecules/" + \
                   formula + "/property.json"
            with open(path, "r") as fr:
                self._properties = json.load(fr)
        except FileNotFoundError:
            raise ValueError(
                "Data of {} molecule is not found.".format(formula))

        # Read data of thermodynamics function
        try:
            path = os.path.dirname(__file__) + "/molecules/" + \
                   formula + "/shomate_nist.dat"
            # TODO: Implement flexibly with abstract class
            self._thermodynamics_function = \
                ShomateThermodynamicsFunction.from_nist_table(path)
        except FileNotFoundError:
            raise ValueError(
                "Thermodynamic data of {} molecule"
                " is not found.".format(formula))

    def __str__(self):
        return self.value

    @property
    def composition(self):
        return Composition(str(self))

    @property
    def data(self):
        return self._data

    def heat_capacity(self, temperature):
        """
        Data is from NIST Chemistry WebBook, SRD 69
        Args:
            temperature:

        Returns:

        """
        return self._thermodynamics_function.heat_capacity(temperature)

    def standard_enthalpy(self, temperature):
        """
        From room temperature.
        Data is from NIST Chemistry WebBook, SRD 69
        Args:
            temperature:

        Returns:

        """
        return self._thermodynamics_function.standard_enthalpy(temperature)

    def standard_entropy(self, temperature):
        """
        Data is from NIST Chemistry WebBook, SRD 69
        Args:
            temperature:

        Returns:

        """
        return self._thermodynamics_function.standard_entropy(temperature)




