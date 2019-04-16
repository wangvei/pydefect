from __future__ import print_function
import copy
from collections import OrderedDict, defaultdict
import os

import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.composition import Composition

from pydefect.analysis.chempotdiag.gas import Gas

# TODO: Once __eq__ is implemented, __hash__ have to be implemented?
# TODO: Now probably we should use pymatgen.composition.


class NotYetCalculatedError(Exception):
    pass


class Compound:
    """
        Object for compound
    """

    def __init__(self, name, elements, composition, energy, gas=None):
        """
        Create a Compound object.

        Args:
            name (str): Name of compound.
            elements (list of string/list of Element(pmg) ):
                Order of elements. It is used to deal composition vector.
            composition (list of float/numpy.array):
                Composition of compound.
                Order of elements is ruled by self.elements.
            energy (float): Energy of compound.
            gas (Gas): If compound is gas, related gas data (like Gas.O2)
                will be used when temperature and pressure are considered.
        """
        self._name = name
        self._elem_comp = OrderedDict()
        number_of_atoms = sum(composition)
        standardized_composition = [c / number_of_atoms for c in composition]
        for elem, comp in zip(elements, standardized_composition):
            self._elem_comp[elem] = comp
        self._energy = energy / number_of_atoms
        self._gas = gas

    @classmethod
    def from_vasp_calculation_files(cls, poscar_path, output_path,
                                    fmt="outcar", is_molecule_gas=True,
                                    elements=None, name=None,
                                    energy_shift=0):
        """
        Args:
            poscar_path (str):
            output_path (str):
            fmt (str): Format of output file.
                       Now "outcar" or "vasprun" modes are implemented.
            is_molecule_gas (bool): Read molecule data as gas.
                                    If path of directory contains "molecule",
                                    gas data are added.
                                    Then zero-point energy and
                                    free energy are considered.
            elements (None/list/str): Considered elements,
                                            like ["Mg", "O"].
                                            If not provided, automatically
                                            determined from poscar.
            name (None/str): Name of the compound.
            energy_shift (float):
        Returns:
            (CompoundsList) CompoundsList object from vasp files.
        """
        pmg_composition = Poscar.from_file(poscar_path).structure.composition
        # HACK: energy type is energy(sigma->0),
        # then we cannot use pymatgen in natural way.
        # TODO: consistent with pydefect?
        gas = None
        if is_molecule_gas and "molecule" in poscar_path:
            for g in Gas:
                if str(g) == pmg_composition.reduced_formula:
                    gas = g
                    break
            if gas is None:
                raise ValueError("{} molecule was not found in database".
                                 format(poscar_path))
        if fmt == "outcar":
            with open(output_path) as fr:
                for line in reversed(fr.readlines()):
                    if line.find("energy(sigma->0)") > 1:
                        energy = float(line.split()[-1]) + energy_shift
                        break
        elif fmt == "vasprun":
            v = Vasprun(output_path)
            energy = v.ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"] +\
                energy_shift
        else:
            raise ValueError("Invalid format {} .".format(fmt))
        if name is None:
            name = pmg_composition.reduced_formula
        if elements is None:
            elements = list(pmg_composition.as_dict.keys())
        composition_vector = np.zeros(len(elements))
        for i, element in enumerate(elements):
            composition_vector[i] = pmg_composition.as_dict[element]
        return cls(name, elements, composition_vector, energy, gas=gas)

    @property
    def name(self):
        """
            (str) Name of compound.
        """
        return self._name

    @property
    def composition(self):
        """
            (numpy.array) Composition of compound.
        """
        return np.array([ec[1] for ec in self._elem_comp.items()])

    @property
    def elements(self):
        """
            (list of str) Ordered list of elements of the composition vector.
        """
        return [ec[0] for ec in self._elem_comp.items()]

    @property
    def dim(self):
        """
            (int) Number of considered atoms.
        """
        return len(self._elem_comp)

    @property
    def energy(self):
        """
            (float) Energy of compound.
        """
        return self._energy

    @property
    def gas(self):
        """
            (Gas) If compound is gas, return related gas data, otherwise None.
        """
        return self._gas

    @gas.setter
    def gas(self, gas):
        """

        Args:
            gas (Gas):

        """
        if not isinstance(gas, Gas):
            raise TypeError("gas must be Gas class, but actually {}".
                            format(type(gas)))
        self._gas = gas

    def gas_energy_shift(self, temperature, pressure):
        """

        Args:
            temperature(float): (K)
            pressure(float): (Pa)

        Returns: (eV/atom)

        """
        if self._gas is None:
            return 0
        else:
            return self._gas.energy_shift(temperature, pressure)

    def free_energy(self, temperature=None, pressure=None):
        """

        Args:
            temperature(float): (K)
            pressure(float): (Pa)

        Returns: (eV/atom)

        """
        if temperature or pressure:
            gas_energy_shift = self.gas_energy_shift(temperature, pressure)
        else:
            gas_energy_shift = 0
        return self.energy + gas_energy_shift

    def standardize_energy(self, standard_energy):
        self._energy = self._energy - standard_energy

    def set_elements(self, elements):
        """
        Change elements of the object.
        This also changes order of composition vector.
        If self.order is subset of input order,
        then new species of atoms are added.
        If self.order is superset of input order and
        compositions of extra elements are nearly zero,
        then those atoms are removed.
        Otherwise, raise ValueError.
        Args:
            elements (list of str): New elements list.
        """
        for e, _ in self._elem_comp.items():
            # case [c_a, c_b, c_c(=nonzero)] -> [c_b, c_a]
            # should be raise error
            if e not in elements and self._elem_comp[e] > 1e-5:
                raise ValueError("Invalid elements list. "
                                 "Composition of removed element {} "
                                 "must be nearly zero ( <=1e-5 ), "
                                 "but actually {} ."
                                 .format(e, self._elem_comp[e]))
        new_elem_comp = OrderedDict()
        for e in elements:
            if e in self.elements:
                new_elem_comp[e] = self._elem_comp[e]
            else:
                new_elem_comp[e] = 0
        self._elem_comp = new_elem_comp

    def __repr__(self):
        return ("Name: {0} , "
                "Elements: {1} , "
                "Composition: {2} , "
                "Energy: {3} , "
                "Gas: {4}"
                .format(self.name, "-".join(self.elements), self.composition,
                        self.energy, self.gas))

    def __eq__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Args:
            other (Compound): Compared compound.
        Returns (bool): If self == other.
        """
        if not isinstance(other, Compound):
            raise TypeError("Compound class can not be"
                            " compared with other class.")
        if self.name != other.name:
            return False
        elif self.energy != other.energy:
            return False
        elif any([e1 != e2 for e1, e2
                  in zip(self.elements, other.elements)]):
            return False
        elif any([c1 != c2 for c1, c2
                  in zip(self.composition, other.composition)]):
            return False
        return True

    def __ne__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Args:
            other (Compound): Compared compound.
        Returns (bool): If self != other.
        """
        return not self == other

    def __lt__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Args:
            other (Compound): Compared compound.
        Returns (bool): If self < other.
        """
        if not isinstance(other, Compound):
            raise TypeError("Compound class can not be"
                            " compared with other class.")
        if self.name != other.name:
            return self.name < other.name
        elif self.energy != other.energy:
            return self.energy < other.energy
        elif any([e1 != e2 for e1, e2 in zip(self.elements, other.elements)]):
            for e1, e2 in zip(self.elements, other.elements):
                if e1 != e2:
                    return e1 < e2
        elif any([c1 != c2 for c1, c2
                  in zip(self.composition, other.composition)]):
            for c1, c2 in zip(self.composition, other.composition):
                if c1 != c2:
                    return c1 < c2
        else:  # all elements are strictly the same.
            return False

    def __gt__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Args:
            other (Compound): Compared compound.
        Returns (bool): If self > other.
        """
        return (not self < other) and self != other

    def __ge__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Args:
            other (Compound): Compared compound.
        Returns (bool): If self >= other.
        """
        return not self < other

    def __le__(self, other):
        """
        This is method for sorting.
        Standard of comparison has not important meaning.
        Args:
            other (Compound): Compared compound.
        Returns (bool): If self > other.
        """
        return not self > other

    def almost_equal(self, other, tol=1e-5):
        """
        Check if two compounds are almost same.
        Args:
            other (Compound): Compared compound.
            tol (float): Tolerance for numeric error of energy and composition.
        Returns (bool): If self almost equals to other.
        """
        if not isinstance(other, Compound):
            raise TypeError("Compound class can not be"
                            " compared with other class.")
        other = copy.deepcopy(other)
        try:
            other.set_elements(self.elements)
        except ValueError:
            return False
        if self.name != other.name:
            return False
        elif abs(self.energy - other.energy) > tol:
            return False
        elif any([e1 != e2 for e1, e2 in zip(self.elements, other.elements)]):
            return False
        elif any([c1 != c2 for c1, c2
                  in zip(self.composition, other.composition)]):
            return False
        return True


class DummyCompoundForDiagram(Compound):
    """
        Object for dummy compound for drawing boundary of diagram.
        It is intended that this class will be used
        for scipy.spatial.halfspaces.
        It will mean halfspace expressed by (x > boundary_energy),
        but halfspace has to be notated by (Ax + b < 0),
        then composition is -1.
    """

    # HACK: This __init__ is needed to allow composition -1,
    # without standardization.
    def __init__(self, name, elements, composition, energy):
        """
        Constructor of dummy boundary.
        Args:
            name (str):
            elements (list of str):
            composition (list of float or numpy.array):
            energy (float):
        """
        self._name = name
        self._elem_comp = OrderedDict()
        for elem, comp in zip(elements, composition):
            self._elem_comp[elem] = comp
        self._energy = energy
        self._gas = None

    @classmethod
    def construct_boundary(cls, elements, element, energy):
        name = str(element) + " boundary"
        composition = [0] * len(elements)
        index = list(elements).index(element)
        composition[index] = -1
        return cls(name, elements, composition, energy)

    @classmethod
    def from_vasp_calculation_files(cls, poscar_path, output_path,
                                    fmt="outcar", elements=None, name=None):
        raise TypeError("Dummy compound can not be constructed from "
                        "DFT calculation.")


class CompoundsList(list):
    _STR_NOT_STANDARDIZED = "Not yet standardized"
    _TYPE_ERROR_MESSAGE = "CompoundsArray must be contains only Compound."

    def __init__(self, *args, pressure=None, temperature=None, **kw):
        """

        Args:
            *args:
            pressure (dict): e.g. {O2: 1} (Pa) When None,
                             default pressure (1e+5 Pa) will be applied.
            temperature (float): (K)
            **kw:
        """
        list.__init__(self, *args, **kw)
        is_different_elements = False
        for c in self:
            if any([e1 != e2
                    for e1, e2 in zip(c.elements, self[0].elements)]):
                is_different_elements = True
        if is_different_elements:
            all_elements = set()
            for c in self:
                all_elements |= set(c.elements)
            for c in self:
                c.set_elements(all_elements)
        else:
            for c in self:
                c.set_elements(self[0].elements)
        for c in self:
            if not isinstance(c, Compound):
                raise TypeError(CompoundsList._TYPE_ERROR_MESSAGE)
        self._element_energies = CompoundsList._STR_NOT_STANDARDIZED
        default_pressure = 1e+5
        if pressure is None:
            self._pressure = defaultdict(lambda: default_pressure, {})
        else:
            self._pressure = defaultdict(lambda: default_pressure, pressure)
        self._temperature = temperature

    def __str__(self):
        return super(CompoundsList, self).__str__() + \
               " at T={}, P={}".format(self._temperature, self._pressure)

    # HACK: If there is not this method, order of elements will change
    #       when deep copied.
    # TODO: Implement the code of (not-deep)copy.
    def __deepcopy__(self, memodict={}):
        return_list = CompoundsList([], temperature=self._temperature,
                                    pressure=self._pressure)
        for c in self:
            return_list.append(c)
        return_list.set_elements(self.elements)
        return_list._element_energies = copy.copy(self._element_energies)
        return return_list

    def append(self, compound):
        # TODO: What should we do when standard_energy is different?
        if not isinstance(compound, Compound):
            raise TypeError(CompoundsList._TYPE_ERROR_MESSAGE)
        if self.elements != compound.elements:
            all_elements = list(set(self.elements) | set(compound.elements))
            for c in self:
                c.set_elements(all_elements)
            copy_compound = copy.deepcopy(compound)
            copy_compound.set_elements(all_elements)
            compound = copy_compound
        list.append(self, compound)

    def extend(self, compounds):
        # TODO: What should we do when standard_energy is different?
        try:
            compounds_list = CompoundsList(compounds)
        except TypeError:
            raise TypeError("Failed to construct CompoundsList from "
                            "given argument {}".format(compounds))
        if self.elements != compounds_list.elements:
            all_elements \
                = list(set(self.elements) | set(compounds_list.elements))
            for c in self:
                c.set_elements(all_elements)
            for c in compounds_list:
                c.set_elements(all_elements)
        list.extend(self, compounds_list)

    def __add__(self, compounds):
        # TODO: What should we do when standard_energy is different?
        try:
            compounds_list = CompoundsList(compounds)
        except TypeError:
            raise TypeError("Failed to construct CompoundsList from "
                            "given argument {}".format(compounds))
        if self.elements != compounds_list.elements:
            new_elements \
                = list(set(self.elements) | set(compounds_list.elements))
            compounds_list.set_elements(new_elements)
        return CompoundsList(list.__add__(self, compounds))

    def __setitem__(self, key, compound):
        if not isinstance(compound, Compound):
            raise TypeError(CompoundsList._TYPE_ERROR_MESSAGE)
        list.__setitem__(self, key, compound)

    @property
    def elements(self):
        if len(self) == 0:
            return []
        return self[0].elements

    @property
    def pressure(self):
        return self._pressure

    @property
    def temperature(self):
        return self._temperature

    def set_elements(self, elements):
        for i, _ in enumerate(self):
            self[i].set_elements(elements)

    @property
    def dim(self):
        return len(self.elements)

    @property
    def element_energies(self):
        if self._element_energies is CompoundsList._STR_NOT_STANDARDIZED:
            raise NotYetCalculatedError("Failed to get element energies "
                                        "because this object has "
                                        "not yet standardized")
        return self._element_energies

    def standardize_energies(self):
        """
        Standardize (let energies of stable simple substance = 0) energy_list.
        """
        comp_list = np.array([c.composition for c in self])
        # energy_list = np.array([c.energy for c in self]) # temporary
        energy_list = np.array(self.free_energies())
        num_composition = len(comp_list[0])

        # TODO: remove are_standard_energies.
        # TODO: If calculations of any simple elements are missed,
        # should be raise error.
        is_simple_substance_found = [False] * num_composition
        element_energies = np.zeros(num_composition)
        # search minimum energies of substances
        for c, e in zip(comp_list, energy_list):
            index = np.where(abs(c - 1) < 1e-8)[0]
            # if index is not empty, c is a simple substance
            if len(index) == 0:  # not simple substance
                continue
            if len(index) >= 2:  # composition is like [1.0, 1.0, 1.0]
                raise ValueError("Composition is {}, it is strange because "
                                 "more than one composition are 1.0. "
                                 "Maybe standardization of composition "
                                 "is failed".format(c))
            if not is_simple_substance_found[index[0]]:
                element_energies[index[0]] = e
                is_simple_substance_found[index[0]] = True
            else:
                if e < element_energies[index[0]]:
                    # More stable simple substance is found.
                    element_energies[index[0]] = e
        not_found_simple = [element for flag, element
                            in zip(is_simple_substance_found, self.elements)
                            if not flag]
        if not_found_simple:
            raise ValueError("Standardization of energies of compounds failed "
                             "because calculation of simple element {} was "
                             "not found.".format(not_found_simple))

        standard_energy_list = np.inner(comp_list, element_energies)
        for i, _ in enumerate(self):
            self[i].standardize_energy(standard_energy_list[i])
        self._element_energies = element_energies

    def gas_energy_shifts(self):
        """
        Free energy energy_shift of compounds if the compound is gas,
        otherwise the shifts are zero.
        Args:
            default_pressure (float):

        Returns (list):
        """
        if self._pressure is None and self._temperature is None:
            return [0] * len(self)
        energy_shifts = []
        for c in self:
            if c.gas:
                energy_shift = \
                    c.gas_energy_shift(self._temperature,
                                       self._pressure[c.gas.formula])
            else:
                energy_shift = 0
            energy_shifts.append(energy_shift)
        return energy_shifts

    def free_energies(self):
        return [c.energy + es
                for c, es in zip(self, self.gas_energy_shifts())]

    def get_indices_and_compounds(self, compound_name):
        """
        Find object of Compound from self by name(str) of compound.
        Args:
            compound_name (str):
        Returns (None/list):
            Matched compound data. If no compounds match, return None.
        """
        compound_name = Composition(compound_name).reduced_formula
        matched_list = [(i, c) for i, c in enumerate(self)
                        if Composition(c.name).reduced_formula == compound_name]
        if len(matched_list) == 0:
            return None
        else:
            return matched_list

    @classmethod
    def from_file(cls, file_name, temperature=None, pressure=None):
        """
        Create a object of CompDat from file.
        Args:
            file_name (str): File name of information of energies.
                The style of file is like below
                     Mg-O
                      Mg  1  0   -2.11327587
                      Mg  1  0   -1.2342344
                       O  0  1   -2.46256238
                    Mg2O  2  1  -14.46256238 ...
                Note that the energies of element substances are essential.
            temperature (float):
            pressure (dict):
        """
        with open(file_name) as f:
            lines = f.readlines()
        elements = lines[0].strip().split('-')
        compounds_list = []
        num_elements = len(elements)
        for l in lines[1:]:
            line = l.split()
            name = line[0]
            num_atoms = np.array([float(i) for i in line[1:num_elements + 1]])
            num_total_atoms = np.sum(num_atoms)
            composition = num_atoms / num_total_atoms
            energy = float(line[num_elements + 1]) / num_total_atoms
            compounds_list.append(Compound(name, elements, composition, energy))
        compounds_list = CompoundsList(compounds_list)
        compounds_list.set_elements(elements)
        return cls(compounds_list, temperature=temperature, pressure=pressure)

    # TODO: make unittest
    def to(self, file_name=None):
        """
        CompoundsList to str. If filename is provided, this method outputs file.
        Args:
            file_name (str): Name of file to output.
                             If it is not provided, string will be returned.
        Returns:
            None/str: If file_name is provided, returns None.
                         Otherwise, returns str.
                The style of file or str must be consistent with from_file.
        """
        return_str = "-".join(self.elements)
        for compound in self:
            return_str += "\n"
            return_str += " {0}".format(compound.name)
            for c in compound.composition:
                return_str += " {0}".format(c)
            return_str += " {0}".format(compound.energy)
        if file_name:
            with open(file_name, 'w') as f:
                f.write(return_str)
            return None
        else:
            return return_str

    @classmethod
    def from_vasp_calculations_files(
            cls, poscar_paths, output_paths, fmt="outcar",
            temperature=None, pressure=None, energy_shift_dict={}):
        """
        Args:
            poscar_paths (list of str):
            output_paths (list of str):
            fmt (str): "outcar" or "vasprun".
            temperature (float): (K)
            pressure (dict): (Pa)
            energy_shift_dict (dict): unit: (eV), e.g. {N2_molecule/OUTCAR: 1}
        Returns:
            (CompoundsList) CompoundsList object from vasp files.
        """
        if len(poscar_paths) != len(output_paths):
            raise ValueError("Number of Paths of poscar ({}) and "
                             "number of paths of output ({}) differs"
                             .format(len(poscar_paths), len(output_paths)))

        # to match energy_shift_dict key,
        # output_paths is transformed to abs_path
        # (for example, it is preferable to match "OUTCAR" and "./OUTCAR")
        energy_shift_dict = {os.path.abspath(k): v
                             for k, v in energy_shift_dict.items()}
        energy_shift_dict = defaultdict(float, energy_shift_dict)
        abs_output_paths = [os.path.abspath(op) for op in output_paths]
        diff = set(energy_shift_dict) - set(abs_output_paths)
        if diff:
            raise ValueError("{} is invalid path.".format(diff))
        # considered_element = set()
        compounds_list = cls([], temperature=temperature, pressure=pressure)
        for poscar_path, output_path in zip(poscar_paths, abs_output_paths):
            try:
                compound = Compound.\
                    from_vasp_calculation_files(poscar_path,
                                                output_path,
                                                fmt=fmt,
                                                energy_shift=energy_shift_dict[
                                                     output_path])
                compounds_list.append(compound)
            except:
                import warnings
                warnings.warn("Could not read {} or {}, skipped".
                              format(poscar_path, output_path))
                continue
        # return cls(compounds_list, temperature=temperature, pressure=pressure)
        return compounds_list

    # TODO: document, unittest
    def almost_equal(self, other, tol=1e-5):
        if not isinstance(other, CompoundsList):
            raise TypeError("Can not be compared because argument is "
                            "{}, not CompoundsList.".format(type(other)))
        if len(self) != len(other):
            return False
        else:
            for c1, c2 in zip(sorted(self), sorted(other)):
                if not c1.almost_equal(c2, tol=tol):
                    return False
        return True

