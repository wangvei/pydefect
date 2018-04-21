#!/usr/bin/env python

from enum import Enum, unique
from functools import reduce
from itertools import product, groupby
import json
import math
from operator import itemgetter

import matplotlib
import matplotlib.pyplot as plt
from monty.json import MontyEncoder
from monty.serialization import loadfn
import numpy as np
from pymatgen.core.lattice import Lattice
import scipy
import scipy.constants as sconst
from scipy.stats import mstats

from pydefect.core.supercell_dft_results import distance_list, defect_center


"""
This module provides functions used to correct error of defect formation energy
due to finite supercell-size effect.
"""

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "December 8, 2017"


def calc_distance_two_planes(lattice_vectors):
    # (a_i \times a_j) \ddot a_k / |a_i \times  a_j|
    distance = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_a_j = np.cross(lattice_vectors[i - 2], lattice_vectors[i - 1])
        a_k = lattice_vectors[i]
        distance[i] = abs(np.dot(a_i_a_j, a_k)) / np.linalg.norm(a_i_a_j)
    return distance


def calc_max_sphere_radius(lattice_vectors):
    # Maximum radius of a sphere fitting inside the unit cell.
    return max(calc_distance_two_planes(lattice_vectors)) / 2.0


class Ewald:

    def __init__(self, lattice_matrix, dielectric_tensor, ewald_param,
                 prod_cutoff_fwhm=25.0,
                 num_real_lattice=None, num_reciprocal_lattice=None):
        self._lattice_matrix = lattice_matrix
        self._reciprocal_lattice_matrix =\
            Lattice(self._lattice_matrix).reciprocal_lattice.matrix
        self._dielectric_tensor = dielectric_tensor
        self._ewald_param = ewald_param
        self._prod_cutoff_fwhm = prod_cutoff_fwhm
        if num_real_lattice is None:
            num_real_lattice = \
                self._num_real_lattice = sum(1 for _ in
                                             self.generate_neighbor_lattices())
        self._num_real_lattice = num_real_lattice
        if num_reciprocal_lattice is None:
            num_reciprocal_lattice = \
                sum(1 for _ in self.generate_neighbor_lattices(
                    include_self=False, is_reciprocal=False))
        self._num_reciprocal_lattice = num_reciprocal_lattice

    @property
    def lattice_matrix(self):
        return self._lattice_matrix

    @property
    def reciprocal_matrix(self):
        return self._reciprocal_lattice_matrix

    @property
    def dielectric_tensor(self):
        return self._dielectric_tensor

    @property
    def ewald_param(self):
        return self._ewald_param

    @property
    def max_r_vector_norm(self):
        root_det_dielectric = math.sqrt(np.linalg.det(self._dielectric_tensor))
        vol = Lattice(self._lattice_matrix).volume
        cube_root_vol = math.pow(vol, 1/3)
        ewald = self._ewald_param / cube_root_vol * root_det_dielectric
        return self._prod_cutoff_fwhm / ewald

    @property
    def max_g_vector_norm(self):
        root_det_dielectric = math.sqrt(np.linalg.det(self._dielectric_tensor))
        vol = Lattice(self._lattice_matrix).volume
        cube_root_vol = math.pow(vol, 1/3)
        ewald = self._ewald_param / cube_root_vol * root_det_dielectric
        return 2 * ewald * self._prod_cutoff_fwhm

    @property
    def num_real_lattice(self):
        return self._num_real_lattice

    @property
    def num_reciprocal_lattice(self):
        return self._num_reciprocal_lattice

    def as_dict(self):
        d = {"lattice_matrix": self._lattice_matrix,
             "dielectric_tensor": self._dielectric_tensor,
             "ewald_param": self._ewald_param,
             "prod_cutoff_fwhm": self._prod_cutoff_fwhm,
             "num_real_lattice": self._num_real_lattice,
             "num_reciprocal_lattice": self._num_reciprocal_lattice}
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["lattice_matrix"], d["dielectric_tensor"],
                   d["ewald_param"], d["prod_cutoff_fwhm"],
                   d["num_real_lattice"], d["num_reciprocal_lattice"])

    def to_json_file(self, filename):
        """
        Dump a json file.

        Args:
            filename (str):

        Returns:

        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @classmethod
    def from_optimization(cls,
                          structure,
                          dielectric_tensor,
                          initial_value=None,
                          convergence=1.05,
                          prod_cutoff_fwhm=25.0):
        """
        Get optimized ewald parameter.
        Args:
            structure (pmg structure): structure
            dielectric_tensor (numpy 3x3 array): dielectric tensor
            initial_value (float): Initial guess of parameter.
            convergence (float):
                If 1/convergence < n_(real)/n_(reciprocal) < convergence,
                where n_(real) and n_(reciprocal) is number of real lattices
                and reciprocal lattices, finishes optimization and
                returns ewald_param.
            prod_cutoff_fwhm (float):
                product of cutoff radius of G-vector and gaussian FWHM.
                Increasing this value, calculation will be more accurate, but
                slower.
        Returns (float):
            Optimized ewald_param.
        """
        root_det_dielectric = math.sqrt(np.linalg.det(dielectric_tensor))
        real_lattice = structure.lattice.matrix
        reciprocal_lattice = \
            structure.lattice.reciprocal_lattice.matrix
        cube_root_vol = \
            math.pow(structure.lattice.volume, 1/3)
        if initial_value is not None:
            ewald_param = initial_value
        else:
            # determine initial ewald parameter to satisfy following:
            # max_int(Real) = max_int(Reciprocal)
            # in generate_neighbor_lattices function.
            # Left term:
            # max_int(Real) = 2 * x * Y  / l_r where x, Y, and l_r are ewald,
            # prod_cutoff_fwhm, and axis length of real lattice, respectively.
            # Right term:
            # max_int(reciprocal) = Y  / (x * l_g)
            # where l_g is axis length of reciprocal lattice, respectively.
            # Then, x = sqrt(l_g / l_r / 2)
            # gmean : geometric mean, like (a1 * a2 * a3)^(1/3)
            l_r = mstats.gmean([np.linalg.norm(v) for v in real_lattice])
            l_g = mstats.gmean([np.linalg.norm(v) for v in reciprocal_lattice])
            ewald_param \
                = np.sqrt(l_g / l_r / 2) * \
                cube_root_vol / root_det_dielectric
        while True:
            ewald = ewald_param / cube_root_vol * root_det_dielectric
            # count number of real lattice
            max_r_vector_norm = prod_cutoff_fwhm / ewald
            num_real_lattice = \
                sum(1 for _ in
                    cls._generate_neighbor_lattices(real_lattice,
                                                    max_r_vector_norm,
                                                    include_self=True))
            # count number of reciprocal lattice
            max_g_vector_norm = 2 * ewald * prod_cutoff_fwhm
            num_reciprocal_lattice = \
                sum(1 for _ in
                    cls._generate_neighbor_lattices(reciprocal_lattice,
                                                    max_g_vector_norm,
                                                    include_self=False))
            diff_real_reciprocal = num_real_lattice / num_reciprocal_lattice
            if 1 / convergence < diff_real_reciprocal < convergence:
                return cls(real_lattice, dielectric_tensor, ewald_param,
                           prod_cutoff_fwhm,
                           num_real_lattice, num_reciprocal_lattice)
            else:
                ewald_param *= diff_real_reciprocal ** 0.17

    def generate_neighbor_lattices(self, include_self=True, shift=np.zeros(3),
                                   is_reciprocal=False):

        if not is_reciprocal:
            lattice = self.lattice_matrix
            max_norm = self.max_r_vector_norm
        else:
            lattice = self.reciprocal_matrix
            max_norm = self.max_g_vector_norm
        for v in Ewald._generate_neighbor_lattices(lattice, max_norm,
                                                   include_self, shift):
            yield v

    # for searching ewald parameter
    @staticmethod
    def _generate_neighbor_lattices(lattice_vectors,
                                    max_length,
                                    include_self,
                                    shift=np.zeros(3)):
        """
        Generator of a set of lattice vectors within the max length.
        Note that angles between any two axes are assumed to be between 60 and
        120 deg.
        Args:
            lattice_vectors (np.ndarray): 3x3 matrix.
            max_length (float): Max length to search lattice set.
            include_self (bool): Flag whether (0, 0, 0) will be yield.
            shift (np.ndarray): Lattice_vector + shift will be yielded.
                                Should be specify by Cartesian vector.
                                Defaults to zero vector.
        Yields (np.ndarray): Cartesian vector of lattice point.
        """
        #print(lattice_vectors)
        max_int = [int(max_length / np.linalg.norm(lattice_vectors[i])) + 1
                   for i in range(3)]
        for index in product(range(-max_int[0], max_int[0] + 1),
                             range(-max_int[1], max_int[1] + 1),
                             range(-max_int[2], max_int[2] + 1)):
            if (not include_self) and index == (0, 0, 0):
                continue
            cart_vector = np.dot(lattice_vectors.transpose(), np.array(index))
            if np.linalg.norm(cart_vector) < max_length:
                yield cart_vector + shift


@unique
class CorrectionMethod(Enum):
    extended_fnv = "Extended FNV"

    def str(self):
        return str(self.value)


class Correction:

    def __init__(self, method, ewald, lattice_energy, diff_ave_pot, alignment,
                 symbols_without_defect, distances_from_defect,
                 difference_electrostatic_pot, model_pot,
                 manually_set_energy=None):
        """
        Args:
            method (CorrectionMethod):
            ewald (Ewald):
            lattice_energy (float):
            diff_ave_pot (float):
            alignment (float):
            symbols_without_defect (list of str):
            distances_from_defect (list of float):
            model_pot (list of float):
            difference_electrostatic_pot (list of float):
            manually_set_energy (float or None):
        """

        # error check just in case (should be removed in the future)
        if not len(symbols_without_defect) == len(distances_from_defect) \
                == len(difference_electrostatic_pot) == len(model_pot):
            raise IndexError(
                "Lengths of symbols({0}), distances({1}), "
                "electro_static_pot({2}), model_pot differs({3}) differ.".
                format(len(symbols_without_defect), len(distances_from_defect),
                       len(difference_electrostatic_pot), len(model_pot)))
        self._method = method
        self._ewald = ewald
        self._lattice_energy = lattice_energy
        self._diff_ave_pot = diff_ave_pot
        self._alignment = alignment
        self._symbols_without_defect = symbols_without_defect
        self._distances_from_defect = list(distances_from_defect)
        self._difference_electrostatic_pot = list(difference_electrostatic_pot)
        self._model_pot = list(model_pot)
        self._manually_set_energy = manually_set_energy

    @property
    def method(self):
        return self._method

    @property
    def ewald(self):
        return self._ewald

    @property
    def manually_set_energy(self):
        return self._manually_set_energy

    @manually_set_energy.setter
    def manually_set_energy(self, value):
        self._manually_set_energy = value

    @property
    def lattice_energy(self):
        return self._lattice_energy

    @property
    def diff_ave_pot(self):
        return self._diff_ave_pot

    @property
    def alignment(self):
        return self._alignment

    @property
    def symbols_without_defect(self):
        return list(self._symbols_without_defect)

    @property
    def distances_from_defect(self):
        return list(self._distances_from_defect)

    @property
    def difference_electrostatic_pot(self):
        return list(self._difference_electrostatic_pot)

    @property
    def model_pot(self):
        return list(self._model_pot)

    @property
    def max_sphere_radius(self):
        lattice_vectors = self._ewald.lattice_matrix
        return calc_max_sphere_radius(lattice_vectors)

    # TODO: if remote connect and can't be display figure by X,
    # matplotlib maybe raise exception.
    def plot_distance_vs_potential(self, file_name=None):
        property_without_defect = list(zip(self._symbols_without_defect,
                                           self._distances_from_defect,
                                           self._difference_electrostatic_pot))
        property_without_defect = sorted(property_without_defect,
                                         key=itemgetter(0))
        # points_dictionary is like {"Mg": [-0.22, -0.7,...], # "O":[...]}
        points_dictionary = {}
        for k, g in groupby(property_without_defect, key=itemgetter(0)):
            values = [(x, y) for _, x, y in g]
            points_dictionary[k] = values
        # print(points_dictionary)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # elemental electrostatic potential
        for i, (symbol, points) in enumerate(points_dictionary.items()):
            # print("symbol = ")
            # print(symbol)
            # print("points = ")
            # print(points)
            x_set = np.array([point[0] for point in points])
            y_set = np.array([point[1] for point in points])
            # print(symbol, x_set, y_set)
            # set color
            gradation = i / len(points_dictionary.items())
            color = tuple(np.array([0, gradation, 1-gradation]))
            ax.scatter(x_set, y_set, c=color, marker="x", label=symbol)
        # model potential
        ax.scatter(self._distances_from_defect, self._model_pot,
                   c=(1, 0, 0), marker=".", label="model potential")
        # difference between model potential and electrostatic potential
        diff_model_electrostatic = \
            np.array(self._difference_electrostatic_pot) -\
            np.array(self._model_pot)
        ax.scatter(self._distances_from_defect, diff_model_electrostatic,
                   c=(1, 1, 1), marker="o", label="model - electrostatic",
                   linewidth="0.6", edgecolors="black"
                   )
        # potential difference
        point_x = [self.max_sphere_radius, self._distances_from_defect[-1]]
        point_y = [self.diff_ave_pot, self.diff_ave_pot]
        ax.plot(point_x, point_y,
                c=(0, 0, 0), label="potential difference")
        ax.legend(loc="upper left")
        plt.title("Distance vs potential")
        if file_name:
            plt.savefig(file_name, format="eps")
            plt.close(fig)
        else:
            plt.show()

    def as_dict(self):
        """
        Returns (dict):
        """
        d = {"method": str(self._method),
             "ewald": self._ewald,
             "lattice_energy": self._lattice_energy,
             "diff_ave_pot": self._diff_ave_pot,
             "alignment": self._alignment,
             "symbols_without_defect": self._symbols_without_defect,
             "distances_from_defect": list(self._distances_from_defect),
             "difference_electrostatic_pot": list(self._difference_electrostatic_pot),
             "model_pot": list(self._model_pot),
             "manually_set_energy": self._manually_set_energy}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectEntry class object from a dictionary.
        """
        method = CorrectionMethod(d["method"])

        return cls(method, d["ewald"], d["diff_ave_pot"], d["alignment"],
                   d["manually_set_energy"])

    def to_json_file(self, filename):
        """
        Dump a json file.

        Args:
            filename (str):

        Returns:

        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @classmethod
    def compute_correction(cls, defect_entry, defect_dft, perfect_dft,
                           unitcell_dft, ewald=None,
                           method=CorrectionMethod.extended_fnv):
        """
        Args
            defect_entry(DefectEntry):
            defect_dft(SupercellDftResults):
            perfect_dft(SupercellDftResults):
            unitcell_dft(UnitcellDftResults):
            ewald(Ewald):
            method(CorrectionMethod):
        """
        # search ewald
        if ewald is None:
            ewald = Ewald.from_optimization(
                perfect_dft.final_structure,
                unitcell_dft.total_dielectric_tensor)
        if method == CorrectionMethod.extended_fnv:
            return cls.compute_alignment_by_extended_fnv(defect_entry,
                                                         defect_dft,
                                                         perfect_dft,
                                                         unitcell_dft,
                                                         ewald)
        else:
            raise ValueError("Method named {0} is not implemented.".
                             format(method))

    @classmethod
    def compute_alignment_by_extended_fnv(cls,
                                          defect_entry,
                                          defect_dft,
                                          perfect_dft,
                                          unitcell_dft,
                                          ewald):

        """
        Corrects error of formation energy of point defect due to
        finite-supercell effect.
            defect_entry(DefectEntry):
            defect_dft(SupercellDftResults):
            perfect_dft(SupercellDftResults):
            unitcell_dft(UnitcellDftResults):
            ewald_param(Ewald):
        Returns (float): Corrected energy by extended FNV method.
        """

        dielectric_tensor = unitcell_dft.total_dielectric_tensor
        perfect_structure = perfect_dft.final_structure
        diff_ep = [-ep for ep in
                   defect_dft.relative_potential(perfect_dft, defect_entry)
                   if ep is not None]
        # print("defect_dft.final_structure")
        # print(defect_dft.final_structure)
        atomic_position_without_defect =\
            [defect_dft.final_structure.frac_coords[i]
             for i, j in enumerate(defect_entry.atom_mapping_to_perfect)
             if j is not None]
        symbols_without_defect =\
            [defect_dft.final_structure.sites[i].specie.symbol
             for i, j in enumerate(defect_entry.atom_mapping_to_perfect)
             if j is not None]

        charge = defect_entry.charge
        axis = np.array(perfect_structure.lattice.matrix)
        volume = perfect_structure.lattice.volume
        # print("defect_entry.name = ")
        # print(defect_entry.name)
        defect_coords = defect_center(defect_entry, defect_dft.final_structure)
        # print("defect_coords")
        # print(defect_coords)
        # print("distance_list")
        # print(distance_list(defect_dft.final_structure, defect_coords))
        distances_from_defect = \
            [distance_list(defect_dft.final_structure, defect_coords)[i]
             for i, j in enumerate(defect_entry.atom_mapping_to_perfect)
             if j is not None]
        # print("distances_from_defect")
        # print(distances_from_defect)

        # TODO: check ewald or ewald_param?
        # model potential and lattice energy
        coeff = charge * sconst.elementary_charge * \
            1.0e10 / sconst.epsilon_0  # [V]
        model_pot = [None for _ in atomic_position_without_defect]
        root_det_epsilon = np.sqrt(np.linalg.det(ewald.dielectric_tensor))
        epsilon_inv = np.linalg.inv(ewald.dielectric_tensor)
        cube_root_vol = math.pow(volume, 1/3)
        ewald_param = ewald.ewald_param / cube_root_vol * root_det_epsilon
        # print("ewald_param = ", ewald_param)
        diff_pot = -0.25 / volume / ewald_param ** 2  # [1/A]
        for i, r in enumerate(atomic_position_without_defect):
            # Ewald real part
            # \sum erfc(ewald*\sqrt(R*\epsilon_inv*R))
            # / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
            summation = 0
            shift =\
                defect_dft.final_structure.lattice.\
                get_cartesian_coords(r-defect_coords)
            # print("matrix, r-defect_coords, shift")
            # print(defect_dft.final_structure.lattice.matrix)
            # print(r-defect_coords)
            # print(shift)

            for v in ewald.generate_neighbor_lattices(shift=shift):
                # Skip the potential caused by the defect itself
                # r = R - atomic_pos_wrt_defect
                if np.linalg.norm(v) < 1e-8:
                    continue
                root_r_inv_epsilon_r = \
                    np.sqrt(reduce(np.dot, [v.T, epsilon_inv, v]))
                summation += \
                    scipy.special.erfc(ewald_param *
                                       root_r_inv_epsilon_r) / \
                    root_r_inv_epsilon_r
            real_part = summation / (4 * np.pi * root_det_epsilon)
            # Ewald reciprocal part
            # \sum exp(-g*\epsilon*g/(4*ewald**2)) / g*\epsilon*g [1/A]
            summation = 0
            for g in ewald.generate_neighbor_lattices(
                    include_self=False, is_reciprocal=True):
                g_epsilon_g = reduce(np.dot, [g.T, dielectric_tensor, g])
                summation += \
                    np.exp(- g_epsilon_g / 4.0 / ewald_param ** 2)\
                    / g_epsilon_g * np.cos(np.dot(g, r))  # [A^2]
            reciprocal_part = summation / volume
            model_pot[i] = (real_part + reciprocal_part + diff_pot) * coeff
            # print("atom index = ", i)
            # print("atomic_pos = ", r)
            # print("shift = ", shift)
            # print("real_part = ", real_part)
            # print("reciprocal_part = ", reciprocal_part)
            # print("diff_pot = ", diff_pot)
            # print("coeff = ", coeff)
            # print("")

        # defect site
        # TODO: Can be this part included above loop?

        # Ewald real part
        # \sum erfc(ewald*\sqrt(R*\epsilon_inv*R))
        #              / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
        root_det_epsilon = np.sqrt(np.linalg.det(ewald.dielectric_tensor))
        epsilon_inv = np.linalg.inv(ewald.dielectric_tensor)
        summation = 0
        for v in ewald.generate_neighbor_lattices():
            # Skip the potential caused by the defect itself
            # r = R - atomic_pos_wrt_defect
            if np.linalg.norm(v) < 1e-8:
                continue
            root_r_inv_epsilon_r = \
                np.sqrt(reduce(np.dot, [v.T, epsilon_inv, v]))
            summation += \
                scipy.special.erfc(ewald_param *
                                   root_r_inv_epsilon_r) / \
                root_r_inv_epsilon_r
        real_part = summation / (4 * np.pi * root_det_epsilon)

        # Ewald reciprocal part
        # \sum exp(-g*\epsilon*g/(4*ewald**2)) / g*\epsilon*g [1/A]
        summation = 0
        for g in ewald.generate_neighbor_lattices(
                include_self=False, is_reciprocal=True):
            g_epsilon_g = reduce(np.dot, [g.T, dielectric_tensor, g])
            summation += \
                np.exp(- g_epsilon_g / 4.0 / ewald_param ** 2) /\
                g_epsilon_g * np.cos(np.dot(g, np.zeros(3)))  # [A^2]
        reciprocal_part = summation / volume

        # self potential
        det_epsilon = np.linalg.det(ewald.dielectric_tensor)
        self_pot =\
            - ewald_param / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))

        # print("atom index = at site")
        # print("real_part = ", real_part)
        # print("reciprocal_part = ", reciprocal_part)
        # print("diff_pot = ", diff_pot)
        # print("self_pot = ", self_pot)
        # print("coeff = ", coeff)
        # print("")

        model_pot_defect_site \
            = (real_part + reciprocal_part + diff_pot + self_pot) * coeff
        lattice_energy = model_pot_defect_site * charge / 2
        # return model_pot, model_pot_defect_site, lattice_energy

        # print("model pot on site = {0}".format(model_pot_defect_site))
        # print("lattice energy = {0}".format(lattice_energy))
        # calc ave_pot_diff
        distance_threshold = calc_max_sphere_radius(axis)
        pot_diff = []

        # error check just in case (should be removed in the future)
        if not len(distances_from_defect) == len(diff_ep) == len(model_pot):
            raise IndexError(
                "Lengths of distances_from_defect({0}), diff_ep({1}), "
                "model_pot differs({2}) differ.".
                format(len(distances_from_defect), len(diff_ep), len(model_pot)))

        for (d, a, m) in zip(distances_from_defect, diff_ep, model_pot):
            if d > distance_threshold:
                pot_diff.append(a - m)
        ave_pot_diff = float(np.mean(pot_diff))
        # print("potential difference = {0}".format(ave_pot_diff))
        alignment = -ave_pot_diff * charge
        # print("alignment-like term = {0}".format(alignment))
        # return alignment
        return cls(CorrectionMethod.extended_fnv,
                   ewald, lattice_energy, ave_pot_diff, alignment,
                   symbols_without_defect, distances_from_defect, diff_ep,
                   model_pot)

