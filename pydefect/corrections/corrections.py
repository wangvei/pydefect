#!/usr/bin/env python

from abc import ABC, abstractmethod
from copy import deepcopy
from functools import reduce
from itertools import product, groupby
import json
from math import sqrt, pow
from operator import itemgetter

import matplotlib.pyplot as plt

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core import Structure

from scipy.special import erfc
from scipy.constants import elementary_charge, epsilon_0
from scipy.stats import mstats

from pydefect.util.structure_tools import defect_center, distance_list
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults

"""
This module provides functions used to correct defect formation energies
with finite supercell-size dependencies.
"""

__author__ = "Akira Takahashi, Yu Kumagai"
__maintainer__ = "Akira Takahashi, Yu Kumagai"


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


# def calc_neighbor_lattices(max_norm, lattice):
#     lattices = []
#     max_int = [int(max_norm / np.linalg.norm(lattice[i])) + 1
#                for i in range(3)]
#     for index in product(range(-max_int[0], max_int[0] + 1),
#                          range(-max_int[1], max_int[1] + 1),
#                          range(-max_int[2], max_int[2] + 1)):
#         cart_vector = np.dot(lattice.transpose(), np.array(index))
#         if np.linalg.norm(cart_vector) < max_norm:
#             lattices.append(cart_vector)

# return lattices


# for searching ewald parameter
def create_neighbor_lattices(lattice_vectors, max_length):
    """ Return a set of lattice vectors within the max length.

    Note that angles between any two axes are assumed to be between 60 and
    120 deg.

    Args:
        lattice_vectors (np.ndarray): 3x3 matrix.
        max_length (float): Max length to search lattice set.

    Return (np.ndarray): Cartesian vector of lattice points.
    """
    max_int = [int(max_length / np.linalg.norm(lattice_vectors[i])) + 1
               for i in range(3)]
    neighbor_lattices = []
    for index in product(range(-max_int[0], max_int[0] + 1),
                         range(-max_int[1], max_int[1] + 1),
                         range(-max_int[2], max_int[2] + 1)):
        cart_vector = np.dot(lattice_vectors.transpose(), np.array(index))
        if np.linalg.norm(cart_vector) < max_length:
            neighbor_lattices.append(cart_vector.tolist())

    return neighbor_lattices


class Ewald(MSONable):

    def __init__(self,
                 lattice: Lattice,
                 dielectric_tensor: np.array,
                 ewald_param: float,
                 prod_cutoff_fwhm: float,
                 real_neighbor_lattices: list,
                 reciprocal_neighbor_lattices: list):
        """
        Args:
            lattice (Lattice):
            dielectric_tensor (3x3 np.array):
            ewald_param (float):
            prod_cutoff_fwhm (float):
            real_neighbor_lattices (list):
            reciprocal_neighbor_lattices (list):
        """

        self.lattice = lattice
        self.volume = lattice.volume
        self.reciprocal_lattice_matrix = \
            self.lattice.reciprocal_lattice.matrix
        self.dielectric_tensor = dielectric_tensor
        self.ewald_param = ewald_param
        self.prod_cutoff_fwhm = prod_cutoff_fwhm

        # include self lattice site
        self.real_neighbor_lattices = list(real_neighbor_lattices)
        self.reciprocal_neighbor_lattices = list(reciprocal_neighbor_lattices)

        # self.real_neighbor_lattices = \
        #     calc_neighbor_lattices(self.max_r_vector_norm,
        #                            self.lattice)

        # self.reciprocal_neighbor_lattices = \
        #     calc_neighbor_lattices(self.max_g_vector_norm,
        #                            self.reciprocal_lattice_matrix)

    @property
    def max_r_vector_norm(self):
        root_det_dielectric = sqrt(np.linalg.det(self.dielectric_tensor))
        cube_root_vol = pow(self.volume, 1 / 3)
        ewald = self.ewald_param / cube_root_vol * root_det_dielectric

        return self.prod_cutoff_fwhm / ewald

    @property
    def max_g_vector_norm(self):
        root_det_dielectric = sqrt(np.linalg.det(self.dielectric_tensor))
        cube_root_vol = pow(self.volume, 1 / 3)
        ewald = self.ewald_param / cube_root_vol * root_det_dielectric

        return 2 * ewald * self.prod_cutoff_fwhm

    def to_json_file(self, filename: str):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def from_optimization(cls,
                          structure: Structure,
                          dielectric_tensor: np.array,
                          initial_ewald_param: float = None,
                          convergence: float = 1.05,
                          prod_cutoff_fwhm: float = 20.0):
        """ Get optimized ewald parameter.

        Args:
            structure (Structure):
                Structure
            dielectric_tensor (numpy 3x3 array):
                dielectric tensor
            initial_ewald_param (float):
                Initial guess of parameter.
            convergence (float):
                If 1/convergence < n_(real)/n_(reciprocal) < convergence,
                where n_(real) and n_(reciprocal) is number of real lattices
                and reciprocal lattices, finishes optimization and
                returns ewald_param.
            prod_cutoff_fwhm (float):
                product of cutoff radius of G-vector and gaussian FWHM.
                Increasing this value, calculation will be more accurate, but
                slower.

        Returns:
            Optimized ewald_param (float).
        """
        root_det_dielectric = sqrt(np.linalg.det(dielectric_tensor))
        real_lattice = structure.lattice.matrix
        reciprocal_lattice = structure.lattice.reciprocal_lattice.matrix
        cube_root_vol = pow(structure.lattice.volume, 1 / 3)

        if initial_ewald_param is not None:
            ewald_param = initial_ewald_param
        else:
            # determine initial ewald parameter to satisfy following:
            # max_int(Real) = max_int(Reciprocal)
            # in neighbor_lattices function.

            # Left term:
            # max_int(Real) = 2 * x * Y  / l_r where x, Y, and l_r are ewald,
            # prod_cutoff_fwhm, and axis length of real lattice, respectively.

            # Right term:
            # max_int(reciprocal) = Y  / (x * l_g)
            # where l_g is axis length of reciprocal lattice, respectively.
            # Then, x = sqrt(l_g / l_r / 2)
            # gmean : geometric mean,  (a1 * a2 * a3)^(1/3)

            l_r = mstats.gmean([np.linalg.norm(v) for v in real_lattice])
            l_g = mstats.gmean([np.linalg.norm(v) for v in reciprocal_lattice])
            ewald_param \
                = np.sqrt(l_g / l_r / 2) * cube_root_vol / root_det_dielectric

        # TODO: add smaller prod_cutoff_fwhm for heavier calculations
        for i in range(10):
            ewald = ewald_param / cube_root_vol * root_det_dielectric
            # count number of real lattice
            max_r_vector_norm = prod_cutoff_fwhm / ewald
            real_neighboring_lattices = \
                create_neighbor_lattices(real_lattice, max_r_vector_norm)

            # count number of reciprocal lattice
            max_g_vector_norm = 2 * ewald * prod_cutoff_fwhm
            reciprocal_neighboring_lattices = \
                create_neighbor_lattices(reciprocal_lattice,
                                         max_g_vector_norm)

            real_to_reciprocal_ratio = len(real_neighboring_lattices) \
                                       / len(reciprocal_neighboring_lattices)

            if 1 / convergence < real_to_reciprocal_ratio < convergence:
                return cls(structure.lattice, dielectric_tensor, ewald_param,
                           prod_cutoff_fwhm, real_neighboring_lattices,
                           reciprocal_neighboring_lattices)
            else:
                ewald_param *= real_to_reciprocal_ratio ** 0.17
        else:
            raise ValueError("The initial ewald param may not be adequate.")

    def neighbor_lattices(self,
                          include_self: bool = True,
                          shift: np.array = np.zeros(3),
                          is_reciprocal_space: bool = False):
        """
        Args:
             include_self (bool):
                Whether to include the central lattice point itself.
             shift (np.array):
             is_reciprocal_space (bool):
        """

        if is_reciprocal_space:
            lattices = np.array(deepcopy(self.reciprocal_neighbor_lattices))
        else:
            lattices = np.array(deepcopy(self.real_neighbor_lattices))

        if include_self is False:
            np.delete(lattices, int(len(lattices) - 1) / 2)

        return lattices + shift


class Correction(ABC):
    @abstractmethod
    def correction_energy(self):
        pass

    # def manually_added_correction_energy(self):
    #     return


# class PointCharge:

# def __init__(self, ewald, lattice_energy): pass


class ExtendedFnvCorrection(Correction, MSONable):

    def __init__(self, ewald, lattice_energy, diff_ave_pot,
                 alignment_correction_energy, symbols_without_defect,
                 distances_from_defect, difference_electrostatic_pot, model_pot,
                 manually_added_correction_energy=0):
        """
        Args:
            ewald (Ewald):
            lattice_energy (float):
            diff_ave_pot (float):
            alignment_correction_energy (float):
            symbols_without_defect (list of str):
            distances_from_defect (list of float):
            model_pot (list of float):
            difference_electrostatic_pot (list of float):
            manually_added_correction_energy (float or None):
        """

        # error check just in case (should be removed in the future)
        if not len(symbols_without_defect) == len(distances_from_defect) == \
               len(difference_electrostatic_pot) == len(model_pot):
            raise IndexError(
                "Lengths of symbols({}), distances({}), "
                "electro_static_pot({}), model_pot differs({}) differ.".
                    format(len(symbols_without_defect),
                           len(distances_from_defect),
                           len(difference_electrostatic_pot),
                           len(model_pot)))
        self.ewald = ewald
        self.lattice_energy = lattice_energy
        self.diff_ave_pot = diff_ave_pot
        self.alignment_correction_energy = alignment_correction_energy
        self.symbols_without_defect = symbols_without_defect
        self.distances_from_defect = list(distances_from_defect)
        self.difference_electrostatic_pot = list(difference_electrostatic_pot)
        self.model_pot = list(model_pot)
        self.manually_added_correction_energy = manually_added_correction_energy

    @property
    def point_charge_correction_energy(self):
        return - self.lattice_energy

    @property
    def correction_energy(self):
        return self.point_charge_correction_energy \
               + self.alignment_correction_energy \
               + self.manually_added_correction_energy

    @property
    def max_sphere_radius(self):
        lattice_vectors = self.ewald.lattice
        return calc_max_sphere_radius(lattice_vectors)

    # TODO: if remote connect and can't be displayed figure by X,
    # matplotlib should raise exception.
    def plot_distance_vs_potential(self, file_name=None):
        property_without_defect = list(zip(self.symbols_without_defect,
                                           self.distances_from_defect,
                                           self.difference_electrostatic_pot))
        property_without_defect = sorted(property_without_defect,
                                         key=itemgetter(0))
        # points_dictionary is like {"Mg": [-0.22, -0.7,...], # "O":[...]}
        points_dictionary = {}
        for k, g in groupby(property_without_defect, key=itemgetter(0)):
            values = [(x, y) for _, x, y in g]
            points_dictionary[k] = values
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # elemental electrostatic potential
        for i, (symbol, points) in enumerate(points_dictionary.items()):
            x_set = np.array([point[0] for point in points])
            y_set = np.array([point[1] for point in points])
            # set color
            gradation = i / len(points_dictionary.items())
            color = tuple(np.array([0, gradation, 1 - gradation]))
            ax.scatter(x_set, y_set, c=color, marker="x", label=symbol)
        # model potential
        ax.scatter(self.distances_from_defect, self.model_pot,
                   c=(1, 0, 0), marker=".", label="model potential")
        # difference between model potential and electrostatic potential
        diff_model_electrostatic = \
            np.array(self.difference_electrostatic_pot) - \
            np.array(self.model_pot)
        ax.scatter(self.distances_from_defect, diff_model_electrostatic,
                   c=(1, 1, 1), marker="o", label="model - electrostatic",
                   linewidth="0.6", edgecolors="black"
                   )
        # potential difference
        point_x = [self.max_sphere_radius, max(self.distances_from_defect)]
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

    def to_json_file(self, filename):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    @classmethod
    def compute_correction(cls, defect_entry, defect_dft, perfect_dft,
                           unitcell_dft, ewald=None):
        """
        Args
            defect_entry(DefectEntry):
            defect_dft(SupercellCalcResults):
            perfect_dft(SupercellCalcResults):
            unitcell_dft(UnitcellCalcResults):
            ewald(Ewald):
        """
        # search ewald
        if ewald is None:
            ewald = Ewald.from_optimization(
                perfect_dft.final_structure,
                unitcell_dft.total_dielectric_tensor)
        return cls.compute_alignment_by_extended_fnv(defect_entry,
                                                     defect_dft,
                                                     perfect_dft,
                                                     unitcell_dft,
                                                     ewald)

    @classmethod
    def compute_alignment_by_extended_fnv(cls,
                                          defect_entry: DefectEntry,
                                          defect_dft: SupercellCalcResults,
                                          perfect_dft: SupercellCalcResults,
                                          unitcell_dft: UnitcellCalcResults,
                                          ewald: Ewald):

        """
        Estimate correction energy of point defect formation energy calculated
        using finite-size supercell.

        Args:
            defect_entry (DefectEntry):
            defect_dft (SupercellCalcResults):
            perfect_dft (SupercellCalcResults):
            unitcell_dft (UnitcellCalcResults):
            ewald (Ewald):

        """

        dielectric_tensor = unitcell_dft.total_dielectric_tensor
        perfect_structure = perfect_dft.final_structure

        diff_ep = [-ep for ep in
                   defect_dft.relative_potential(perfect_dft, defect_entry)
                   if ep is not None]

        atomic_position_without_defect = \
            [defect_dft.final_structure.frac_coords[i]
             for i, j in enumerate(defect_entry.atom_mapping_to_perfect)
             if j is not None]

        symbols_without_defect = \
            [defect_dft.final_structure.sites[i].specie.symbol
             for i, j in enumerate(defect_entry.atom_mapping_to_perfect)
             if j is not None]

        charge = defect_entry.charge
        volume = perfect_structure.lattice.volume
        defect_coords = defect_center(defect_entry, defect_dft.final_structure)
        distances_from_defect = \
            [distance_list(defect_dft.final_structure, defect_coords)[i]
             for i, j in enumerate(defect_entry.atom_mapping_to_perfect)
             if j is not None]

        # TODO: check ewald or ewald_param?
        # model potential and lattice energy
        coeff = charge * elementary_charge * 1e10 / epsilon_0  # [V]
        model_pot = []
        root_det_epsilon = np.sqrt(np.linalg.det(ewald.dielectric_tensor))
        epsilon_inv = np.linalg.inv(ewald.dielectric_tensor)
        cube_root_vol = pow(volume, 1 / 3)
        ewald_param = ewald.ewald_param / cube_root_vol * root_det_epsilon
        diff_pot = -0.25 / volume / ewald_param ** 2  # [1/A]

        for r in atomic_position_without_defect:
            # Ewald real part
            # \sum erfc(ewald*\sqrt(R*\epsilon_inv*R))
            # / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
            summation = 0
            shift = defect_dft.final_structure.lattice. \
                get_cartesian_coords(r - defect_coords)

            for v in ewald.neighbor_lattices(shift=shift):
                # Skip the potential caused by the defect itself
                # r = R - atomic_pos_wrt_defect
                if np.linalg.norm(v) < 1e-8:
                    continue
                root_r_inv_epsilon_r = \
                    np.sqrt(reduce(np.dot, [v.T, epsilon_inv, v]))
                summation += erfc(ewald_param * root_r_inv_epsilon_r) / \
                             root_r_inv_epsilon_r
            real_part = summation / (4 * np.pi * root_det_epsilon)

            # Ewald reciprocal part
            # \sum exp(-g*\epsilon*g/(4*ewald**2)) / g*\epsilon*g [1/A]
            summation = 0
            for g in ewald.neighbor_lattices(
                    include_self=False, is_reciprocal_space=True):
                g_epsilon_g = reduce(np.dot, [g.T, dielectric_tensor, g])
                summation += \
                    np.exp(- g_epsilon_g / 4.0 / ewald_param ** 2) \
                    / g_epsilon_g * np.cos(np.dot(g, r))  # [A^2]

            reciprocal_part = summation / volume

            model_pot.append((real_part + reciprocal_part + diff_pot) * coeff)

        # Madelung potential energy
        # TODO: Can be this part included above loop?

        # Ewald real part
        # \sum erfc(ewald*\sqrt(R*\epsilon_inv*R))
        #              / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
        root_det_epsilon = np.sqrt(np.linalg.det(ewald.dielectric_tensor))
        epsilon_inv = np.linalg.inv(ewald.dielectric_tensor)
        summation = 0
        for v in ewald.neighbor_lattices():
            # Skip the potential caused by the defect itself
            # r = R - atomic_pos_wrt_defect
            if np.linalg.norm(v) < 1e-8:
                continue
            root_r_inv_epsilon_r = \
                np.sqrt(reduce(np.dot, [v.T, epsilon_inv, v]))
            summation += erfc(ewald_param * root_r_inv_epsilon_r) / \
                         root_r_inv_epsilon_r
        real_part = summation / (4 * np.pi * root_det_epsilon)

        # Ewald reciprocal part
        # \sum exp(-g*\epsilon*g/(4*ewald**2)) / g*\epsilon*g [1/A]
        summation = 0
        for g in ewald.neighbor_lattices(
                include_self=False, is_reciprocal_space=True):
            g_epsilon_g = reduce(np.dot, [g.T, dielectric_tensor, g])
            summation += \
                np.exp(- g_epsilon_g / 4.0 / ewald_param ** 2) / \
                g_epsilon_g * np.cos(np.dot(g, np.zeros(3)))  # [A^2]
        reciprocal_part = summation / volume

        # self potential
        det_epsilon = np.linalg.det(ewald.dielectric_tensor)
        self_pot = - ewald_param / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))

        model_pot_defect_site = \
            (real_part + reciprocal_part + diff_pot + self_pot) * coeff
        lattice_energy = model_pot_defect_site * charge / 2

        # calc ave_pot_diff
        distance_threshold = \
            calc_max_sphere_radius(np.array(perfect_structure.lattice.matrix))
        pot_diff = []

        # error check just in case (should be removed in the future)
        if not len(distances_from_defect) == len(diff_ep) == len(model_pot):
            raise IndexError(
                "Lengths of distances_from_defect({0}), diff_ep({1}), "
                "model_pot differs({2}) differ.".
                    format(len(distances_from_defect), len(diff_ep),
                           len(model_pot)))

        for (d, a, m) in zip(distances_from_defect, diff_ep, model_pot):
            if d > distance_threshold:
                pot_diff.append(a - m)
        ave_pot_diff = float(np.mean(pot_diff))
        alignment = -ave_pot_diff * charge

        return cls(ewald, lattice_energy, ave_pot_diff, alignment,
                   symbols_without_defect, distances_from_defect, diff_ep,
                   model_pot)
