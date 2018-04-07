#!/usr/bin/env python

import math
import numpy as np
import scipy
import scipy.constants as sconst
from scipy.stats import mstats
from functools import reduce
from itertools import product
from pymatgen.core.structure import Structure

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


# search ewald parameter
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


def compute_ewald_param(dielectric_tensor,
                        structure,
                        initial_value=None,
                        convergence=1.05,
                        prod_cutoff_fwhm=25.0):
    """
    Get optimized ewald parameter.
    Once optimized parameter is calculated (usually slow),
    the value is cached and from next calculation,
    return the cached value if you don't specify
    computes_again = True.
    Args:
        dielectric_tensor (np.ndarray) : 3x3 matrix
        structure (Structure): pymatgen structure
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
    reciprocal_lattice = structure.lattice.reciprocal_lattice.matrix
    cube_root_vol = math.pow(structure.lattice.volume, 1 / 3)
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
        num_real_lattice = 0
        max_r_vector_norm = prod_cutoff_fwhm / ewald
        for _ in _generate_neighbor_lattices(real_lattice,
                                             max_r_vector_norm,
                                             include_self=True):
            num_real_lattice += 1
        # count number of reciprocal lattice
        num_reciprocal_lattice = 0
        max_g_vector_norm = 2 * ewald * prod_cutoff_fwhm
        for _ in _generate_neighbor_lattices(reciprocal_lattice,
                                             max_g_vector_norm,
                                             include_self=False):
            num_reciprocal_lattice += 1
        diff_real_reciprocal = num_real_lattice / num_reciprocal_lattice
        if 1 / convergence < diff_real_reciprocal < convergence:
            return ewald
        else:
            ewald_param *= diff_real_reciprocal ** 0.17


def calc_ewald_real_pot(ewald, dielectric_tensor,
                        real_lattice_vector, max_length, shift):
    """
    \sum erfc(ewald*\sqrt(R*\epsilon_inv*R)) 
                 / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(dielectric_tensor))
    epsilon_inv = np.linalg.inv(dielectric_tensor)
    summation = 0
    for r in _generate_neighbor_lattices(real_lattice_vector, max_length,
                                         True, shift=shift):
        # Skip the potential caused by the defect itself
        # r = R - atomic_pos_wrt_defect
        if np.linalg.norm(r) < 1e-8:
            continue
        root_r_inv_epsilon_r = np.sqrt(reduce(np.dot, [r.T, epsilon_inv, r]))
        summation += scipy.special.erfc(ewald * root_r_inv_epsilon_r) /\
            root_r_inv_epsilon_r
    return summation / (4 * np.pi * root_det_epsilon)


def _calc_ewald_reciprocal_potential(ewald,
                                     dielectric_tensor,
                                     reciprocal_lattice_vector,
                                     max_length,
                                     volume,
                                     atomic_pos_wrt_defect=np.zeros(3)):
    """
    \sum exp(-g*\epsilon*g/(4*ewald**2)) / g*\epsilon*g [1/A]
    """
    summation = 0
    for g in _generate_neighbor_lattices(reciprocal_lattice_vector,
                                         max_length, False):
        g_epsilon_g = reduce(np.dot, [g.T, dielectric_tensor, g])  # [1/A^2]
        summation += np.exp(- g_epsilon_g / 4.0 / ewald ** 2) / \
            g_epsilon_g * np.cos(np.dot(g, atomic_pos_wrt_defect))  # [A^2]
    return summation / volume


def _calc_ewald_self_potential(ewald, dielectric_tensor):  # [1/A]
    det_epsilon = np.linalg.det(dielectric_tensor)
    return -ewald / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))


def _calc_ewald_diff_potential(ewald, volume):
    return -0.25 / volume / ewald ** 2  # [1/A]


def calc_model_pot_and_lat_energy(ewald,
                                  charge,
                                  atomic_pos_wo_defect,
                                  defect_pos,
                                  dielectric_tensor,
                                  volume,
                                  set_R_vectors,
                                  set_G_vectors,
                                  real_lattice_vector,
                                  reciprocal_lattice_vector,
                                  real_max_length,
                                  reciprocal_max_length):
    atomic_pos_wrt_defect \
        = [np.dot(real_lattice_vector, v - defect_pos)
           for v in atomic_pos_wo_defect]
    coeff = charge * sconst.elementary_charge * \
        1.0e10 / sconst.epsilon_0  # [V]
    model_pot = [None for _ in atomic_pos_wo_defect]
    for i, r in enumerate(atomic_pos_wo_defect):
        set_r_vectors = set_R_vectors - atomic_pos_wrt_defect[i]
        real_part\
            = calc_ewald_real_pot(ewald,
                                  dielectric_tensor,
                                  real_lattice_vector,
                                  real_max_length,
                                  shift=-atomic_pos_wrt_defect[i])
        reciprocal_part \
            = _calc_ewald_reciprocal_potential(ewald,
                                               dielectric_tensor,
                                               reciprocal_lattice_vector,
                                               reciprocal_max_length,
                                               volume,
                                               atomic_pos_wrt_defect=r)
        diff_pot = _calc_ewald_diff_potential(ewald, volume)
        model_pot[i] \
            = (real_part + reciprocal_part + diff_pot) * coeff
    real_part = calc_ewald_real_pot(ewald, 
                                    dielectric_tensor,
                                    real_lattice_vector,
                                    real_max_length)
    reciprocal_part = _calc_ewald_reciprocal_potential(ewald,
                                                       dielectric_tensor,
                                                       set_G_vectors,
                                                       volume,
                                                       np.zeros(3))
    diff_pot = _calc_ewald_diff_potential(ewald, volume)
    self_pot = _calc_ewald_self_potential(ewald, dielectric_tensor)
    model_pot_defect_site \
        = (real_part + reciprocal_part + diff_pot + self_pot) * coeff
    lattice_energy = model_pot_defect_site * charge / 2
    return model_pot, model_pot_defect_site, lattice_energy


def get_distance_two_planes(lattice_vectors):
    # (a_i \times a_j) \ddot a_k / |a_i \times  a_j| 
    distance = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_times_a_j = np.cross(lattice_vectors[i - 2], lattice_vectors[i - 1])
        a_k = lattice_vectors[i]
        distance[i] = abs(np.dot(a_i_times_a_j, a_k)) / \
            np.linalg.norm(a_i_times_a_j)
    return distance


def get_max_sphere_radius(lattice_vectors):
    # Maximum radius of a sphere fitting inside the unit cell.
    return max(get_distance_two_planes(lattice_vectors)) / 2.0


def calc_average_potential_difference(distance,
                                      abinitio_pot,
                                      model_pot,
                                      lattice_vectors):
    distance_threshold = get_max_sphere_radius(lattice_vectors)
    pot_diff = []
    for (d, a, m) in zip(distance, abinitio_pot, model_pot):
        if d > distance_threshold:
            pot_diff.append(a - m)
    return np.mean(pot_diff)


def compute_alignment_by_extended_fnv(dielectric_tensor,
                                      ewald_param,
                                      perfect_structure,
                                      perfect_electrostatic_potential,
                                      defect_structure,
                                      defect_electrostatic_potential,
                                      defect_pos,
                                      atomic_position_without_defect,
                                      distances_from_defect,
                                      charge):

    """
    Corrects error of formation energy of point defect due to finite-supercell
    effect.
    Args:
        dielectric_tensor (np.ndarray): 3x3 matrix
        ewald_param (float):
        perfect_structure (Structure):
        perfect_electrostatic_potential (np.ndarray):
        defect_structure (Structure):
        defect_electrostatic_potential (np.ndarray):
        defect_pos (np.ndarray):
        atomic_position_without_defect (np.ndarray):
        distances_from_defect (list of float):
        charge (int): Charge of defect

    Returns (float): Corrected energy by extended FNV method.

    """
    axis = np.array(perfect_structure.lattice.matrix)
    volume = perfect_structure.lattice.volume
    ref_pot = perfect_electrostatic_potential

    diff_pot = -(defect_electrostatic_potential - ref_pot)

    # TODO: check ewald or ewald_param?
    # potential.sh 3
    model_pot, model_pot_site, lattice_energy \
        = calc_model_pot_and_lat_energy(ewald_param,
                                        charge,
                                        atomic_position_without_defect,
                                        defect_pos,
                                        dielectric_tensor,
                                        volume,
                                        axis)
    print("model pot on site = {0}".format(model_pot_site))
    print("lattice energy = {0}".format(lattice_energy))
    ave_pot_diff\
        = calc_average_potential_difference(distances_from_defect,
                                            diff_pot,
                                            model_pot,
                                            axis)
    print("potential difference = {0}".format(ave_pot_diff))
    alignment = -1.0 * ave_pot_diff * charge
    print("alignment-like term = {0}".format(alignment))
    return alignment

