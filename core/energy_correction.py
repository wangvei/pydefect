#!/usr/bin/env python

import argparse
import glob
import json
import numpy as np
import scipy
import scipy.constants as sconst
from functools import reduce
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from pydefect.taka_my_classes.perfect import Unitcell, Supercell
from pydefect.input_maker.input_maker import Defect

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "December 8, 2017"


def correct_energy(perfect_unitcell,
                   perfect_supercell,
                   defect_input,
                   defect_raw_energy,
                   defect_electrostatic_potential,
                   defect_structure,
                   method = "ExtendedFNV"):
    """
    Args:
        perfect_unitcell (Unitcell):
        perfect_supercell (Supercell)
        defect_input (Defect):
        defect_raw_energy (float):
        defect_electrostatic_potential (list of float):
        defect_structure (Structure):
        method (str): method to correct energy.
                      Now "no_correction" and "extended_FNV" is available.

    Returns (float): Corrected energy.
    """
    implemented_method = ["extended_FNV", "no_correction"]
    if method not in implemented_method:
        raise ValueError("{} method is not implemented.".format(method))
    alignment = "Not defined"
    if method == "no_correction":
        alignment = 0
    elif method == "extended_FNV":
        alignment = \
            _compute_alignment_by_extended_fnv(perfect_unitcell,
                                               perfect_supercell,
                                               defect_input,
                                               defect_electrostatic_potential,
                                               defect_structure)
    return defect_raw_energy + alignment


def calc_ewald_real_pot(ewald, dielectric_tensor,
                        set_r_vectors):
    """
    \sum erfc(ewald*\sqrt(R*\epsilon_inv*R)) 
                 / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(dielectric_tensor))
    epsilon_inv = np.linalg.inv(dielectric_tensor)
    each = np.zeros(len(set_r_vectors))
    for i, R in enumerate(set_r_vectors):
        # Skip the potential caused by the defect itself
        # r = R - atomic_pos_wrt_defect
        if np.linalg.norm(R) < 1e-8: 
            continue
        root_r_inv_epsilon_r = np.sqrt(reduce(np.dot, [R.T, epsilon_inv, R]))
        each[i] = scipy.special.erfc(ewald * root_r_inv_epsilon_r) /\
            root_r_inv_epsilon_r
    return np.sum(each) / (4 * np.pi * root_det_epsilon)


def calc_ewald_reciprocal_potential(ewald,
                                    dielectric_tensor,
                                    set_g_vectors,
                                    volume,
                                    atomic_pos_wrt_defect=np.zeros(3)):
    """
    \sum exp(-g*\epsilon*g/(4*ewald**2)) / g*\epsilon*g [1/A]
    """
    each = np.zeros(len(set_g_vectors))
    for i, g in enumerate(set_g_vectors):
        g_epsilon_g = reduce(np.dot, [g.T, dielectric_tensor, g])  # [1/A^2]
        each[i] = np.exp(- g_epsilon_g / 4.0 / ewald ** 2) / \
            g_epsilon_g * np.cos(np.dot(g, atomic_pos_wrt_defect))  # [A^2]
    return np.sum(each) / volume


def calc_ewald_self_potential(ewald, dielectric_tensor):  # [1/A]
    det_epsilon = np.linalg.det(dielectric_tensor)
    return -ewald / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))


def calc_ewald_diff_potential(ewald, volume):
    return -0.25 / volume / ewald ** 2  # [1/A]


def calc_model_pot_and_lat_energy(ewald,
                                  charge,
                                  atomic_pos_wo_defect,
                                  defect_pos,
                                  dielectric_tensor,
                                  volume,
                                  set_R_vectors,
                                  set_G_vectors,
                                  axis):
    atomic_pos_wrt_defect \
        = [np.dot(axis, v - defect_pos) for v in atomic_pos_wo_defect]
    coeff = charge * sconst.elementary_charge * \
        1.0e10 / sconst.epsilon_0  # [V]
    model_pot = [None for _ in atomic_pos_wo_defect]
    for i, r in enumerate(atomic_pos_wo_defect):
        set_r_vectors = set_R_vectors - atomic_pos_wrt_defect[i]
        real_part\
            = calc_ewald_real_pot(ewald,
                                  dielectric_tensor,
                                  set_r_vectors)
        reciprocal_part \
            = calc_ewald_reciprocal_potential(ewald,
                                              dielectric_tensor,
                                              set_G_vectors,
                                              volume,
                                              atomic_pos_wrt_defect = r)
        diff_pot = calc_ewald_diff_potential(ewald, volume)
        model_pot[i] \
            = (real_part + reciprocal_part + diff_pot) * coeff
    real_part = calc_ewald_real_pot(ewald, 
                                    dielectric_tensor, 
                                    set_R_vectors)
    reciprocal_part = calc_ewald_reciprocal_potential(ewald,
                                                      dielectric_tensor,
                                                      set_G_vectors,
                                                      volume,
                                                      np.zeros(3))
    diff_pot = calc_ewald_diff_potential(ewald, volume)
    self_pot = calc_ewald_self_potential(ewald, dielectric_tensor)
    model_pot_defect_site \
        = (real_part + reciprocal_part + diff_pot + self_pot) * coeff
    lattice_energy = model_pot_defect_site * charge / 2
    return model_pot, model_pot_defect_site, lattice_energy


def get_distance_two_planes(lattice_vectors):
    # (a_i \times a_j) \ddot a_k / |a_i \times  a_j| 
    distance = np.zeros(3, dtype = float)
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
                                      axis):
    distance_threshold = get_max_sphere_radius(axis)
    pot_diff = []
    for (d, a, m) in zip(distance, abinitio_pot, model_pot):
        if d > distance_threshold:
            pot_diff.append(a - m)
    return np.mean(pot_diff)


def _compute_alignment_by_extended_fnv(perfect_unitcell,
                                       perfect_supercell,
                                       defect_input,
                                       defect_electrostatic_potential,
                                       defect_structure):

    """
    Corrects error of formation energy of point defect due to finite-supercell
    effect.
    Args:
        perfect_unitcell (Unitcell):
        perfect_supercell (Supercell):
        defect_input (Defect):
        defect_electrostatic_potential (list of float):
        defect_structure (Structure):

    Returns (float): Corrected energy by extended FNV method.

    """
    axis = np.array(perfect_unitcell.structure.lattice.matrix)
    defect_coords = [site.coords for site in defect_structure.sites]
    volume = defect_structure.lattice.volume

    # read defect_pos (fractional coordination)
    # TODO: implement after complete Defect class
    read_dpos = defect_dict["defect_position"]
    if isinstance(read_dpos, int):
        defect_pos = np.array(defect_coords[read_dpos - 1])
    elif isinstance(read_dpos, list):
        defect_pos = np.array(read_dpos)
    else:
        raise TypeError("Could not read defect_pos.")
    ref_pot = np.array(correct_dict["reference_potential"])
    if len(defect_electrostatic_potential) == len(ref_pot):
        defect_index = read_dpos
        atomic_position_without_defect =\
            np.delete(defect_coords, defect_index - 1, 0)
    elif len(defect_electrostatic_potential) == len(ref_pot) + 1:
        defect_index = read_dpos
        defect_electrostatic_potential =\
            np.delete(defect_electrostatic_potential, defect_index - 1)
        atomic_position_without_defect =\
            np.delete(defect_coords, defect_index - 1, 0)
    elif len(defect_electrostatic_potential) == len(ref_pot) - 1:
        defect_index = 192  # in outcar, index 1~, temporary hard coding
        ref_pot = np.delete(ref_pot, defect_index - 1)
        print("warning: defect index should be automatic implement")
        atomic_position_without_defect = defect_coords
    else:
        raise ValueError("This code can not applied to more than one defect.")

    diff_pot = -(defect_electrostatic_potential - ref_pot)

    distances_from_defect = np.array(defect_dict["distance_from_defect"])

    dielectric_tensor = perfect_unitcell.dielectric_tensor

    charge = defect_input.charge
    # TODO: check ewald or ewald_param?
    ewald = perfect_unitcell.get_ewald_param()
    set_R_vectors = np.array(correct_json["set_R_vector"])
    set_G_vectors = np.array(correct_json["set_G_vector"])
    # potential.sh 3
    model_pot, model_pot_site, lattice_energy \
        = calc_model_pot_and_lat_energy(ewald,
                                        charge,
                                        atomic_position_without_defect,
                                        defect_pos,
                                        dielectric_tensor,
                                        volume,
                                        set_R_vectors,
                                        set_G_vectors,
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

