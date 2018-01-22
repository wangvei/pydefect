#!/usr/bin/env python

import sys
sys.path.append("/Users/takahashi/psi/bitbucket/chempotdiag/")
import chem_pot
import argparse
import glob
import json
import math
import numpy as np
import scipy
import scipy.constants as sconst
import scipy.stats.mstats as mstats
from analyze_defect import append_to_json
from functools import reduce
from itertools import product
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, get_el_sp
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from enum import Enum

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "December 8, 2017"

def make_lattice_set(lattice_vectors, max_length, include_self): 
    """
    Return a set of lattice vectors within the max length.
    Note that angles between any two axes are assumed to be between 60 and 
    120 deg.
    """
    max_int = \
    [int(max_length / np.linalg.norm(lattice_vectors[i])) + 1 for i in range(3)]
    vectors = []
    for index in product(range(-max_int[0], max_int[0]+1),
                         range(-max_int[1], max_int[1]+1),
                         range(-max_int[2], max_int[2]+1)):
        if (not include_self) and index == (0, 0, 0): 
            continue
        vector = np.dot(lattice_vectors.transpose(), np.array(index))
        norm = np.linalg.norm(vector)
        if norm < max_length: 
            vectors.append(vector)
    return np.array(vectors)

def determine_ewald_param(def_structure,
                          root_det_dielectric):
    PROD_CUTOFF_FWHM = 25.0 #product of cutoff radius of G-vector and gaussian FWHM.
    real_lattice = def_structure.lattice.matrix
    cube_root_vol = math.pow(def_structure.lattice.volume, 1.0/3)
    reciprocal_lattice = def_structure.lattice.reciprocal_lattice.matrix
#determine inital ewald parameter to satisfy following:
# max_int(Real) = max_int(Reciprocal) in make_lattice_set function.
# Left term:
# max_int(Real) = 2 * x * Y  / l_r where x, Y, and l_r are ewald,
# PROD_CUTOFF_FWHM, and axis length of real lattice, respectively.
# Right term:
# max_int(reciprocal) = Y  / (x * l_g) where l_g is axis length of reciprocal lattice, respectively.
#Then, x = sqrt(l_g / l_r / 2)
#To calculate ewald_param (not ewald), need to consider cube_root_vol and root_det_dielectric
    l_r = mstats.gmean([np.linalg.norm(v) for v in real_lattice])
    l_g = mstats.gmean([np.linalg.norm(v) for v in reciprocal_lattice])
    ewald_param = np.sqrt( l_g / l_r / 2 ) * cube_root_vol / root_det_dielectric

    while True:
        ewald = ewald_param / cube_root_vol * root_det_dielectric
        max_G_vector_norm = 2 * ewald * PROD_CUTOFF_FWHM
        set_G_vectors = make_lattice_set(reciprocal_lattice, 
                                         max_G_vector_norm, 
                                         include_self=False)
        max_R_vector_norm = PROD_CUTOFF_FWHM / ewald
        set_R_vectors = make_lattice_set(real_lattice, 
                                         max_R_vector_norm,
                                         include_self=True)
        print("Number of R, G vectors = {0}, {1}".format(len(set_R_vectors), len(set_G_vectors)))
        diff_real_recipro = (float(len(set_R_vectors)) / len(set_G_vectors)) 
        if (diff_real_recipro > 1/1.05) and (diff_real_recipro < 1.05):
            ewald_param = ewald * cube_root_vol / root_det_dielectric
            return ewald, set_R_vectors, set_G_vectors
        else:
            ewald_param *= diff_real_recipro  ** 0.17

def calc_ewald_real_pot(ewald, dielectric_tensor,
                        set_R_vectors
                        ):
                        #atomic_pos_wrt_defect):
    """
    \sum erfc(ewald*\sqrt(R*\epsilon_inv*R)) 
                 / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(dielectric_tensor))
    epsilon_inv = np.linalg.inv(dielectric_tensor)
    each = np.zeros(len(set_R_vectors))
    for i, R in enumerate(set_R_vectors):
        # Skip the potential caused by the defect itself
        #r = R - atomic_pos_wrt_defect
        if np.linalg.norm(R) < 1e-8: 
            continue
        root_R_epsilonI_R = np.sqrt(reduce(np.dot, [R.T, epsilon_inv, R]))
        each[i] = scipy.special.erfc(ewald * root_R_epsilonI_R) \
                             / root_R_epsilonI_R

    return np.sum(each) / (4 * np.pi * root_det_epsilon)

def calc_ewald_recipro_pot(ewald,
                           dielectric_tensor,
                           set_G_vectors,
                           volume,
                           atomic_pos_wrt_defect=np.zeros(3)):
    """
    \sum exp(-G*\epsilon*G/(4*ewald**2)) / G*\epsilon*G [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(dielectric_tensor))
    each = np.zeros(len(set_G_vectors))
    for i, G in enumerate(set_G_vectors):
        G_epsilon_G = reduce(np.dot, [G.T, dielectric_tensor, G]) # [1/A^2]
        each[i] = np.exp(- G_epsilon_G / 4.0 / ewald ** 2) \
                   / G_epsilon_G * np.cos(np.dot(G, atomic_pos_wrt_defect)) # [A^2]
    return np.sum(each) / volume

def calc_ewald_self_potential(ewald, dielectric_tensor): # [1/A]
    det_epsilon = np.linalg.det(dielectric_tensor)
    root_det_epsilon = np.sqrt(det_epsilon)
    return -ewald / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))

def calc_ewald_diff_potential(ewald, volume):
    return -0.25 / volume / ewald ** 2 # [1/A]

def calc_model_pot_and_lat_energy(ewald,
                                  charge,
                                  atomic_pos_wo_defect,
                                  defect_pos,
                                  dielectric_tensor,
                                  volume,
                                  set_R_vectors,
                                  set_G_vectors,
                                  axis):
    #atomic_pos_wrt_defect = [v-defect_pos for v in atomic_pos_wo_defect] 
    atomic_pos_wrt_defect_frac = [v-defect_pos for v in atomic_pos_wo_defect] #try
    atomic_pos_wrt_defect = [np.dot(axis, v - defect_pos) for v in atomic_pos_wo_defect] #try
    coeff = charge * sconst.elementary_charge * \
            1.0e10 / sconst.epsilon_0 #[V]
    model_pot = [None for i in atomic_pos_wo_defect]
    for i, r in enumerate(atomic_pos_wo_defect):
        set_r_vectors = set_R_vectors - atomic_pos_wrt_defect[i]
        real_part\
            = calc_ewald_real_pot(ewald,
                                  dielectric_tensor,
                                  set_r_vectors)
        recipro_part \
            = calc_ewald_recipro_pot(ewald,
                                     dielectric_tensor, 
                                     set_G_vectors, 
                                     volume,
                                     atomic_pos_wrt_defect = r)
        diff_pot = calc_ewald_diff_potential(ewald, volume)
        model_pot[i] \
            = (real_part + recipro_part + diff_pot) * coeff
    real_part = calc_ewald_real_pot(ewald, 
                                    dielectric_tensor, 
                                    set_R_vectors)
                                    #np.zeros(3))
    recipro_part = calc_ewald_recipro_pot(ewald,
                                          dielectric_tensor,
                                          set_G_vectors,
                                          volume,
                                          np.zeros(3))
    diff_pot = calc_ewald_diff_potential(ewald, volume)
    self_pot = calc_ewald_self_potential(ewald, dielectric_tensor)
    model_pot_defect_site \
        = (real_part + recipro_part + diff_pot + self_pot) * coeff
    lattice_energy = model_pot_defect_site * charge / 2
    return model_pot, model_pot_defect_site, lattice_energy

def get_distance_two_planes(lattice_vectors):
    # (a_i \times a_j) \ddot a_k / |a_i \times  a_j| 
    distance = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_times_a_j = np.cross(lattice_vectors[i-2], lattice_vectors[i-1])
        a_k = lattice_vectors[i]
        distance[i] = abs(np.dot(a_i_times_a_j, a_k)) \
                         / np.linalg.norm(a_i_times_a_j)
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

class _DefType(Enum):
    VACANCY = 1
    SUBSTITUTIONAL = 2
    INTERSTITIAL = 3

def correct_energy(dirname, defect_dict, correct_dict):

    axis = defect_dict["axis"]
    elements = [ Element(e_name) for e_name in defect_dict["elements"] ]
    def_coord_frac = defect_dict["frac_coords"]
    def_structure = Structure(axis, elements, def_coord_frac)
    volume = def_structure.lattice.volume

    #read defect_pos (fractional coordination)
    read_dpos = defect_dict["defect_position"]
    if (isinstance(read_dpos, int)):
        defect_pos = np.array(def_coord_frac[read_dpos - 1])
    elif (isinstance(read_dpos, list)):
        defect_pos = np.array(read_dpos)
    else:
        raise TypeError("Could not read defect_pos.")

    ref_pot = np.array(correct_dict["reference_potential"])
    def_pot = np.array(defect_dict["atomic_site_pot"])
    if len(def_pot) == len(ref_pot):
        def_type = _DefType.SUBSTITUTIONAL
        defect_index = read_dpos
        atomic_pos_wo_defect = np.delete(def_coord_frac, defect_index-1, 0)
    elif len(def_pot) == len(ref_pot)+1:
        def_type = _DefType.INTERSTITIAL
        defect_index = read_dpos
        def_pot = np.delete(def_pot, defect_index-1)
        atomic_pos_wo_defect = np.delete(def_coord_frac, defect_index-1, 0)
    elif len(def_pot) == len(ref_pot)-1:
        def_type = _DefType.VACANCY
        defect_index = 192 #in outcar, index 1~, temporary hard coding
        ref_pot = np.delete(ref_pot, defect_index-1)
        print("warning: defect index should be automatic implement")
        atomic_pos_wo_defect = def_coord_frac
    else:
        raise ValueError("This code can not applied to more than one defect.")

    diff_pot = -(def_pot - ref_pot)

    distances_from_defect = np.array(defect_dict["distance_from_defect"])

    dielectric_tensor = np.array(correct_dict["dielectric_tensor"])\
                      + np.array(correct_dict["dielectric_ionic_tensor"])
    root_det_dielectric = np.sqrt(np.linalg.det(dielectric_tensor))

    charge = defect_dict["charge"]

    ewald = correct_json["ewald"]
    set_R_vectors = np.array(correct_json["set_R_vector"])
    set_G_vectors = np.array(correct_json["set_G_vector"])
#potential.sh 3
    model_pot, model_pot_site, lattice_energy \
        =  calc_model_pot_and_lat_energy(ewald,
                                         charge,
                                         atomic_pos_wo_defect,
                                         defect_pos,
                                         dielectric_tensor,
                                         volume,
                                         set_R_vectors,
                                         set_G_vectors,
                                         axis)
    print("model pot on site = {0}".format(model_pot_site))
    print("lattice energy = {0}".format(lattice_energy))
    ave_pot_diff\
        =calc_average_potential_difference(distances_from_defect,
                                           diff_pot,
                                           model_pot,
                                           axis)
    print("potential difference = {0}".format(ave_pot_diff))
    alignment = -1.0 * ave_pot_diff * charge
    print("alignment-like term = {0}".format(alignment))

def plot_fermi_vs_formation():


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--defect_dir", dest="defect_dir", type=str,
                        help="Directry name of calculation of structure with defect.\
                              If you want to correct energy of one of defect calculations, specify with this option.\
                              Otherwise, correction will done with all results of defect calculations.\
                              ")
    opts = parser.parse_args()

    #read files
    with open("./correction.json", 'r') as f:
        correct_json=json.load(f)

    #determine ewald parameter (old potential.sh 0)
    #axis = defect_json["axis"]
    #elements = [ Element(e_name) for e_name in defect_json["elements"] ]
    #def_coord_frac = defect_json["frac_coords"]
    #def_structure = Structure(axis, elements, def_coord_frac)
    structure_for_ewald = Poscar.from_file("defect/Va_O1_2/POSCAR-final").structure
    print("warning: structure for ewald must be extracted from correction_file (only axis is needed)" )   
    dielectric_tensor = np.array(correct_json["dielectric_tensor"]) \
                      + np.array(correct_json["dielectric_ionic_tensor"])
    root_det_dielectric = np.sqrt(np.linalg.det(dielectric_tensor))
    if not "ewald" in correct_json:
        ewald, set_r, set_g = determine_ewald_param(structure_for_ewald, root_det_dielectric)
        append_to_json("correction.json", "ewald", ewald)
        append_to_json("correction.json", "set_R_vector", set_r.tolist())
        append_to_json("correction.json", "set_G_vector", set_g.tolist())
        correct_json["ewald"] = ewald
        correct_json["set_R_vector"] = set_r
        correct_json["set_G_vector"] = set_g

    #main
    if opts.defect_dir:
        with open(opts.defect_dir, 'r') as f:
            defect_json=json.load(f)
        correct_energy(opts.defect_dir,
                       defect_json,
                       correct_json)
    else:
        dirs = glob.glob("./defect/*_*_*/")
        if not dirs:
            print("Warning: No directory matched name defect_*_*_*.")
        for dirname in dirs:
            with open(dirname + "/defect.json", 'r') as f:
                defect_json=json.load(f)
            correct_energy(dirname,
                           defect_json,
                           correct_json)

