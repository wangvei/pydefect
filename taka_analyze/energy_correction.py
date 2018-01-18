#!/usr/bin/env python

import sys
import argparse
import glob
import json
import math
import numpy as np
import scipy
from functools import reduce
from itertools import product
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, get_el_sp
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun, Outcar
#from defect_property import DefectProperty
from enum import Enum, auto

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "December 8, 2017"

def __init__(self, defect_property, perfect_property):
    self.defect_property = defect_property
    self.perfect_property = perfect_property
    self.defect_index = 1 #TODO temporary magic number, but must be automatically read!!!
    if self.perfect_property.structure.lattice != self.defect_property.structure.lattice:
        sys.exit("Lattice vectors of perfect structure and those of structure with defect is different. Did you set ISIF correctly?")

@property
def defect_property(self):
    return self.defect_property

@property
def perfect_property(self):
    return self.perfect_property

def calc_distances_from_defect(self):# TODO: Check if this algorithm is consistent with previous version.
    """
    Determine the shortest distance at each atom from the image defects.
    The shortest distance between atom and the defect is searched for iteratively.
    This algorithm may not be perfect. 

    (By Yu Kumagai.)
    """

    self.lattice_vect = self.defect_property.structure.lattice.matrix
    self.volume = self.defect_property.structure.lattice.volume
    self.defect_pos = self.defect_property.structure.cart_coords[self.defect_index-1]
    atomic_pos = self.defect_property.structure.cart_coords
    self.atomic_pos_wo_defect = np.delete(atomic_pos, self.defect_index-1, 0)
#use n_atom for each symbol? (in old defect_structure.py, defect_position)
    self.min_distances = [None for i in range(len(self.atomic_pos_wo_defect))]
    candidate_distances = [[] for i in range(len(self.atomic_pos_wo_defect))]
    for i in range(len(self.atomic_pos_wo_defect)):
        for index in product(range(-1,2), range(-1,2), range(-1,2)):

            # TODO: Check if transpose is correct.
            index_vector = np.dot(np.array(index), self.lattice_vect.transpose())

            atomic_pos_wrt_defect = (self.atomic_pos_wo_defect[i] + index_vector) - self.defect_pos
            distance = np.linalg.norm(atomic_pos_wrt_defect)
            candidate_distances[i].append(distance)
    self.min_distances[i] = min(candidate_distances[i])

def make_lattice_set(lattice_vectors, max_length, include_self): 
    """
    Return a set of lattice vectors within the max length.
    Note that angles between any two axes are assumed to be between 60 and 
    120 deg.
    """
    max_int = \
    [int(max_length / np.linalg.norm(lattice_vectors[i])) + 1 for i in range(3)]
    vectors = []
    for index in product(
            range(-max_int[0], max_int[0]+1),
            range(-max_int[1], max_int[1]+1),
            range(-max_int[2], max_int[2]+1)
            ):
        if (not include_self) and index == (0, 0, 0): 
            continue
        vector = np.dot(lattice_vectors.transpose(), np.array(index))
        norm = np.linalg.norm(vector)
        if norm < max_length: 
            vectors.append(vector)
    return np.array(vectors)

def determine_ewald_param(def_structure,
                          root_det_dielectric,
                          ):
    INIT_EWALD_PARAM = 0.010111311355097064
    PROD_CUTOFF_FWHM = 25.0 #product of cutoff radius of G-vector and gaussian FWHM.
                            #increasing this value, both accuracy and computational cost will be increased
    real_lattice = def_structure.lattice.matrix
    volume = def_structure.lattice.volume
    reciprocal_lattice_vect = def_structure.lattice.reciprocal_lattice.matrix
    ewald_param = INIT_EWALD_PARAM
    while True:
        print("hoge")
        ewald = ewald_param / math.pow(volume, float(1)/3) * root_det_dielectric
        max_G_vector_norm = 2 * ewald * PROD_CUTOFF_FWHM
        set_G_vectors = make_lattice_set(reciprocal_lattice_vect, 
                                         max_G_vector_norm, 
                                         include_self=False)
        max_R_vector_norm = PROD_CUTOFF_FWHM / ewald
        set_R_vectors = make_lattice_set(real_lattice, 
                                         max_R_vector_norm,
                                         include_self=True)
        print(f"ewald_param= {ewald_param}")
        print(f"Number of R vectors = {len(set_R_vectors)}")
        print(f"Number of G vectors = {len(set_G_vectors)}")
        diff_real_recipro = (float(len(set_R_vectors)) / len(set_G_vectors)) 
        if (diff_real_recipro > 1/1.05) and (diff_real_recipro < 1.05):
            return ewald_param
        else:
            ewald_param *= diff_real_recipro  ** 0.17

def calc_ewald_real_pot(self, atomic_pos_wrt_defect):
    """
    \sum erfc(ewald*\sqrt(R*\epsilon_inv*R)) 
                 / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(self.dielectric_tensor))
    epsilon_inv = np.linalg.inv(self.dielectric_tensor)
    ewald = self.ewald_param / math.pow(self.volume, float(1)/3) * root_det_epsilon
    each = np.zeros(len(self.set_R_vectors))
    for i, R in enumerate(self.set_R_vectors):
        # Skip the potential caused by the defect itself
        r = R - atomic_pos_wrt_defect
        if np.linalg.norm(r) < 1e-8: 
            continue
        root_R_epsilonI_R = np.sqrt(reduce(np.dot, [r.T, epsilon_inv, r]))
        each[i] = scipy.special.erfc(ewald * root_R_epsilonI_R) \
                             / root_R_epsilonI_R

    return np.sum(each) / (4 * np.pi * root_det_epsilon)

def calc_ewald_recipro_pot(self, atomic_pos_wrt_defect=np.zeros(3)):
    """
    \sum exp(-G*\epsilon*G/(4*ewald**2)) / G*\epsilon*G [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(self.dielectric_tensor))
    ewald = self.ewald_param / math.pow(self.volume, float(1)/3) * root_det_epsilon
    each = np.zeros(len(self.set_G_vectors))

    for i, G in enumerate(self.set_G_vectors):
        G_epsilon_G = reduce(np.dot, [G.T, self.dielectric_tensor, G]) # [1/A^2]
        each[i] = np.exp(- G_epsilon_G / 4.0 / ewald ** 2) \
                   / G_epsilon_G * np.cos(np.dot(G, atomic_pos_wrt_defect)) # [A^2]

    return np.sum(each) / self.volume

def calc_ewald_self_potential(self): # [1/A]
    det_epsilon = np.linalg.det(self.dielectric_tensor)
    root_det_epsilon = np.sqrt(det_epsilon)
    ewald = self.ewald_param / math.pow(self.volume, float(1)/3) * root_det_epsilon
    return -ewald / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))

def calc_ewald_diff_potential(self):
    root_det_epsilon = np.sqrt(np.linalg.det(self.dielectric_tensor))
    ewald = self.ewald_param / math.pow(self.volume, float(1)/3) * root_det_epsilon
    return -0.25 / self.volume / ewald ** 2 # [1/A]

def calc_model_potential(self):
    self.model_pot = [None for i in self.atomic_pos_wo_defect]
    coeff = \
            self.defect_property.charge * \
            scipy.constants.elementary_charge * \
            1.0e10 / scipy.constants.epsilon_0 #[V]
    for i, R in enumerate(self.atomic_pos_wo_defect):
        real_part = self.calc_ewald_real_pot(R-self.defect_pos)
        recipro_part = self.calc_ewald_recipro_pot(R)
        diff_pot = self.calc_ewald_diff_potential()
        self.model_pot[i] =\
                (real_part + recipro_part + diff_pot)*coeff
        #print(real_part, recipro_part, diff_pot, coeff)
        #print(self.defect_property.charge, scipy.constants.elementary_charge, scipy.constants.epsilon_0)
    real_part = self.calc_ewald_real_pot(np.zeros(3))
    recipro_part = self.calc_ewald_recipro_pot(np.zeros(3))
    diff_pot = self.calc_ewald_diff_potential()
    self_pot = self.calc_ewald_self_potential()
    model_pot_defect_site = (real_part + recipro_part + diff_pot + self_pot) * coeff
    print(real_part, recipro_part, diff_pot, self_pot)
    self.lattice_energy = model_pot_defect_site * self.defect_property.charge / 2




class _DefType(Enum):
    VACANCY = auto()
    SUBSTITUTIONAL = auto()
    INTERSTITIAL = auto()

def correct_energy(dirname):
    with open(dirname + "/defect.json", 'r') as f:
        defect_json=json.load(f)
    with open("./correction.json", 'r') as f:
        correct_json=json.load(f)


    axis = defect_json["axis"]
    elements = [ Element(e_name) for e_name in defect_json["elements"] ]
    def_coord_frac = defect_json["frac_coords"]
    def_structure = Structure(axis, elements, def_coord_frac)

    #readdefect_pos (fractional coordination)
    if (isinstance(defect_json["defect_position"], int)):
        defect_pos = np.array()
    elif (isinstance(defect_json["defect_position"], list)):
        defect_pos = np.array(defect_json["defect_position"])
    else:
        raise TypeError("Could not read defect_pos.")

    ref_pot = np.array(correct_json["reference_potential"])
    def_pot = np.array(defect_json["atomic_site_pot"])
    if len(def_pot) == len(ref_pot):
        def_type = _DefType.SUBSTITUTIONAL
    elif len(def_pot) == len(ref_pot)+1:
        defect_index = defect_json["defect_position"]
        def_pot = np.delete(def_pot, defect_index - 1)
        def_type = _DefType.INTERSTITIAL
    elif len(def_pot) == len(ref_pot)-1:
        def_type = _DefType.VACANCY
        defect_index = 9 #in outcar, index 1~, temporary hard coding
        ref_pot = np.delete(ref_pot, defect_index - 1)
        print("warning: should be automatic implement")
    else:
        raise ValueError("This code can not applied to more than one defect.")
    diff_pot = -(def_pot - ref_pot)

    distances_from_defect = np.array(defect_json["distance_from_defect"])

    dielectric_tensor = np.array(correct_json["dielectric_tensor"])
    root_det_dielectric = np.sqrt(np.linalg.det(dielectric_tensor))

#potential.sh step 0 and 3
    determine_ewald_param(def_structure, root_det_dielectric)
    sys.exit("temporary exit here")
    self.calc_model_potential()
    print(self.model_pot)
    print(self.lattice_energy)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--defect_dir", dest="defect_dir", type=str,
                        help="Directry name of calculation of structure with defect.\
                              If you want to correct energy of one of defect calculations, specify with this option.\
                              Otherwise, correction will done with all results of defect calculations.\
                              ")
    opts = parser.parse_args()

    if opts.defect_dir:
        correct_energy(opts.defect_dir)
    else:
        dirs = glob.glob("./defect/*_*_*/")
        if not dirs:
            print("Warning: No directory matched name defect_*_*_*.")
        for dirname in dirs:
            correct_energy(dirname)


