#!/usr/bin/env python

import sys
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
from defect_property import DefectProperty

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "December 8, 2017"

class DefectCorrector:
    """
    Object for correction of defect energy. (Kumagai method)
    Args:
        defect_property  (DefectProperty)
        perfect_property (DefectProperty)
    """
    #temporary determined
    __e11 = 23.665292999999998
    __e22 = 31.206830999999998
    __e33 = 41.606976000000003
    dielectric_tensor = np.array([[__e11, 0, 0],
                                  [0, __e22, 0],
                                  [0, 0, __e33]])
    def __init__(self, defect_property, perfect_property):
        self.__defect_property = defect_property
        self.__perfect_property = perfect_property
        self.__defect_index = 1 #TODO temporary magic number, but must be automatically read!!!
        if self.perfect_property.structure.lattice != self.defect_property.structure.lattice:
            sys.exit("Lattice vectors of perfect structure and those of structure with defect is different. Did you set ISIF correctly?")

    @property
    def defect_property(self):
        return self.__defect_property

    @property
    def perfect_property(self):
        return self.__perfect_property

    def __calc_distances_from_defect(self):# TODO: Check if this algorithm is consistent with previous version.
        """
        Determine the shortest distance at each atom from the image defects.
        The shortest distance between atom and the defect is searched for iteratively.
        This algorithm may not be perfect. 
        (By Yu Kumagai.)
        """

        self.__lattice_vect = self.defect_property.structure.lattice.matrix
        self.__volume = self.defect_property.structure.lattice.volume
        self.__defect_pos = self.defect_property.structure.cart_coords[self.__defect_index-1]
        atomic_pos = self.defect_property.structure.cart_coords
        self.__atomic_pos_wo_defect = np.delete(atomic_pos, self.__defect_index-1, 0)
#use n_atom for each symbol? (in old defect_structure.py, defect_position)
        self.__min_distances = [None for i in range(len(self.__atomic_pos_wo_defect))]
        candidate_distances = [[] for i in range(len(self.__atomic_pos_wo_defect))]
        for i in range(len(self.__atomic_pos_wo_defect)):
            for index in product(range(-1,2), range(-1,2), range(-1,2)):

                # TODO: Check if transpose is correct.
                index_vector = np.dot(np.array(index), self.__lattice_vect.transpose())

                atomic_pos_wrt_defect = (self.__atomic_pos_wo_defect[i] + index_vector) - self.__defect_pos
                distance = np.linalg.norm(atomic_pos_wrt_defect)
                candidate_distances[i].append(distance)
        self.__min_distances[i] = min(candidate_distances[i])

    def __make_lattice_set(self, lattice_vectors, max_length, include_self): 
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

    def __determine_ewald_param(self):
        INIT_EWALD_PARAM = 0.010111311355097064
        PROD_CUTOFF_FWHM = 25.0 #product of cutoff radius of G-vector and gaussian FWHM.
                                #increasing this value, both accuracy and computational cost will be increased

        real_lattice = self.__lattice_vect
        reciprocal_lattice_vect = \
        self.perfect_property.structure.lattice.reciprocal_lattice.matrix
        root_det_dielectric = np.sqrt(np.linalg.det(self.dielectric_tensor))
        self.__ewald_param = INIT_EWALD_PARAM
        while True:
            ewald = self.__ewald_param / math.pow(self.__volume, float(1)/3) * root_det_dielectric
            max_G_vector_norm = 2 * ewald * PROD_CUTOFF_FWHM
            self.__set_G_vectors = self.__make_lattice_set(reciprocal_lattice_vect, 
                                             max_G_vector_norm, 
                                             include_self=False)
            max_R_vector_norm = PROD_CUTOFF_FWHM / ewald
            self.__set_R_vectors = self.__make_lattice_set(real_lattice, 
                                             max_R_vector_norm,
                                             include_self=True)
            print("ewald_param= %f" % self.__ewald_param)
            print("Number of R vectors = %f" % len(self.__set_R_vectors))
            print("Number of G vectors = %f" % len(self.__set_G_vectors))
            diff_real_recipro = (float(len(self.__set_R_vectors)) / len(self.__set_G_vectors)) 
            if (diff_real_recipro > 1/1.05) and (diff_real_recipro < 1.05):
                return None
            else:
                self.__ewald_param *= diff_real_recipro  ** 0.17

    def __calc_ewald_real_pot(self, atomic_pos_wrt_defect):
        """
        \sum erfc(ewald*\sqrt(R*\epsilon_inv*R)) 
                     / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
        """
        root_det_epsilon = np.sqrt(np.linalg.det(self.dielectric_tensor))
        epsilon_inv = np.linalg.inv(self.dielectric_tensor)
        ewald = self.__ewald_param / math.pow(self.__volume, float(1)/3) * root_det_epsilon
        each = np.zeros(len(self.__set_R_vectors))
        for i, R in enumerate(self.__set_R_vectors):
            # Skip the potential caused by the defect itself
            r = R - atomic_pos_wrt_defect
            if np.linalg.norm(r) < 1e-8: 
                continue
            root_R_epsilonI_R = np.sqrt(reduce(np.dot, [r.T, epsilon_inv, r]))
            each[i] = scipy.special.erfc(ewald * root_R_epsilonI_R) \
                                 / root_R_epsilonI_R

        return np.sum(each) / (4 * np.pi * root_det_epsilon)

    def __calc_ewald_recipro_pot(self, atomic_pos_wrt_defect=np.zeros(3)):
        """
        \sum exp(-G*\epsilon*G/(4*ewald**2)) / G*\epsilon*G [1/A]
        """
        root_det_epsilon = np.sqrt(np.linalg.det(self.dielectric_tensor))
        ewald = self.__ewald_param / math.pow(self.__volume, float(1)/3) * root_det_epsilon
        each = np.zeros(len(self.__set_G_vectors))

        for i, G in enumerate(self.__set_G_vectors):
            G_epsilon_G = reduce(np.dot, [G.T, self.dielectric_tensor, G]) # [1/A^2]
            each[i] = np.exp(- G_epsilon_G / 4.0 / ewald ** 2) \
                       / G_epsilon_G * np.cos(np.dot(G, atomic_pos_wrt_defect)) # [A^2]

        return np.sum(each) / self.__volume

    def __calc_ewald_self_potential(self): # [1/A]
        det_epsilon = np.linalg.det(self.dielectric_tensor)
        root_det_epsilon = np.sqrt(det_epsilon)
        ewald = self.__ewald_param / math.pow(self.__volume, float(1)/3) * root_det_epsilon
        return -ewald / (2.0 * np.pi * np.sqrt(np.pi * det_epsilon))

    def __calc_ewald_diff_potential(self):
        root_det_epsilon = np.sqrt(np.linalg.det(self.dielectric_tensor))
        ewald = self.__ewald_param / math.pow(self.__volume, float(1)/3) * root_det_epsilon
        return -0.25 / self.__volume / ewald ** 2 # [1/A]

    def __calc_model_potential(self):
        self.__model_pot = [None for i in self.__atomic_pos_wo_defect]
        coeff = \
                self.defect_property.charge * \
                scipy.constants.elementary_charge * \
                1.0e10 / scipy.constants.epsilon_0 #[V]
        for i, R in enumerate(self.__atomic_pos_wo_defect):
            real_part = self.__calc_ewald_real_pot(R-self.__defect_pos)
            recipro_part = self.__calc_ewald_recipro_pot(R)
            diff_pot = self.__calc_ewald_diff_potential()
            self.__model_pot[i] =\
                    (real_part + recipro_part + diff_pot)*coeff
            #print(real_part, recipro_part, diff_pot, coeff)
            #print(self.defect_property.charge, scipy.constants.elementary_charge, scipy.constants.epsilon_0)
        real_part = self.__calc_ewald_real_pot(np.zeros(3))
        recipro_part = self.__calc_ewald_recipro_pot(np.zeros(3))
        diff_pot = self.__calc_ewald_diff_potential()
        self_pot = self.__calc_ewald_self_potential()
        model_pot_defect_site = (real_part + recipro_part + diff_pot + self_pot) * coeff
        print(real_part, recipro_part, diff_pot, self_pot)
        self.__lattice_energy = model_pot_defect_site * self.defect_property.charge / 2

    def correct(self):
        #potential.sh step 1
        ref_pot = self.perfect_property.atomic_site_pot
        defect_pot = self.defect_property.atomic_site_pot
        #This part (removing defect_index of atomic_site_pot) must be written for defect type (Vac, Int, Sub).
        #Currently can be applied to only substitutional. 
        if len(ref_pot) != len(defect_pot):
            sys.exit("Sorry, current version can not be applied to neither vacancy nor interstitial defect.")
        ref_pot = np.delete(ref_pot, self.__defect_index-1) #reference_potential.txt
        defect_pot = np.delete(defect_pot, self.__defect_index-1) 
        diff_pot = - (defect_pot - ref_pot) # vasp_potential.txt
        print(diff_pot)
        #potential.sh step 2
        #This part (removing defect_index when calculate distances) must be written 
        #for all defect types (Vac, Int, Sub).
        #Currently can be applied to only substitutional. 
        if len(ref_pot) != len(defect_pot):
            sys.exit("Sorry, current version can not be applied to neither vacancy nor interstitial defect.")
        self.__calc_distances_from_defect() #defect_structure.txt

        #potential.sh step 0 and 3
        self.__determine_ewald_param()
        self.__calc_model_potential()
        print(self.__model_pot)
        print(self.__lattice_energy)

