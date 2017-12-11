#!/usr/bin/env python

import sys
import math
import numpy as np
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
        ewald_param = INIT_EWALD_PARAM
        volume = np.linalg.det(real_lattice)
        while True:
            ewald = ewald_param / math.pow(volume, float(1)/3) * root_det_dielectric
            max_G_vector_norm = 2 * ewald * PROD_CUTOFF_FWHM
            set_G_vectors = self.__make_lattice_set(reciprocal_lattice_vect, 
                                             max_G_vector_norm, 
                                             include_self=False)
            max_R_vector_norm = PROD_CUTOFF_FWHM / ewald
            set_R_vectors = self.__make_lattice_set(real_lattice, 
                                             max_R_vector_norm,
                                             include_self=True)
            print("ewald_param= %f" % ewald_param)
            print("Number of R vectors = %f" % len(set_R_vectors))
            print("Number of G vectors = %f" % len(set_G_vectors))
            diff_real_recipro = (float(len(set_R_vectors)) / len(set_G_vectors)) 
            if (diff_real_recipro > 1/1.05) and (diff_real_recipro < 1.05):
                return ewald_param
            else:
                ewald_param *= diff_real_recipro  ** 0.17

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
        ewald_param = self.__determine_ewald_param()
        print(ewald_param)


