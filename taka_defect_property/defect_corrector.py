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

    @property
    def defect_property(self):
        return self.__defect_property

    @property
    def perfect_property(self):
        return self.__perfect_property


    def __determine_ewald_param(self, real_lattice, reciprocal_lattice):
        INIT_EWALD_PARAM = 0.010111311355097064
        PROD_CUTOFF_FWHM = 25.0 #product of cutoff radius of G-vector and gaussian FWHM.
                                #increasing this value, both accuracy and computational cost will be increased
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
                if index == (0, 0, 0): 
                    if not include_self:
                        continue
                vector = np.dot(lattice_vectors.transpose(), np.array(index))
                norm = np.linalg.norm(vector)
                if norm < max_length: 
                    vectors.append(vector)
            return np.array(vectors)

        root_det_dielectric = np.sqrt(np.linalg.det(self.dielectric_tensor))
        ewald_param = INIT_EWALD_PARAM
        volume = real_lattice.volume
        real_lattice_vect = real_lattice.matrix
        reciprocal_lattice_vect = reciprocal_lattice.matrix
        while True:
            ewald = ewald_param / math.pow(volume, float(1)/3) * root_det_dielectric
            max_G_vector_norm = 2 * ewald * PROD_CUTOFF_FWHM
            set_G_vectors = make_lattice_set(reciprocal_lattice_vect, 
                                             max_G_vector_norm, 
                                             include_self=False)
            max_R_vector_norm = PROD_CUTOFF_FWHM / ewald
            set_R_vectors = make_lattice_set(real_lattice_vect, 
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
        #potential.sh 0 0
        real_lattice = self.perfect_property.structure.lattice
        reciprocal_lattice = self.perfect_property.structure.lattice.reciprocal_lattice
        ewald_param = self.__determine_ewald_param(real_lattice, reciprocal_lattice)
        print(ewald_param)
        #potential.sh step 1
        ref_pot = self.perfect_property.atomic_site_pot
        defect_pot = self.defect_property.atomic_site_pot
        print(ref_pot)
        print(defect_pot)
        defect_pot[0]="NISHIYA!!!!!!!!!!!!!!!!!!!!!!!!!"
        print(self.defect_property.atomic_site_pot)
        print(defect_pot)
        

    

