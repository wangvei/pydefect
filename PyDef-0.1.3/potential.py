#!/usr/bin/env python

"""
This script allows you to calculate the potential caused by the point charges 
under periodic boundary condition. Dielectric responce is considered with 
dielectric tensor.

Ref: Y. Kumagai and F. Oba, Phys. Rev. B 89, 195205 (2014).
When using this code in order to correct the defect formation energies,
cite it, please.

Input:
  1. Name of the "POSCAR"-format file in vasp. (-p or --poscar)
  2. Charge state q (-c or --charge).
  3. The position of point charge in charge state q in fractional coordinate
     (-d or --defect). You can choose whether to remove defect site or not 
     for estimating the potential.
  4. Dielectric tensor (-e or --epsilon). When allowing the atomic 
     relaxation, this should be the sum of electronic and ionic parts.
  5. Optional accuracy parameter (--accuracy).
"""

import math
import numpy as np
import parse_poscar as ppos
import scipy.special
import sys
from defect_structure import defect_position
from optparse import OptionParser
from physical_constants import ElementaryCharge, VacuumPermittivity

def ewald_real_part_potential(ewald, set_R, epsilon):
    """
    -- input --
    ewald   : Ewalt parameter [1/A]
    set_R   : A set of lattice vector [[Rx1,Ry1,Rz1],...] [A]
    epsilon : Dielectric tensor

    -- output --
    \sum erfc(ewald*\sqrt(R*\epsilon_inv*R)) 
                 / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
    """
    root_det_epsilon = np.sqrt(np.linalg.det(epsilon))
    epsilon_inv = np.linalg.inv(epsilon)
    each = np.zeros(len(set_R))
 
    for i, R in enumerate(set_R):
        # Skip the potential caused by the defect itself
        if np.linalg.norm(R) < 1e-8: 
            continue
        root_R_epsilonI_R = np.sqrt(reduce(np.dot, [R.T, epsilon_inv, R]))
        each[i] = scipy.special.erfc(ewald * root_R_epsilonI_R) \
                             / root_R_epsilonI_R

    return np.sum(each) / (4 * np.pi * root_det_epsilon)

def ewald_reciprocal_part_potential(
                            ewald, set_G, epsilon, volume, position=[0,0,0]):
    """
    -- input --
    ewald   : Ewalt parameter [1/A]
    set_G   : A set of reciprocal lattice vector [[Gx1,Gy1,Gz1],...] [1/A]
    epsilon : dielectric tensor

    -- output --
    \sum exp(-G*\epsilon*G/(4*ewald**2)) / G*\epsilon*G [1/A]
    """
    each = np.zeros(len(set_G))

    for i, G in enumerate(set_G):
        G_epsilon_G = reduce(np.dot, [G.T, epsilon, G]) # [1/A^2]
        each[i] = np.exp(- G_epsilon_G / 4.0 / ewald ** 2) \
                   / G_epsilon_G * np.cos(np.dot(G, position)) # [A^2]

    return np.sum(each) / volume

def ewald_self_potential(ewald, epsilon): # [1/A]
    return -ewald / (2.0 * np.pi * np.sqrt(np.pi * np.linalg.det(epsilon)))

def ewald_diff_potential(ewald, volume):
    return -0.25 / volume / ewald ** 2 # [1/A]

def lattice_set(lattice_vectors, max_length, self): 
    """
    Return a set of lattice vectors within the max length.
    Note that angles between any two axes are assumed to be between 60 and 
    120 deg.
    """
    max_integer = \
    [int(max_length / np.linalg.norm(lattice_vectors[i])) + 1 for i in range(3)]
    vectors = []

    for i in range(-max_integer[0], max_integer[0] + 1):
        for j in range(-max_integer[1], max_integer[1] + 1):
            for k in range(-max_integer[2], max_integer[2] + 1):
                index = [i, j, k]
                if index == [0, 0, 0]: 
                    if self == False: continue
                vector = np.dot(lattice_vectors.transpose(), index)
                norm = np.linalg.norm(vector)
                if norm < max_length: 
                    vectors.append(vector)

    return np.array(vectors)

parser = OptionParser()
parser.add_option("-p", "--poscar", 
                  dest="poscar", 
                  type="string", 
                  default="POSCAR", 
                  help="POSCAR-type file for the atomic positions.", 
                  metavar="FILE")
parser.add_option("-c", "--charge", 
                  dest="charge", 
                  type="float", 
                  help="Point charge at the defect position [e].")
parser.add_option("-d", "--defect",
                  dest="defect",
                  help="Input the position of the defect in fractional \
                        coordinate or atomic number. \
                        E.g. '0.5 0.5 0.5' or 82")
parser.add_option("-e", "--epsilon",  
                  dest="epsilon", 
                  help="Dielectric constant tensor.\
                        E.g.1 '10' :only a single component (cubic) \
                        E.g.2 '10 12 14' :only diagonal components \
                        E.g.3 '42 0 -1 0 42 0 -1 0 16' : in general")
parser.add_option("--omega",
                  dest="omega", 
                  default=0.01, 
                  type="float", 
                  help="Determine Ewalt parameter.")
parser.add_option("--accuracy",
                  dest="accuracy", 
                  default=25, 
                  type="float", 
                  help="Cutoff radius of G-vector * gaussian FWHM.")
parser.add_option("--check",
                  dest="check", 
                  action="store_true",
                  help="Check the number of G vectors and R vectors.")

(opts, args) = parser.parse_args()

charge = opts.charge
lattice_vectors, symbols, num_atoms, atomic_pos = ppos.parsePOSCAR(opts.poscar)
volume = ppos.get_volume(lattice_vectors)

try:
    epsilon_each =opts.epsilon.split()
    if len(epsilon_each) == 9:
        epsilon = \
            np.array([float(epsilon_each[x]) for x in range(9)]).reshape(3,3)
    elif len(epsilon_each) == 3:
        epsilon_each = [float(epsilon_each[x]) for x in range(3)]
        epsilon = np.array([[epsilon_each[0], 0.0, 0.0],
                            [0.0, epsilon_each[1], 0.0],
                            [0.0, 0.0, epsilon_each[2]]])
    elif len(epsilon_each) == 1:
        epsilon_each = float(epsilon_each[0])
        epsilon = np.array([[epsilon_each, 0.0, 0.0],
                            [0.0, epsilon_each, 0.0],
                            [0.0, 0.0, epsilon_each]])

except:
    print "The number of the dielectric components are not 1, 3, or 9"
    sys.exit(1)

# Get the reciprocal lattice vectors with 2pi. [[norm, G], [norm, G], ...]
reciprocal_lattice_vectors = ppos.get_reciprocal_lattice(lattice_vectors)

# Get a set of R and G vectors for summing the real space and reciprocal space.
root_det_epsilon = np.sqrt(np.linalg.det(epsilon))

# Roughly, omega = 0.05 * average of lattice constant
if opts.omega:
    omega = opts.omega
else:
    omega = 0.05 * np.average([np.linalg.norm([l]) for l in lattice_vectors])

accuracy = opts.accuracy
# In the paper, gamma is used instead of ewald to indicate Ewald parameter.
ewald = omega / math.pow(volume, float(1)/3) * root_det_epsilon

max_G_vector_norm = 2 * ewald * accuracy
set_G_vectors = \
         lattice_set(reciprocal_lattice_vectors, max_G_vector_norm, self=False)
max_R_vector_norm = accuracy / ewald
set_R_vectors = lattice_set(lattice_vectors, max_R_vector_norm, self=True)

if opts.check:
    print "omega=", omega
    print "Number of R vectors =", len(set_R_vectors)
    print "Number of G vectors =", len(set_G_vectors)
    while float(len(set_R_vectors)) / len(set_G_vectors) > 1.05 or \
          float(len(set_G_vectors)) / len(set_R_vectors) > 1.05:
        omega *= (float(len(set_R_vectors)) / len(set_G_vectors))  ** 0.17
        print "omega=", omega
        ewald = omega / math.pow(volume, float(1)/3) * root_det_epsilon
        max_G_vector_norm = 2 * ewald * accuracy
        set_G_vectors = \
         lattice_set(reciprocal_lattice_vectors, max_G_vector_norm, self=False)
        max_R_vector_norm = accuracy / ewald
        set_R_vectors = \
         lattice_set(lattice_vectors, max_R_vector_norm, self=True)
        print "Number of R vectors =", len(set_R_vectors)
        print "Number of G vectors =", len(set_G_vectors)

    sys.exit(0)

defect_pos, num_atoms, atomic_pos, remove = \
                      defect_position(opts.defect, num_atoms, atomic_pos)

"""
Units are as follows:
  G [1/A]
  R [A]
  r [A]
  ewald [1/A]
  charge [|e|]
  ElementaryCharge [C/|e|]
  VacuumPermittivity [F/m] = [C/V/m]
  volume [A^3]    

Unit conversion
1 [|e|/A] = ElementaryCharge * 1.0e10 / VacuumPermittivity [V]

** CAUTION ** 
In SI unit, eqs.(8) and (14) in PRB 89, 195205 (2014) should be divided by 4Pi.
"""

atomic_pos_wrt_defect = \
     ppos.get_cartesian_coordinates(lattice_vectors, atomic_pos - defect_pos)

for i in range(num_atoms.sum()):
    set_r_vectors = set_R_vectors - atomic_pos_wrt_defect[i]
    phi = (ewald_real_part_potential(ewald, set_r_vectors, epsilon)\
        +  ewald_reciprocal_part_potential(ewald, set_G_vectors, 
                                     epsilon, volume, atomic_pos_wrt_defect[i])\
        +  ewald_diff_potential(ewald, volume)) \
        *  charge * ElementaryCharge *  1.0e10 / VacuumPermittivity # [V]

    print phi

phi_defect_site = \
      (ewald_real_part_potential(ewald, set_R_vectors, epsilon) \
     + ewald_reciprocal_part_potential(ewald, set_G_vectors, epsilon, volume) \
     + ewald_diff_potential(ewald, volume) \
     + ewald_self_potential(ewald,epsilon)) \
     * charge * ElementaryCharge * 1.0e10 / VacuumPermittivity # [V]

# Lattice energy must be divided by 2 for cancelling out double-count.
lattice_energy = phi_defect_site * charge / 2 # [eV] 

print " E_lattice: %10.7f [eV]" % (lattice_energy)
sys.exit(0)
