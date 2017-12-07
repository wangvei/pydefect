#!/usr/bin/env python

import numpy as np
import optparse
import re
import sys

def parsePOSCAR(poscar_name):
    poscar = open(poscar_name)
    # skip the comment at first line
    poscar.readline()
    # normalization factor
    normalization = float(poscar.readline().split()[0])
    # lattice vectors
    lattice_vectors = []
    for i in range(3): 
        lattice_vectors.append([float(x) for x in poscar.readline().split()])

    lattice_vectors = np.array(lattice_vectors) * normalization
    # symbols and numbers of atoms
    line = poscar.readline().split()
    if line[0].isdigit():  # determine whether first item is integer or not.
        num_atoms = np.array([int(x) for x in line])
        symbols = ["X" for i in range(len(num_atoms))]
    else:
        symbols = line
        line = poscar.readline().split()
        num_atoms = np.array([int(x) for x in line])

    line = poscar.readline().split()
    if line[0][0] == "S":  # read one more line if selective dynamics is on.
        poscar.readline()

    # fractional coordinates of atoms
    atomic_pos = np.zeros((num_atoms.sum(), 3), dtype=float)
    for i in range(num_atoms.sum()): 
        atomic_pos[i, :] = \
                        [float(x) for x in poscar.readline().split()[0:3]]
    poscar.close()
    return lattice_vectors, symbols, num_atoms, atomic_pos

def printPOSCAR(lattice_vectors,
                num_atoms, 
                atomic_pos,
                header="", 
                normalization=1.0,
                symbols=None):
    print header
    print normalization
    for i in range(3):
        print "%22.16f %22.16f %22.16f" % tuple(lattice_vectors[i])
    if symbols: 
        print "%3s " * len(symbols) % tuple(symbols)
    print "%3i " * len(num_atoms) % tuple(num_atoms)
    print "Direct"
    for i in range(num_atoms.sum()):
        print "%22.16f %22.16f %22.16f" % tuple([p - np.floor(p) for p in atomic_pos[i]])
#        print "%22.16f %22.16f %22.16f" % tuple(atomic_pos[i])

def get_distance(lattice_vectors, fractional_vectors):
    return \
        np.linalg.norm(np.dot(lattice_vectors.transpose(), fractional_vectors))

def get_cartesian_coordinates(lattice_vectors, fractional_vectors):
    return np.array(
           [np.dot(lattice_vectors.transpose(), f) for f in fractional_vectors])

def get_reciprocal_lattice_wo_2Pi(lattice_vectors):  # wo 2*Pi
    reciprocal_lattice = np.array( 
            [np.cross(lattice_vectors[i-2], lattice_vectors[i-1]) 
                        / get_volume(lattice_vectors) for i in range(3)])
    return reciprocal_lattice

def get_reciprocal_lattice(lattice_vectors):
    return 2.0 * np.pi * get_reciprocal_lattice_wo_2Pi(lattice_vectors)

def get_half_maximum_lattice(lattice_vectors):
    return np.max(get_lattice_constants(lattice_vectors)) / 2.0

def get_volume(lattice_vectors):
    return np.linalg.det(lattice_vectors)

def get_distance_two_planes(lattice_vectors):
    # (a_i \times a_j) \ddot a_k / |a_i \times  a_j| 
    distance = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_times_a_j = np.cross(lattice_vectors[i-2], lattice_vectors[i-1])
        a_k = lattice_vectors[i]
        distance[i] = abs(np.dot(a_i_times_a_j, a_k)) \
                         / np.linalg.norm(a_i_times_a_j)
    return distance

def get_sphere_radius(lattice_vectors):
    # Maximum radius of a sphere fitting inside the unit cell.
    return min(get_distance_two_planes(lattice_vectors)) / 2.0

def get_lattice_constants(lattice_vectors):
    return [np.linalg.norm(lattice_vectors[i]) for i in range(3)]

def get_angles(lattice_vectors):
    angles = np.zeros(3, dtype=float)

    for r in range(3):
        angles[r] = np.dot(lattice_vectors[r-2], lattice_vectors[r-1]) \
                  / np.linalg.norm(lattice_vectors[r-2]) \
                  / np.linalg.norm(lattice_vectors[r-1])
        angles[r] = np.arccos(angles[r]) * 180.0/np.pi
    return angles

parser = optparse.OptionParser()
parser.add_option("-p", "--poscar", 
                  dest="poscar", 
                  type="string", 
                  default="POSCAR", 
                  help="POSCAR file for the atomic positions.", 
                  metavar="FILE")
parser.add_option("-r", "--radius", 
                  dest="radius", 
                  action="store_true",
                  default=False,
                  help="Print a radius of a sphere fitting inside the cell.")

if __name__ == "__main__":
    (opts, args) = parser.parse_args()
    (lattice_vectors, symbols, num_atoms, atomic_pos) \
                                    = parsePOSCAR(opts.poscar)
    if opts.radius:
        print "Maximum radius of spheres on lattice point occupying space = %10.5f [A]"\
                           % (get_half_maximum_lattice(lattice_vectors))
        sys.exit(0)
    lattice_constants = get_lattice_constants(lattice_vectors)
    lattice_angles = get_angles(lattice_vectors)
    print get_reciprocal_lattice(lattice_vectors) 
    print "lattice matrix [A]"
    print "%12.7f %12.7f %12.7f" % tuple(lattice_vectors[0])
    print "%12.7f %12.7f %12.7f" % tuple(lattice_vectors[1])
    print "%12.7f %12.7f %12.7f" % tuple(lattice_vectors[2])
    print "lattice constatns [A]"
    print "%12.7f %12.7f %12.7f" % tuple(lattice_constants)
    print "lattice angles [deg.]"
    print "%12.4f %12.4f %12.4f" % tuple(lattice_angles)
    print "lattice volume [A^3]"
    print "%12.4f" % get_volume(lattice_vectors)
#    print interatomic_distances(lattice_vectors,atomic_pos,[0,0,0])
