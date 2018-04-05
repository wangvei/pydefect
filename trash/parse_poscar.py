#!/usr/bin/env python

import argparse
import defect_structure
import mathematics as lmath
import numpy as np
import re
import spglib as spg
import sys


# def cell_for_spglib(poscar_name, cell):
#
def make_spg_cell(poscar_name):
    lattice_vectors, symbols, num_atoms, atom_pos = parsePOSCAR(poscar_name)
    distinguish_atoms = \
        [i for i in range(len(num_atoms)) for j in range(num_atoms[i])]

    cell = tuple([lattice_vectors, atom_pos, distinguish_atoms])

    return cell


def symmetryPOSCAR(poscar_name, use_magmoms=False,
                   symprec=1e-5, angle_tolerance=-1.0):
    cell = make_spg_cell(poscar_name)

    return spg.get_symmetry_dataset(cell)


def get_atom_mapping(poscar_name, use_magmoms=False,
                     symprec=1e-5, angle_tolerance=-1.0):
    spglib_atom_mapping = \
        symmetryPOSCAR(poscar_name, use_magmoms=use_magmoms, symprec=symprec,
                       angle_tolerance=angle_tolerance)["equivalent_atoms"]

    atom_map = {}

    for i, a in enumerate(spglib_atom_mapping):

        atom_index = a + 1

        if not atom_index in atom_map:
            atom_map[atom_index] = []

        atom_map[atom_index].append(i + 1)

    return atom_map


def parsePOSCAR(poscar_name, potcar_name=None):
    poscar = open(poscar_name)
    # skip the comment at the first line
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
    if line[0][0] == "S":  # Read 1 line if the selective dynamics is turned on.
        poscar.readline()

    # Fractional coordinates of atoms
    atom_pos = np.zeros((num_atoms.sum(), 3), dtype=float)
    for i in range(num_atoms.sum()):
        atom_pos[i, :] = \
            [float(x) for x in poscar.readline().split()[0:3]]
    poscar.close()
    # lattice_vectors, num_atoms and atom_pos are numpy arrays.
    return lattice_vectors, symbols, num_atoms, atom_pos


def printPOSCAR(lattice_vectors, num_atoms, atom_pos, header="",
                symbols=None, normalization=1.0):
    print header
    print normalization

    for i in range(3):
        print "%22.16f %22.16f %22.16f" % tuple(lattice_vectors[i])

    if symbols:
        print "%3s " * len(symbols) % tuple(symbols)
    print "%3i " * len(num_atoms) % tuple(num_atoms)
    print "Direct"

    for i in range(num_atoms.sum()):
        print "%22.16f %22.16f %22.16f" % \
              tuple(p - np.floor(p) for p in atom_pos[i])


def writePOSCAR(lattice_vectors, num_atoms, atom_pos, file_name, header="",
                symbols=None, normalization=1.0, note=None):
    f = open(file_name, 'w')
    print >> f, header
    print >> f, normalization

    for i in range(3):
        print >> f, "%22.16f %22.16f %22.16f" % tuple(lattice_vectors[i])

    if symbols:
        print >> f, "%3s " * len(symbols) % tuple(symbols)

    print >> f, "%3i " * len(num_atoms) % tuple(num_atoms)
    print >> f, "Direct"

    x = 0
    for i in range(len(num_atoms)):
        for j in range(num_atoms[i]):
            if note is None:
                print >> f, "%22.16f %22.16f %22.16f %3s" % \
                            tuple(p - np.floor(p) for p in atom_pos[x]) + symbols[i]
            else:
                print >> f, "%22.16f %22.16f %22.16f %3s %2s" % \
                            tuple([p - np.floor(p) for p in atom_pos[x]] + [symbols[i], note[x]])

            x += 1

    #        print >> f, "%22.16f %22.16f %22.16f" % \
    #                            tuple(p - np.floor(p) for p in atom_pos[i])

    f.close()


def write_empty_sphere(site, file_name):
    f = open(file_name, 'a')
    print >> f, "Empty (spheres)"
    print >> f, "1"
    print >> f, "%22.16f %22.16f %22.16f" % tuple(site)
    f.close()


def make_point_defects(poscar_name, removed_atom, added_coord, added_atom_symbol,
                       output_poscar_name, random_parameters=False, header="", empty=True):
    # removed_atom : atom index
    # added_coord : fractional coordinates:
    # added_atom_symbol : element name
    # random_parameters["cutoff"]: Cutoff [A] to determine atoms for considering the randomization.
    # random_parameters["max_disp"]: Maximum distance for the random displacement.
    # -------------------------------------------------------------------------
    # removed_atom yes & added_coord yes & added_atom_symbol yes -> vac + int
    # removed_atom  no & added_coord yes & added_atom_symbol yes -> interstitial
    # removed_atom yes & added_coord  no & added_atom_symbol yes -> substitute
    # removed_atom yes & added_coord  no & added_atom_symbol  no -> vacancy

    lattice_vectors, symbols, num_atoms, atom_pos = parsePOSCAR(poscar_name)

    # Removed one atom from the original atom set in poscar_name.
    if removed_atom:
        for i in range(len(num_atoms)):
            if sum(num_atoms[0:i + 1]) >= removed_atom:
                num_atoms[i] -= 1
                break

        # Removed or replaced site. This is also used for empty sphere.
        vacant_coord = atom_pos[removed_atom - 1]
        atom_pos = np.delete(atom_pos, removed_atom - 1, axis=0)

    # ------ define defect site ------
    # defect_site: vacancy + intersitial
    if removed_atom and added_coord:
        substituted_coord = added_coord
        defect_coord = [(vacant_coord[i] + added_coord[i]) / 2 for i in range(3)]
    # defect_site: interstitial
    elif added_coord:
        substituted_coord = added_coord
        defect_coord = added_coord
    # defect_site: substitutional
    elif added_atom_symbol:
        substituted_coord = vacant_coord
        defect_coord = vacant_coord
    else:
        # defect_site: vacancy
        defect_coord = vacant_coord
    # --------------------------------

    if added_atom_symbol:
        # native defects
        if added_atom_symbol in symbols:
            num_atoms[symbols.index(added_atom_symbol)] += 1
            add_index = np.sum(num_atoms[0:symbols.index(added_atom_symbol)])
        # foreign species
        else:
            symbols = [added_atom_symbol] + symbols
            num_atoms = np.insert(num_atoms, 0, 1, axis=0)
            add_index = 0

        atom_pos = np.insert(atom_pos, add_index, substituted_coord, axis=0)

    if random_parameters is not False:
        atom_pos, displaced_atom_index = random_displaces(lattice_vectors,
                 atom_pos, defect_coord, random_parameters["cutoff"], random_parameters["max_disp"])
    else:
        displaced_atom_index = ["" for i in range(len(atom_pos))]

        # "interstitial" or "substitutional"
    if (not removed_atom and added_coord and added_atom_symbol) \
            or (removed_atom and not added_coord and added_atom_symbol):
        header += " DefectSupercell index " + str(add_index + 1)
    # the others
    else:
        header += " DefectSupercell position " + " ".join(map(str, defect_coord))

    writePOSCAR(lattice_vectors, num_atoms, atom_pos, output_poscar_name,
                header, symbols=symbols, note=displaced_atom_index)

    # Add an empty sphere in cases of "vacancy + intersitial" or "vacancy".


#    if (removed_atom and not added_coord and not added_atom_symbol and empty) \
#    or (removed_atom and     added_coord and     added_atom_symbol and empty):
#        write_empty_sphere(vacant_coord, output_poscar_name)

def random_displace(lattice_vectors, distance):
    random_disp_cartesian_vector = \
        lmath.random_three_vector() * distance * np.random.random()

    return np.dot(np.linalg.inv(lattice_vectors), random_disp_cartesian_vector)


def random_displaces(lattice_vectors, atom_pos, center, cutoff, distance):
    """
    Add random displacement of [0, *distance*) lenght if the distance from
    *center* to *atom_pos* is shorter than cutoff.
    """
    distances = defect_structure.distances_from_point(
        lattice_vectors, atom_pos, center)
    displaced_atom = []

    for i, d in enumerate(distances):
        if d < cutoff:
            atom_pos[i] += random_displace(lattice_vectors, distance)
            displaced_atom.append("D")
        else:
            displaced_atom.append("")

    # Return the set of fractional coordinates and signature of displaced atom.
    return atom_pos, displaced_atom


def get_distance(lattice_vectors, fractional_vectors):
    return np.linalg.norm(np.dot(lattice_vectors, fractional_vectors))


def get_cartesian_coordinates(lattice_vectors, fractional_vectors):
    return np.array([np.dot(lattice_vectors, f) for f in fractional_vectors])


def get_reciprocal_lattice_wo_2Pi(lattice_vectors):  # wo 2*Pi
    reciprocal_lattice = np.array(
        [np.cross(lattice_vectors[i - 2], lattice_vectors[i - 1])
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
        a_i_times_a_j = np.cross(lattice_vectors[i - 2], lattice_vectors[i - 1])
        a_k = lattice_vectors[i]
        distance[i] = abs(np.dot(a_i_times_a_j, a_k)) \
                      / np.linalg.norm(a_i_times_a_j)
    return distance


def get_min_sphere_radius(lattice_vectors):
    # Minimum radius of a sphere fitting inside the unit cell.
    return min(get_distance_two_planes(lattice_vectors)) / 2.0


def get_max_sphere_radius(lattice_vectors):
    # Maximum radius of a sphere fitting inside the unit cell.
    return max(get_distance_two_planes(lattice_vectors)) / 2.0


def get_lattice_constants(lattice_vectors):
    return [np.linalg.norm(lattice_vectors[i]) for i in range(3)]


def get_angles(lattice_vectors):
    angles = np.zeros(3, dtype=float)

    for r in range(3):
        angles[r] = np.dot(lattice_vectors[r - 2], lattice_vectors[r - 1]) \
                    / np.linalg.norm(lattice_vectors[r - 2]) \
                    / np.linalg.norm(lattice_vectors[r - 1])
        angles[r] = np.arccos(angles[r]) * 180.0 / np.pi

    return angles


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                    type=str, help="POSCAR name.")
parser.add_argument("--nion", dest="nion", action="store_true",
                    help="Print a number of ions in POSCAR file.")
parser.add_argument("-r", "--radius", dest="radius", action="store_true",
                    help="Print a radius of a sphere fitting inside the cell.")
parser.add_argument("-v", "--volume", dest="volume", action="store_true",
                    help="Print volume of the cell.")

if __name__ == "__main__":
    opts = parser.parse_args()
    (lattice_vectors, symbols, num_atoms, atom_pos) \
        = parsePOSCAR(opts.poscar)
    if opts.nion:
        print sum(num_atoms)
    elif opts.radius:
        print "Maximum radius of spheres on lattice point occupying space = %10.5f [A]" \
              % (get_max_sphere_radius(lattice_vectors))
    elif opts.volume:
        print "%12.4f" % get_volume(lattice_vectors)
    else:
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
        #    print symmetryPOSCAR(opts.poscar)
        print "Atom mapping"
        print get_atom_mapping(opts.poscar)
        #    print interatomic_distances(lattice_vectors,atom_pos,[0,0,0])
        sys.exit(0)
