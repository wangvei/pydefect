#!/usr/bin/env python

import numpy as np
import optparse 
import parse_poscar as ppos
import sys  

def defect_position(defect, num_atoms, atom_pos):
    defect=defect.split()

    if len(defect) == 3:
        defect_pos = np.array([float(defect[x]) for x in range(3)])
        # If the defect is vacancy remove is set to False, otherwise True.
        remove = False
    elif len(defect) == 1:
        defect_pos = atom_pos[int(defect[0]) - 1]
        atom_pos = np.delete(atom_pos, int(defect[0]) - 1, 0)
        for i in range(len(num_atoms)):
            if int(defect[0]) <= num_atoms[:i+1].sum():
                num_atoms[i] -= 1
                break
        remove = True

    return defect_pos, num_atoms, atom_pos, remove

def defect_centered_POSCAR(poscar,defect):
    (lattice_vectors, symbols, num_atoms, atom_pos) = \
                                      ppos.parsePOSCAR(poscar)
    # Set the defect position at the center.
    defect_pos, num_atoms, atom_pos, remove = \
                          defect_position(defect, num_atoms, atom_pos)

    atom_pos_wrt_defect = atom_pos - defect_pos + 0.5
    atom_pos_wrt_defect = np.insert(atom_pos_wrt_defect, 
                                      0, [0.5, 0.5, 0.5], axis=0)
    if remove:
        comment="The " + defect + "th atom is centerd."
    else:
        comment="The vacancy at " + repr(defect_pos) + " is centered."

    num_atoms = np.insert(num_atoms,0,1)
    symbols.insert(0, "Def")
    ppos.printPOSCAR(lattice_vectors, num_atoms, 
                     atom_pos_wrt_defect, comment, symbols, 1.0)

def shortest_distance(lattice_vectors, fractional_coord):
    distance = ppos.get_distance(lattice_vectors, fractional_coord)

    # index = [[-1,-1,-1],[-1,-1,0], ...., [1,1,1]]
    index = \
         [[i-1, j-1, k-1] for i in range(3) for j in range(3) for k in range(3)]

    while True:
        candidate = {}                                                       
        for i in index:                                                      
            # [[x1, y1, z1], [x2, y2, z2], ... ,[x9, y9, z9]]
            candidate_position = [fractional_coord[y] - i[y] for y in range(3)]
            # candidate[x1, y1, z1] = distance1
            candidate[tuple(candidate_position)] = \
                    ppos.get_distance(lattice_vectors, candidate_position)   

        if min(candidate.values()) + sys.float_info.epsilon < distance:         
            distance = min(candidate.values())         
            # get the key to the entry which contains the minimum value.     
            fractional_coord = min(candidate, key=candidate.get)        
        else:                                                                
            break 

    return distance

def distances_from_point(lattice_vectors, atom_pos, point):
    """
    Determine the shortest distance at each atom from the image points.
    The shortest one between atom and the defect is seeked for iteratively.
    This algorithm may not be perfect, but practically enough.
    """
    atom_pos_wrt_defect = atom_pos - point

    return np.array([shortest_distance(lattice_vectors, c) 
                                                  for c in atom_pos_wrt_defect])
         
def print_distances_from_defect(poscar, defect):
    lattice_vectors, symbols, num_atoms, atom_pos = ppos.parsePOSCAR(poscar)

    defect_pos, num_atoms, atom_pos, remove = \
                                    defect_position(defect, num_atoms, atom_pos)

    atomic_distance_from_defect = \
                     distances_from_point(lattice_vectors, atom_pos, defect_pos)

    # Print the list of neighboring atomic distances.
    print "-- Defect position: %12.7f %12.7f %12.7f" % tuple(defect_pos)
    print "----------------------------------------------------------"
    print "symbol number  ------- coordinations -------  distances[A]"
    
    x = 0
    for i, n in enumerate(num_atoms):
        for j in range(n):
            print_tuple = \
      (symbols[i], x+1) + tuple(atom_pos[x]) + (atomic_distance_from_defect[x],)
            print "%6s %5i %10.6f %10.6f %10.6f %12.5f" % print_tuple
            x += 1

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-p", "--poscar",
                      dest="poscar",
                      help="Set POSCAR-type file name.",
                      metavar="FILE")
    parser.add_option("-d", "--defect",
                      dest="defect",
                      help="Input the position of the defect in fractional \
                            coordinate or atomic number. \
                            E.g. '0.5 0.5 0.5' or 82")
    parser.add_option("-c", "--centering",
                      dest="centering",
                      action="store_true",                                      
                      default=False,
                      help="Center the defect and print POSCAR type format.")
    parser.add_option("--cutoff",
                      dest="cutoff",
                      type="float",
                      help="Cutoff radius for plotting defect neighbors.")
    (opts, args) = parser.parse_args()

    if opts.centering: 
        defect_centered_POSCAR(opts.poscar,opts.defect)
    else:
        print_distances_from_defect(opts.poscar,opts.defect)  
