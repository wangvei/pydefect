import numpy as np
import optparse 
import parse_poscar as ppos
import sys  

def defect_position(defect, num_atoms, atomic_pos):
    defect=defect.split()
    if len(defect) == 3:
        defect_pos = np.array([float(defect[x]) for x in range(3)])
        # If the defect is vacancy, remove is set to False, otherwise True.
        remove = False
    elif len(defect) == 1:
        defect_pos = atomic_pos[int(defect[0]) - 1]
        atomic_pos = np.delete(atomic_pos, int(defect[0]) - 1, 0)
        for i in range(len(num_atoms)):
            if int(defect[0]) <= num_atoms[:i+1].sum():
                num_atoms[i] -= 1
                break
        remove = True

    return defect_pos, num_atoms, atomic_pos, remove

def defect_centered_POSCAR(poscar,defect):
    (lattice_vectors, symbols, num_atoms, atomic_pos) = \
                                      ppos.parsePOSCAR(poscar)
    # Set the defect position at the center.
    defect_pos, num_atoms, atomic_pos, remove = \
                          defect_position(defect, num_atoms, atomic_pos)

    atomic_pos_wrt_defect = atomic_pos - defect_pos + 0.5
    atomic_pos_wrt_defect = np.insert(atomic_pos_wrt_defect, 
                                      0, [0.5, 0.5, 0.5], axis=0)
    if remove:
        comment="The " + defect + "th atom is centerd."
    else:
        comment="The vacancy at " + repr(defect_pos) + " is centered."

    num_atoms = np.insert(num_atoms,0,1)
    symbols.insert(0, "Def")
    ppos.printPOSCAR(lattice_vectors, num_atoms, 
                     atomic_pos_wrt_defect, comment, 1.0, symbols)

def distances_from_defect(poscar,defect):
    ############################################################################
    # Determine the shortest distance at each atom from the image defects.
    # The shortest distance between atom and the defect is searched for iteratively.
    # This algorithm may not be perfect.
    ############################################################################
    (lattice_vectors, symbols, num_atoms, atomic_pos) \
                                    = ppos.parsePOSCAR(poscar)
    # Set the defect position
    defect_pos, num_atoms, atomic_pos, remove = \
                            defect_position(defect, num_atoms, atomic_pos)

    atomic_pos_wrt_defect = atomic_pos - defect_pos
    atomic_distance_from_defect = \
              [ppos.get_distance(lattice_vectors, a) for a in atomic_pos_wrt_defect]

    # index = [[-1,-1,-1],[-1,-1,0], ...., [1,1,1]]
    index = \
         [[i-1, j-1, k-1] for i in range(3) for j in range(3) for k in range(3)]

    for x, a in enumerate(atomic_pos_wrt_defect):
        while True:
            candidate = {}
            for i in index:
                candidate_position = [a[y] - i[y] for y in range(3)]
                candidate[tuple(candidate_position)] = \
                        ppos.get_distance(lattice_vectors, candidate_position)
    
            if min(candidate.values()) < atomic_distance_from_defect[x]:
                atomic_distance_from_defect[x] = min(candidate.values())
                # get the key to the entry which contains the minimum value.
                atomic_pos_wrt_defect = min(candidate, key=candidate.get) 
            else:
                break
    
    # Print the list of neighboring atomic distances.
    print "-- Defect position: %12.7f %12.7f %12.7f" % tuple(defect_pos)
    print "----------------------------------------------------------"
    print "symbol number  ------- coordinations -------  distances[A]"
    
    x = 0
    for i, n in enumerate(num_atoms):
        for j in range(n):
            print_tuple = \
    (symbols[i], x+1) + tuple(atomic_pos[x]) + (atomic_distance_from_defect[x],)
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
        distances_from_defect(opts.poscar,opts.defect)  
