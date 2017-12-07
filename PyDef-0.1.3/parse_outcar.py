import numpy as np
import optparse 
import sys  

"""

"""
def get_val(outcar_name, tag):
    outcar = open(outcar_name,"r")

    for line in outcar:
        if line.find(tag) != -1: 
            if   tag == "NIONS":  
                val = int(line.split()[11])
                break
            elif tag == "NKPTS":  
                val = int(line.split()[3])
                break
            elif tag == "NBANDS":
                val = int(line.split()[14])
                break
            elif tag == "ISPIN":
                val = int(line.split()[2])
                break
            elif tag == "LNONCOLLINEAR":
                val = line.split()[2]
                break
    outcar.close()

    return val

#def symmetry_equivalent_atom(outcar_name, num_atoms):
#    outcar = open(outcar_name,"r")
#    while True:
#        line = outcar.readline()    
#        if not line: break
#        if line.find('irot  :') != -1:
            

def potential(outcar_name, num_atoms):
    outcar = open(outcar_name,"r")
    num_potential_lines = num_atoms / 5 + 1

    # Since outcar file is read in the field, "while True:" is essential.
    while True:
        line = outcar.readline()
        if not line: break
        if line.find('the norm of the test charge is') != -1:
            potential = [] # scratch the previous potential if exists.
            outcar.readline
            for i in range(num_potential_lines):
                line = outcar.readline().split("-")
                for j in range(min(5, num_atoms - 5 * i)):
                    potential.append(-float(line[j + 1].split()[0]))

    outcar.close()
    return np.array(potential)

def eigenvalues(outcar_name, num_bands, num_kpts, ispin):
    outcar = open(outcar_name,"r")

    while True:
        line = outcar.readline()
        if not line: break
        if line.find('k-point   ') != -1 or line.find('k-point     1 :') != -1:
            eigenvalues = np.zeros((ispin, num_kpts, num_bands), dtype=float)
            occupations = np.zeros((ispin, num_kpts, num_bands), dtype=float)
            outcar.readline()
            for i in range(ispin):
                for k in range(num_kpts):
                    for b in range(num_bands):
                        (eigenvalues[i, k, b], occupations[i, k, b]) =\
                          [float(j) for j in outcar.readline().split()[1:3]]
                    outcar.readline()
                    outcar.readline()
                    outcar.readline()
                outcar.readline()
                outcar.readline()
                     
    outcar.close()
    return eigenvalues, occupations
    
#num_bands = get_val("OUTCAR-finish", "NBANDS")
#num_kpts = get_val("OUTCAR-finish", "NKPTS")
#ispin = get_val("OUTCAR-finish", "ISPIN")
#eigenvalues, occupations = eigenvalues("OUTCAR-finish", num_bands, num_kpts, ispin)
#print eigenvalues, occupations
#print potential("OUTCAR-finish", num_atoms)
#print get_val("OUTCAR-1", "LNONCOLLINEAR")
