import numpy as np
import optparse 
import parse_outcar as pout
import sys  

"""
Extract the atomic site potential written in the last part of OUTCAR file.
"""

parser = optparse.OptionParser()
parser.add_option("-o", "--outcar", 
                  dest="outcar", 
                  default="OUTCAR", 
                  help="OUTCAR file including average potential.",
                  metavar="FILE")

parser.add_option("-r", "--reference", 
                  dest="reference", 
                  default="OUTCAR", 
                  help="Reference potential file composed of 4 columns.\
                  E.g. 'Si 1 255 -10.5123',\
                  where first is atomic species, second is atom index,\
                  third is number of atoms with its index, and fourth is the\
                  reference potential.",
                  metavar="FILE")

parser.add_option("-d", "--defect",
                  dest="defect",
                  default=False,
                  help="The potential of the input number of atom is removed.\
                        E.g. 82. This is used for interstitial, anti-site, or\
                        extrinsic defect.",
                  metavar="INT")

(opts, args) = parser.parse_args()

num_atoms = pout.get_val(opts.outcar, "NIONS")
potential = pout.potential(opts.outcar, num_atoms)

if len(opts.defect.split()) == 1: 
    potential = np.delete(potential, int(opts.defect) - 1)

# Extract the reference potential from reference_potential file
reference = open(opts.reference)
reference_potential = []
for r in reference.readlines():
    num_each_site, each_reference_potential = r.split()[2:4]
    for i in range(int(num_each_site)): 
        reference_potential.append(float(each_reference_potential))

if len(potential) == len(reference_potential):
    for i in range(len(potential)): 
        # Potential is based on the conventional electrostatics.
        # Thus, the sign is opposite to vasp output.
        print " %8.4f" % (-(potential[i] - reference_potential[i]))
else:
    print "Number of atoms are different between",\
                            opts.outcar, " and ", opts.reference, "."
reference.close()
sys.exit(0)
