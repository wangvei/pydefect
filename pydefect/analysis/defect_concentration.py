#!/usr/bin/env python

import numpy as np
import optparse
import parse_eigenval as peig
import sys
from statistical_distribution import fermi_dirac, maxwell_boltzmann
from physical_constants import BoltzmannConstant, EV
# BoltzmannConstant [J/K]
# EV [J/eV]

"""
This program calculates the defect and carrier concentrations iteratively.
"""
def valence_conduction_eigenval(eigenvalues, kpt_weight, fermi_level, 
                                     temperature, ispin=1, doscar=False):
    """
    Fermi level is assumed to locate in between valence band maximum and 
    conduction band minimum. Currently, only non-magnetic host can be handled.
    eigenvalues = [[energy, weight], [energy, weight],...]
    """

    valence_band = []
    conduction_band = []

    for k in eigenvalues:
        v = eigenvalues[k]
        if v < fermi_level:
            valence_band.append([v, kpt_weight[k[0]]])
        else:
            conduction_band.append([v, kpt_weight[k[0]]])
    valence_band = np.array(sorted(valence_band)) # [[energy, weight], ....]
    conduction_band = np.array(sorted(conduction_band))

#    electron_concentration = np.sum(
#      [e[1] * fermi_dirac(e[0], fermi_level, temperature) 
#                            for e in conduction_band])
#    hole_concentration = np.sum( 
#      [e[1] * (1.0 - fermi_dirac(e[0], fermi_level, temperature))
#                            for e in valence_band])

    return valence_band, conduction_band

parser = optparse.OptionParser()
parser.add_option("-d", "--defect_list", 
                  dest="defect_list", 
                  type="string", 
                  default=False, 
                  help="List of defect energies at the VBM. \
                  First line is comment. \
                  Example: \
                  name charge formation_energy@VBM  number_of_site \
                  VZn    2     1.411                   1",
                  metavar="FILE")
parser.add_option("-i", "--intrinsic", 
                  dest="intrinsic", 
                  action="store_true",
                  default=False,
                  help="Estimate the intrinsic carriers.", 
                  metavar="FILE")
parser.add_option("-q", "--quench", 
                  dest="quench", 
                  type="string", 
                  default=False,
                  help="quench.", metavar="FILE")
parser.add_option("-t", "--temperature", 
                  dest="temperature", 
                  type="float", 
                  default="300", 
                  help="Temperature. Default:300K.")
parser.add_option("-v", "--volume", 
                  dest="volume", 
                  type="float", 
                  help="Volume of the u.c. used for counting the number of \
                        sites in Angstrom^3.")
#parser.add_option("-c", "--concentration", 
#                  dest="concentration", 
#                  default=False,
#                  help="")
parser.add_option("-e", 
                  "--eigenval", 
                  dest="eigenval", 
                  type="string", 
                  default="EIGENVAL", 
                  help="EIGENVAL type format file.", 
                  metavar="FILE")
parser.add_option("-f", 
                  "--fermi", 
                  dest="fermi", 
                  type="float", 
                  help="Fermi level wrt the vbm.")
parser.add_option("--offset", 
                  dest="offset", 
                  type="float", 
                  default=0.0, 
                  help="Shift the eigenvalues of EIGENVAL file.\
                      This is useful for using GW EIGENVAL file since \
                      the change of PAW potential shift absolute eigenvalues.")
parser.add_option("--doscar", 
                  dest="doscar", 
                  type="string", 
                  default="DOSCAR", 
                  help="DOSCAR type format file.", 
                  metavar="FILE")

opts, args = parser.parse_args()
nr_electron, kpt_coordinate, kpt_weight, eigenvalues, ispin =\
                                peig.parse_EIGENVAL(opts.eigenval)
temperature = opts.temperature
volume = opts.volume

# Only relative cbm is meaningful
vbm_index, cbm_index = peig.vbm_cbm(nr_electron, kpt_weight, eigenvalues, ispin)
vbm = eigenvalues[vbm_index] 
cbm = eigenvalues[cbm_index]

valence_band, conduction_band = valence_conduction_eigenval(
                     eigenvalues, kpt_weight, (vbm+cbm)/2.0, temperature, ispin)

if opts.fermi:
#    e = opts.fermi
    e = opts.fermi + vbm
# we need to multiply 2 for both spin up and down.
    p = np.sum(2 * [v[1] * (1.0 - fermi_dirac(v[0], e, temperature))
                                                 for v in valence_band])
    n = np.sum(2 * [c[1] * fermi_dirac(c[0], e, temperature)
                                                 for c in conduction_band])
    print "p: %.8e, n: %.8e, net: %.8e" % \
        (p / volume * 1.0e24, n / volume * 1.0e24, (p - n) / volume * 1.0e24) 
    sys.exit(0)

# Collect the defect information
defect_name = []
defect_charge = []
defect_energy = []
defect_num_site = []

defect = open(opts.defect_list, "r")
defect.readline() # skip first supercell_comment line
for i, d in enumerate(defect):
    d_each = d.split()
    if d_each == []: continue
#    print d_each
    defect_name.append(d_each[0])
    defect_charge.append(int(d_each[1]))
    defect_energy.append(float(d_each[2]))
    defect_num_site.append(int(d_each[3]))

if opts.quench:
    fixed_defect_concentration = []
    quench_file = open(opts.quench, "r")
    quench_file.readline() # skip first supercell_comment line
    for i in range(len(defect_name)):
        each_concentration = quench_file.readline().split()[2]
        fixed_defect_concentration.append(
                                   float(each_concentration) * volume * 1.0e-24)
    # Sum up the total concentration for each defect species.
    species_index = {}
    for i, name in enumerate(defect_name):
        if name in species_index:
            species_index[name].append(i) 
        else:
            species_index[name] = [i]

defect_concentration = np.zeros(len(defect_name), dtype=float)
mesh = (cbm - vbm) / 2
e = (vbm + cbm) / 2 # start to find the scf answer from the middle of the band gap.
max_iteration = 100

for iteration in range(max_iteration):
    p = np.sum([v[1] * (1.0 - fermi_dirac(v[0], e, temperature))
                                                 for v in valence_band])
    n = np.sum([c[1] * fermi_dirac(c[0], e, temperature)
                                                 for c in conduction_band])

    for i, d in enumerate(defect_energy):
        each_defect_energy = d + (e - vbm) * defect_charge[i]
        defect_concentration[i] =\
                        maxwell_boltzmann(each_defect_energy, temperature)\
                      * defect_num_site[i]

    if opts.quench:
        for k in species_index.values():
            fixed_concentration_sum = \
                              np.sum([fixed_defect_concentration[i] for i in k])
            distribution_sum = np.sum([defect_concentration[i] for i in k])
            for i in k:
                defect_concentration[i] = \
            fixed_concentration_sum * defect_concentration[i] / distribution_sum

    if opts.intrinsic:
        defect_concentration = np.zeros(len(defect_name), dtype=float)

    sum_charge = p - n + np.dot(defect_concentration, defect_charge)
    """
    In case the Fermi level locates in between vbm and cbm, the common ration 
    0.5 is sufficient. Otherwise, higher common ratio is essential, and so 0.75 
    is set here. 
    """
    mesh *= 0.75
    e = e + np.sign(sum_charge) * mesh
#    print "Sum of chargies %.8e" % (sum_charge / volume * 1.0e24)
#    print "p: %.8e, n: %.8e, net: %.8e" % \
#        (p / volume * 1.0e24, n / volume * 1.0e24, (p - n) / volume * 1.0e24) 

    # This line controles the accuracy.
    if np.abs(sum_charge / np.amax([np.amax(defect_concentration),n,p])) \
                                                                      < 1.0e-4:
#                                                                      < 1.0e-8:
        break

    if iteration == max_iteration - 1:
        print "Scf has not been reached. Bye."
        sys.exit(1)

print "Sum of chargies %.8e" % (sum_charge / volume * 1.0e24)
for i, d in enumerate(defect_concentration):
    print "%s  %3i  %.8e" % (defect_name[i], defect_charge[i], d / volume * 1.0e24)
print "p: %.8e, n: %.8e, net: %.8e" % \
        (p / volume * 1.0e24, n / volume * 1.0e24, (p - n) / volume * 1.0e24) 
        # 1 A^-3 = 1e24 cm^-3
print "Fermi level from VBM %10.7f" % (e - vbm)

sys.exit(0)
