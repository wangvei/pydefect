#!/usr/bin/env python

import numpy as np
import optparse
import sys

parser = optparse.OptionParser()
parser.add_option("-p", "--potential",
                  dest="potential",
                  type="string",
                  help="potential file with the following 8 columns.\
                  symbol # fractional_coordinates(3 columns) \
                    distance_from_a_defect abinitio_potential model_potential",
                  metavar="FILE")
parser.add_option("-d","--distance", 
                  dest="distance", 
                  type="float", 
                  help="Distance from which the potential shift is estimated.")
parser.add_option("-c", "--charge",
                  dest="charge",
                  type="float",
                  help="Charge of gaussian [e].")

(opts, args) = parser.parse_args()

potentials = open(opts.potential)
potential_diff = []
for line in potentials:
    distance, abinitio_potential, model_potential = \
                                [float(i) for i in line.split()[5:8]]
    if distance > opts.distance:
        # potential(ab-initio) - potential(model)
        # Note that the potential written in vasp output is electron's one.
        # Here, all potential is treated within the convetional electrostatics.
        potential_diff.append(abinitio_potential - model_potential)

average_potential_diff = np.mean(potential_diff)

# correction is -q * \delta V
print "Potential difference: %10.7f [eV] \n alignment-like term: %10.7f [eV]" %\
            (average_potential_diff, -1.0 * average_potential_diff * opts.charge)
sys.exit(0)
