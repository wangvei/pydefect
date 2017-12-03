#!/usr/bin/env python

import numpy as np
import optparse 
import physical_constants as pc
import sys  

def parsePOTCAR(potcar_name):
    potcar = open(potcar_name,"r")

    functional = []
    pot_name = []
    dates = []
    orbitals = []
    num_electrons = []    
    rcore = []
    rwigs = []
    enmax = []

    potcar.readline()
    num_electrons.append(float(potcar.readline().split()[0]))

    while True:

        line = potcar.readline().split()

        if line == []: continue

        if line[0] == "End": 
            try:
                potcar.readline()
                num_electrons.append(float(potcar.readline().split()[0]))
                continue
            except:
                potcar.close()
                break

        if line[0] == "VRHFIN":
            orbitals.append(line[2])

        if line[0] == "TITEL":
            if line[2] == "PAW":
                functional.append("LDA")
            elif line[2] == "PAW_PBE":
                functional.append("PBE")

            pot_name.append(line[3]) 
            dates.append(line[4]) 

        elif line[0] == "RCORE":
            rcore.append(float(line[2])) 

        elif line[0] == "RWIGS":
            rwigs.append(float(line[2].split(";")[0])) 

        elif line[0] == "ENMAX":
            enmax.append(float(line[2].split(";")[0])) 

    elements = [i.split("_")[0] for i in pot_name]

    return functional, pot_name, dates, orbitals,\
                                    num_electrons, rcore, rwigs, enmax, elements

parser = optparse.OptionParser()
parser.add_option("--potcar",
                  dest="potcar",
                  type="string",
                  default="POTCAR",
                  help="POTCAR file.",
                  metavar="FILE")
parser.add_option("--rwigs",                                                      
                  dest="rwigs",                                                   
                  action="store_true",                                           
                  default=False,                                                 
                  help="Print RWIGS [au] in POTCAR file.")     
parser.add_option("--addrwigs",                                                      
                  dest="addrwigs",                                                   
                  type="float",                                           
                  default=False,                                                 
                  help="Add one more RWIGS. This can be used for empty sphere.")

if __name__ == "__main__":
    (opts, args) = parser.parse_args()

    if opts.rwigs:                                                                
        print "RWIGS =", 
        for i in parsePOTCAR(opts.potcar)[6]:
            print i,
        if opts.addrwigs: 
            print opts.addrwigs,
        print "" 
        sys.exit(0)  

    a = parsePOTCAR(opts.potcar)
    print "      Functional:", tuple(a[0])
    print " Potential names:", tuple(a[1])
    print "           Dates:", tuple(a[2])
    print "Valence orbitals:", tuple(a[3])
    print "  Num. electrons:", tuple(a[4])
    print "Core radii  [au]:", tuple(a[5])
    print "Core radii   [A]:", tuple([i * pc.Bohr for i in a[5]])
    print "RWIGS radii [au]:", tuple(a[6])
    print "RWIGS radii  [A]:", tuple([i * pc.Bohr for i in a[6]])
    print "      ENMAX [eV]:", tuple(a[7])
    print "      ENMAX [eV]:", tuple(a[7])
    print "        Elements:", tuple(a[8])
    sys.exit(0)
