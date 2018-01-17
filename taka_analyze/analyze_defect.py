#!/bin/env python3                                                                                                                       

import sys
import os
import shutil
import argparse
from pymatgen.io.vasp.outputs import Outcar
import json
from monty.json import MontyEncoder

#outcar_path=sys.argv[1]
#o=Outcar(outcar_path)
#fw = open("defect.json", 'w')
#d = {}
#d["dielectric_tensor"] = o.dielectric_tensor
#d["dielectric_ionic_tensor"] = o.dielectric_ionic_tensor
#json.dump(d, fw, indent=2, cls=MontyEncoder)

def append_to_json(file_path, key, val):
    print("append_to_json")
    print(file_path)
    print(os.path.exists(file_path))
    with open(file_path) as f:
        d = json.load(f)
    print(d)
    if key in d:
        raise ValueError(f"Key( {key} ) is already in file( {file_path} )")
    d[key] = val
    backup_path = file_path + "._backup_by_analyze_defect_py"
    shutil.copy(file_path, backup_path)
    os.remove(file_path)
    fw = open(file_path, 'w')
    json.dump(d, fw, indent=2, cls=MontyEncoder)
    os.remove(backup_path)

def add_refpot_to_json(outcar_path, json_path):
    try:
        outcar = Outcar(outcar_path)
    except:
        print("Failed to read reference potential.\
               Specify OUTCAR of perfect structure by -p option.")
        sys.exit()
        print("read_reference pot") 
    ref_pot = outcar.electrostatic_potential
    append_to_json(json_path, "reference_potential", ref_pot)
    print(ref_pot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--all", action="store_true",
                        help="All procedure mode. (try 4.a-1, 4.a-2, and 4.a-3 of README.txt)")
    parser.add_argument("-r", "--reference_pot", action="store_true",
                        help="Reference_pot mode.\
                              Add a reference potential of perfect supercell to correction.json with\
                              (4.a-1 in README.txt)")
    parser.add_argument("-p", "--perfect_ref", dest="perfect_dir", type=str, default="perfect",
                        help="Directry name of calculation of perfect structure for reference.\
                              Needed To add a reference potential of perfect supercell \
                              to correction.json with.\
                              (4.a-1 in README.txt)")
    parser.add_argument("-i", "--info_defect", action="store_true",
                        help="Information_of_defect mode.\
                              Construct a full information JSON file(defect.json) of \
                              the defect in each directory.\
                              (4.a-2 in README.txt)")
    parser.add_argument("-v", "--visualize_structure", action="store_true",
                        help="Visualize_structure mode.\
                              Make local structure files for visualization.\
                              (4.a-3 in README.txt)")
    parser.add_argument("-d", "--defect_dir", dest="defect_dir", type=str,
                        help="Directry name of calculation of structure with defect.\
                              If you want to analyze one of defect calculations, specify with this option.\
                              Otherwise, procedure will done with all results of defect calculations.\
                              (4.a-2,3 in README.txt)")
    opts = parser.parse_args()
    if opts.all or opts.reference_pot:
        outcar_path = opts.perfect_dir + "/OUTCAR-final"
        json_path = "./correction.json"
        ref_pot = add_refpot_to_json(outcar_path, json_path)
    if opts.all or opts.info_defect:
        pass
        #info_defect()
    if opts.all or opts.visualize_structure:
        pass
        #visualize_structure()
        
    

