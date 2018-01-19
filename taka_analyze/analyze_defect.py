#!/bin/env python3                                                                                                                       

import sys
import os
import shutil
import glob
import argparse
import numpy as np
from copy import deepcopy
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Outcar
import json
from monty.json import MontyEncoder
from itertools import product

def append_to_json(file_path, key, val):
    if not os.path.exists(file_path):
        sys.exit( "File( {0} ) doesn't exist.".format([file_path]))
    with open(file_path) as f:
        d = json.load(f)
    if key in d:
        sys.exit( "Key( {0} ) is already in file( {1} )".format([key,file_path]))
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
    ref_pot = outcar.electrostatic_potential
    append_to_json(json_path, "reference_potential", ref_pot)

def calc_min_distance_and_its_v2coord(v1, v2, axis):
    candidate_list = []
    for x, y, z in product((-1, 0, 1), repeat=3):
        delta_vect = np.dot(axis, np.array([x, y, z]))
        distance = np.linalg.norm(delta_vect+v2-v1)
        candidate_list.append((distance, delta_vect+v2))
    return min(candidate_list, key = lambda t: t[0])

def read_defect_pos(d):
    if(isinstance(d["defect_position"], list)):
        defect_pos_frac = np.array(d["defect_position"])
    elif(isinstance(d["defect_position"], int)):
        defect_pos_frac = poscar_final.structure.frac_coords[d["defect_position"]]
    else:
        raise TypeError("Failed to read defect position. (not list and int)")
    return defect_pos_frac

def complete_defect_json(dirname):
    poscar_initial = Poscar.from_file(dirname + "/POSCAR-initial")
    poscar_final = Poscar.from_file(dirname + "/POSCAR-final")
    coords_initial = poscar_initial.structure.cart_coords
    coords_final = poscar_final.structure.cart_coords
    elements = [e.name for e in poscar_final.structure.species]
    axis = poscar_initial.structure.lattice.matrix
    axis_inv = np.linalg.inv(axis)
    with open(dirname+"/defect.json") as f:
        d = json.load(f)
    defect_pos_frac = read_defect_pos(d)
    defect_pos = np.dot(axis, defect_pos_frac)
    distance_list, displacement_list, angle_list = [], [], []
    with open(dirname+"/structure.txt", 'w') as f:
        f.write("#Defect position (frac)\n")
        #f.write(f"#{defect_pos_frac[0]: .6f} {defect_pos_frac[0]: .6f} {defect_pos_frac[0]: .6f}\n")
        f.write("#{0: .6f} {1: .6f} {2: .6f}\n".format((defect_pos_frac[0],
                                                        defect_pos_frac[1],
                                                        defect_pos_frac[2])))
        f.write("#" + " " * 8 + "-" * 5 + "coordinations (frac)" \
                + "-" * 5 + " " * 3 +"dist.(init)[A]" + " " * 3\
                + "dist.(final)[A]  disp.[A]  angle[deg.]\n")
        for i, (vi, vf, e) in enumerate(zip(coords_initial, coords_final, elements)): # calculate displacement, distance_from_defect, angle
#To calculate displacement, it is sometimes needed to find atoms in neighbor atoms
#due to periodical boundary condition.
#For example, initial position is (0,0,0) can be displaced to (0.99,0.99,0.99).
#Then it is needed to calculate distance between (0,0,0) and (-0.01,-0.01,-0.01).
            disp, neighbor_vf = calc_min_distance_and_its_v2coord(vi, vf, axis)
            distance_init, _ = calc_min_distance_and_its_v2coord(defect_pos, vi, axis)
            distance, _ = calc_min_distance_and_its_v2coord(defect_pos, vf, axis)
            if disp >= 0.1:
            #if disp >= 0:
                v1 = defect_pos - vi
                v2 = neighbor_vf - vi
                cosine = np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
                angle = 180 - np.degrees(np.arccos(np.clip(cosine, -1, 1)))
            else:
                angle = "-"
            ne_vf_frac = np.round(np.dot(axis_inv, neighbor_vf),6)
            #f.write(f"{e:3s} {str(i+1).rjust(3)}  {ne_vf_frac[0]: .5f}   {ne_vf_frac[1]: .5f}   {ne_vf_frac[2]: .5f}\
            f.write("{0:3s} {1}  {2: .5f}   {3: .5f}   {4: .5f}\
       {5:.3f}            {6:.3f}         {7:.3f}         {8}\n".format(\
           (e, str(i+1).rjust(3), ne_vf_frac[0], ne_vf_frac[1], ne_vf_frac[2],
            distance_init, distance, disp, angle)))
            distance_list.append(distance)
            displacement_list.append(disp)
            angle_list.append(angle)
    json_name = dirname + "/defect.json"
    charge = int(dirname.split("_")[2][:-1])
    append_to_json(json_name, "charge", charge)
    outcar = Outcar(dirname + "/OUTCAR-final")
    total_energy = outcar.final_energy
    append_to_json(json_name, "total_energy", total_energy)
    atomic_site_pot = outcar.electrostatic_potential
    append_to_json(json_name, "atomic_site_pot", atomic_site_pot)
    append_to_json(json_name, "axis", axis.tolist())
    coords_final_frac = poscar_final.structure.frac_coords
    append_to_json(json_name, "frac_coords", coords_final_frac.tolist())
    append_to_json(json_name, "elements", elements)
    append_to_json(json_name, "displacement", displacement_list)
    append_to_json(json_name, "distance_from_defect", distance_list)
    append_to_json(json_name, "angle", angle_list)

def make_local_structure_for_visualization(dirname, rate):
    poscar_final = Poscar.from_file(dirname + "/POSCAR-final")
    with open(dirname+"/defect.json") as f:
        d = json.load(f)
    defect_pos_frac = read_defect_pos(d)
    try:
        distance_list = d["distance_from_defect"]
    except:
        sys.exit('Not found "distance_from_defect" in {0}/defect.json '.format((dirname)))
    threshold = min(distance_list) * rate
    poscar_vis = deepcopy(poscar_final) # poscar for visualization
    remove_list = []
    for i, d in enumerate(distance_list):
        if d >= threshold:
            remove_list.append(i)
        else:
            translate = -np.array(defect_pos_frac) + np.array([0.5, 0.5, 0.5])
            poscar_vis.structure.translate_sites(i, translate)
    poscar_vis.structure.remove_sites(remove_list)
    poscar_vis.write_file(dirname+"/POSCAR-local_visualize")


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
                              Used to add a reference potential of perfect supercell \
                              to correction.json with, and calculation of atomic displacement. \
                              (4.a-1, 2 in README.txt)")
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
        print("reference_pot mode")
        outcar_path = opts.perfect_dir + "/OUTCAR-final"
        json_path = "./correction.json"
        ref_pot = add_refpot_to_json(outcar_path, json_path)
    if opts.all or opts.info_defect:
        print("info_defect mode")
        if opts.defect_dir:
            complete_defect_json(opts.defect_dir)
        else:
            dirs = glob.glob("./defect/*_*_*/")
            if not dirs:
                print("Warning: No directory matched name defect_*_*_*.")
            for dirname in dirs:
                complete_defect_json(dirname)
    if opts.all or opts.visualize_structure:
        print("visual_local_structure mode")
        RATE = 1.9
        if opts.defect_dir:
            make_local_structure_for_visualization(opts.defect_dir, RATE)
        else:
            dirs = glob.glob("./defect/*_*_*/")
            if not dirs:
                print("Warning: No directory matched name defect_*_*_*.")
            for dirname in dirs:
                make_local_structure_for_visualization(dirname, RATE)
        
    

