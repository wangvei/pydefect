# -*- coding: utf-8 -*-
import json
import os
import shutil
from glob import glob
from itertools import combinations

from pymatgen import MPRester, Composition, Structure
from pymatgen.io.vasp import Poscar

from pydefect.analysis.chempotdiag.chem_pot_diag import molecule_directory
from pydefect.analysis.chempotdiag.gas import Gas


__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


def make_vasp_inputs_from_mp(elements,
                             dir_path=os.getcwd(),
                             criterion_e_above_hull=float("inf"),
                             api_key=None,
                             gets_poly=False,
                             adds_molecule=True):
    """

    Args:
        elements(list): like ["Cu", "O"]
        dir_path(str):
        criterion_e_above_hull(float):
        api_key(str):
        gets_poly(bool):
        adds_molecule(bool):

    Returns:

    """
    molecule_dir_name = "molecule_pydefect"
    if not os.path.isdir(dir_path):
        raise NotADirectoryError(dir_path + " is not directory.")
    mp_rester = MPRester(api_key)

    # get molecules
    molecules_formula_list = []
    if adds_molecule:
        for g in Gas:
            me = g.composition
            if set([str(e) for e in me.elements]) < set(elements):
                molecules_formula_list.append(me.reduced_formula)
                comp_name = str(me)
                dirname = comp_name + molecule_dir_name
                if not os.path.exists(dirname):
                    os.mkdir(dirname)
                    shutil.copyfile("{}/{}/POSCAR".
                                    format(molecule_directory, str(g)),
                                    dirname+"/POSCAR")

    # get from mp
    for i in range(len(elements)):
        for els in combinations(elements, i+1):
            els_str = "-".join(els)
            all_materials = mp_rester.get_data(els_str, data_type="vasp")
            materials = [material for material in all_materials
                         if material["e_above_hull"] < criterion_e_above_hull]
            comp_stable = {}
            if gets_poly:
                materials_to_output = materials
            else:
                for material in materials:
                    c = Composition(material["full_formula"])
                    if c.reduced_formula in molecules_formula_list:
                        continue
                    if c.reduced_formula not in comp_stable or \
                            material["e_above_hull"] <\
                            comp_stable[c.reduced_formula]["e_above_hull"]:
                        comp_stable[c.reduced_formula] = material
                    else:
                        continue
                materials_to_output = comp_stable.values()
            for material in materials_to_output:
                # remove solids when directory of molecules exist
                length = len(molecule_dir_name)
                exist_molecules_reduced_formulas = \
                    [Composition(path.split("/")[-1][:-length]).reduced_formula
                     for path in glob(r"/*" + molecule_dir_name)]
                if adds_molecule:
                    reduced_formula = \
                        Composition(material["full_formula"]).reduced_formula
                    if reduced_formula in exist_molecules_reduced_formulas:
                        break
                # make directory and files
                mp_id = material["material_id"]
                dirname = material["full_formula"] + "_" + mp_id
                if not os.path.exists(dirname):
                    os.mkdir(dirname)
                    structure = Structure.from_str(material["cif"], "cif")
                    poscar = Poscar(structure)
                    poscar.write_file(os.path.join(dirname, "POSCAR"))
                    # dump json
                    keys_to_get = ["energy_per_atom", "band_gap",
                                   "total_magnetization"]
                    d = {k: material[k] for k in keys_to_get}
                    json_path = os.path.join(dirname, "prior_info.json")
                    with open(json_path, "w") as fw:
                        json.dump(d, fw)


