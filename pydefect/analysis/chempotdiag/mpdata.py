import json
import os
import shutil
from glob import glob
from itertools import combinations

from monty.serialization import loadfn
from pymatgen import MPRester, Composition, Structure
from pymatgen.io.vasp import Poscar

from pydefect.analysis.chempotdiag.chem_pot_diag import molecule_file_names
from pydefect.analysis.chempotdiag.gas import Gas

# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class PriorInfo:

    def __init__(self, energy_per_atom=None, band_gap=None,
                 total_magnetization=None, data_source=None,
                 mag_threshold=0.001, band_gap_threshold=0.1):
        self._energy_per_atom = energy_per_atom
        self._band_gap = band_gap
        self._total_magnetization = total_magnetization
        self._data_source = data_source
        self._mag_threshold = mag_threshold
        self._band_gap_threshold = band_gap_threshold

    def as_dict(self):
        d = {"energy_per_atom": self._energy_per_atom,
             "band_gap": self._band_gap,
             "total_magnetization": self._total_magnetization,
             "data_source": self._data_source}
        return d

    def dump_json(self, filename="prior_info.json"):
        with open(filename, "w") as fw:
            json.dump(self.as_dict(), fw, indent=2)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectEntry class object from a dictionary.
        """
        # The keys need to be converted to integers.

        return cls(d["energy_per_atom"], d["band_gap"],
                   d["total_magnetization"])

    @classmethod
    def load_json(cls, filename="prior_info.json"):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def energy_per_atom(self):
        return self._energy_per_atom

    @property
    def band_gap(self):
        return self._band_gap

    @property
    def total_magnetization(self):
        return self._total_magnetization

    @property
    def is_magnetic(self):
        if self._total_magnetization > self._mag_threshold:
            return True

    @property
    def is_band_gap(self):
        if self._band_gap > self._band_gap_threshold:
            return True


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
    # make directory
    chem_sys = "-".join(elements)
    chem_sys_dir = os.path.join(dir_path, chem_sys)
    if not os.path.exists(chem_sys_dir):
        os.mkdir(os.path.join(dir_path, chem_sys))

    # get molecules
    if adds_molecule:
        molecules_elements = [m.composition for m in Gas]
        for me, file_name in zip(molecules_elements, molecule_file_names):
            if set([str(e) for e in me.elements]) < set(elements):
                comp_name = str(me)
                name_dict = {"N1 H3": "NH3", "N1 O2": "NO2", "O1 H2": "H2O"}
                if comp_name in name_dict.keys():
                    comp_name = name_dict[comp_name]
                dirname = comp_name + molecule_dir_name
                dirname2 = os.path.join(chem_sys_dir, dirname)
                if not os.path.exists(dirname2):
                    os.mkdir(dirname2)
                    shutil.copyfile(file_name+"/POSCAR", dirname2+"/POSCAR")

    # get from mp
    for i in range(len(elements)):
        for els in combinations(elements, i+1):
            els_str = "-".join(els)
            all_materials = mp_rester.get_data(els_str, data_type="vasp")
            materials = [material for material in all_materials
                         if material["e_above_hull"] < criterion_e_above_hull]
            materials_to_output = []
            comp_stable = {}
            if gets_poly:
                materials_to_output = materials
            else:
                for material in materials:
                    c = Composition(material["full_formula"])
                    if gets_poly:
                        materials_to_output.append(material)
                    else:
                        if c not in comp_stable or \
                                material["e_above_hull"] <\
                                comp_stable[c.reduced_formula]["e_above_hull"]:
                            comp_stable[c.reduced_formula] = material
                        else:
                            break
                materials_to_output = comp_stable.values()
            for material in materials_to_output:
                # remove solids when directory of molecules exist
                length = len(molecule_dir_name)
                exist_molecules_reduced_formulas = \
                    [Composition(path.split("/")[-1][:-length]).reduced_formula
                     for path in glob(chem_sys_dir + r"/*" + molecule_dir_name)]
                if adds_molecule:
                    reduced_formula = \
                        Composition(material["full_formula"]).reduced_formula
                    if reduced_formula in exist_molecules_reduced_formulas:
                        break
                # make directory and files
                mp_id = material["material_id"]
                dirname = material["full_formula"] + "_" + mp_id
                dirname2 = os.path.join(chem_sys_dir, dirname)
                if not os.path.exists(dirname2):
                    os.mkdir(dirname2)
                    structure = Structure.from_str(material["cif"], "cif")
                    poscar = Poscar(structure)
                    poscar.write_file(os.path.join(dirname2, "POSCAR"))
                    # dump json
                    keys_to_get = ["energy_per_atom", "band_gap",
                                   "total_magnetization"]
                    d = {k: material[k] for k in keys_to_get}
                    json_path = os.path.join(dirname2, "prior_info.json")
                    with open(json_path, "w") as fw:
                        json.dump(d, fw)


