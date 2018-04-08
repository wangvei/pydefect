# -*- coding: utf-8 -*-

import os
import ruamel.yaml as yaml
import shutil

from pymatgen.io.vasp.inputs import Potcar

from pydefect.core.defect_entry import get_num_atoms_for_elements
from pydefect.core.dft_results import defect_center
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting
from pydefect.input_maker.input_maker import \
    DefectMaker, DefectInputSetMaker, print_already_exist, \
    print_is_being_constructed,  perturb_neighbors

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pydefect.yaml")


def potcar_dir():
    """    
    Return the name of POTCAR file directory.
    SETTINGS_FILE needs to be defined in the same module.
    """
    pydefect_yaml = None
    potcar_director_path = None
    try:
        with open(SETTINGS_FILE, "r") as f:
            pydefect_yaml = yaml.safe_load(f)
    except IOError:
        print('.pydefect.yaml cannot be opened.')

    for k, v in pydefect_yaml.items():
        if k == "DEFAULT_POTCAR":
            potcar_director_path = v

    try:
        return potcar_director_path
    except ValueError:
        print('DEFAULT_POTCAR is not set in .pydefect.yaml')


def make_potcar(dirname, elements, default_potcar_dir):
    """    
    Write POTCAR with a sequence of given elements names at *dirname*.
    So far, only default POTCAR files are supported.    
    """    
    constructed_potcar = os.path.join(dirname, "POTCAR")

    with open(constructed_potcar, 'w') as potcar:
        for e in elements:
            potcar_e = "POTCAR_" + e
            potcar_file_name = os.path.join(default_potcar_dir, potcar_e)

            with open(potcar_file_name) as pot:
                potcar.write(pot.read())


def get_num_electrons_from_potcar(potcar, nions, charge):
    """
    Return total charge from POTCAR file, number of ions, and charge state.
    """
    p = Potcar.from_file(potcar)
    # check only the length of potcar and nions.
    if not len(p) == len(nions):
        raise ValueError("Size of elements in POTCAR file is different")

    return sum([v.nelectrons * nions[i] for i, v in enumerate(p)]) - charge


class VaspDefectInputSetMaker(DefectInputSetMaker):

    def __init__(self, defect_initial_setting, filtering_words=None,
                 particular_defects=None, incar="INCAR", kpoints="KPOINTS"):

        # make self._defect_initial_setting and self._defect_name_set
        super().__init__(defect_initial_setting, filtering_words,
                         particular_defects)

        self._incar = incar
        self._kpoints = kpoints

        self.make_input()

    def _make_perfect_input(self):
        if os.path.exists("perfect"):
            print_already_exist("perfect")
        else:
            print_is_being_constructed("perfect")
            os.makedirs("perfect")
            self._defect_initial_setting.structure.to(
                filename=os.path.join("perfect", "POSCAR"))
            shutil.copyfile(self._incar, os.path.join("perfect", "INCAR"))
            shutil.copyfile(self._kpoints, os.path.join("perfect", "KPOINTS"))
            elements = self._defect_initial_setting.structure.symbol_set
            make_potcar("perfect", elements, potcar_dir())

    def _make_defect_input(self, defect_name):

        #TODO: check if the defect_name is proper or not.
        if os.path.exists(defect_name):
            print_already_exist(defect_name)
        else:
            print_is_being_constructed(defect_name)
            os.makedirs(defect_name)

            # Constructs three POSCAR-type files
            # POSCAR-Initial: POSCAR with a defect
            # POSCAR-DisplacedInitial: POSCAR with perturbation near the defect
            # POSCAR: POSCAR-DisplacedInitial when neighboring atoms are
            #         perturbed, otherwise POSCAR-Initial
            d = DefectMaker(
                defect_name,
                self._defect_initial_setting.structure,
                self._defect_initial_setting.irreducible_sites,
                self._defect_initial_setting.interstitial_coords).defect
            d.to_json_file(os.path.join(defect_name, "defect.json"))
            d.initial_structure.to(
                filename=os.path.join(defect_name, "POSCAR-Initial"))

            if not self._defect_initial_setting.distance == 0.0:
                center = defect_center(d.initial_structure, d)
                perturbed_defect_structure, perturbed_sites = \
                    perturb_neighbors(d.initial_structure,
                                      center,
                                      self._defect_initial_setting.cutoff,
                                      self._defect_initial_setting.distance)
                perturbed_defect_structure.\
                    to(filename=os.path.join(defect_name,
                                             "POSCAR-DisplacedInitial"))
                shutil.copyfile(
                    os.path.join(defect_name, "POSCAR-DisplacedInitial"),
                    os.path.join(defect_name, "POSCAR"))
            else:
                shutil.copyfile(
                    os.path.join(defect_name, "POSCAR-Initial"),
                    os.path.join(defect_name, "POSCAR"))

            # Construct POTCAR file
            elements = d.initial_structure.symbol_set
            make_potcar(defect_name, elements, potcar_dir())

            # Construct INCAR file
            shutil.copyfile(self._incar, os.path.join(defect_name, "INCAR"))
            nions = get_num_atoms_for_elements(d.initial_structure)
            nelect = get_num_electrons_from_potcar(
                 os.path.join(defect_name, "POTCAR"), nions, d.charge)

            with open(os.path.join(defect_name, 'INCAR'), 'a') as fa:
                fa.write('NELECT = ' + str(nelect))

            # copy KPOINTS file
            shutil.copyfile(self._kpoints, os.path.join(defect_name, "KPOINTS"))


def main():
    import argparse
    parser = argparse.ArgumentParser(
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--defect_in", dest="defect_in",
                        default="defect.in", type=str, help="defect.in name.")
    parser.add_argument("-p", "--dposcar", dest="dposcar", default="DPOSCAR",
                        type=str, help="DPOSCAR name.")
    parser.add_argument("--incar", dest="incar", default="INCAR", 
                        type=str, help="INCAR name.")
    parser.add_argument("--kpoints", dest="kpoints", default="KPOINTS", 
                        type=str, help="KPOINTS name.")
    parser.add_argument("--add", dest="add", type=str, default=None, nargs="+",
                        help="Particular defect name added.")

    opts = parser.parse_args()
    defect_initial_setting = DefectInitialSetting.\
        from_defect_in(poscar=opts.dposcar, defect_in_file=opts.defect_in)

    # if opts.add:
    #     for d in opts.add:
    #         a = VaspDefectInputSetMaker(defect_initial_setting,
    #                                     opts.incar,
    #                                     opts.kpoints)
    #         if a.is_directory:
    #             print_already_exist(d)
    #         else:
    #             print_is_being_constructed(d)
    #             a.constructor()
    # else:
    VaspDefectInputSetMaker(defect_initial_setting=defect_initial_setting,
                            filtering_words=opts.add,
                            incar=opts.incar,
                            kpoints=opts.kpoints)


if __name__ == "__main__":
    main()
