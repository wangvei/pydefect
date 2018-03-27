import os
import ruamel.yaml as yaml
import shutil

from pymatgen.io.vasp.inputs import Potcar

from pydefect.input_maker.defect_in import DefectInitialSetting
from input_maker.defect_entry import get_nions
from pydefect.input_maker.input_maker import \
    DefectMaker, DefectInputSetMaker,  print_already_exist, \
    print_is_being_constructed,  perturb_around_a_point

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
        with open(SETTINGS_FILE, "rt") as f:
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
    with open(dirname + '/POTCAR', 'w') as potcar:
        for e in elements:
            potcar_file_name = default_potcar_dir + "/POTCAR_" + e
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

    def __init__(self, defect_initial_setting, particular_defects="", incar="INCAR",
                 kpoints="KPOINTS"):

        # make self._defect_initial_setting and self._defect_name_set
        super().__init__(defect_initial_setting, particular_defects)

        self._incar = incar
        self._kpoints = kpoints

        # Construct defect input files.
        self._make_perfect_input()
        for d in self._defect_name_set:
            self._make_defect_input(d)

    def _make_perfect_input(self):
        dir_name = "perfect/"
        if os.path.exists(dir_name):
            print_already_exist("perfect")
        else:
            print_is_being_constructed("perfect")
            os.makedirs(dir_name)
            self._defect_initial_setting.structure.to(filename=dir_name + "POSCAR")
            shutil.copyfile(self._incar, dir_name + "INCAR")
            shutil.copyfile(self._kpoints, dir_name + "KPOINTS")
            elements = self._defect_initial_setting.structure.symbol_set
            make_potcar(dir_name, elements, potcar_dir())

    def _make_defect_input(self, defect_name):
        # Construct: defect_structure, defect_coords, defect_index
        dir_name = defect_name + "/"

        if os.path.exists(defect_name):
            print_already_exist(defect_name)
        else:
            print_is_being_constructed(defect_name)
            os.makedirs(dir_name)

            # Constructs three POSCAR-type files
            # POSCAR-Initial: POSCAR with a defect
            # POSCAR-DisplacedInitial: POSCAR with perturbation near the defect
            # POSCAR: POSCAR-DisplacedInitial if exist, otherwise POSCAR-Initial
            d = DefectMaker(
                defect_name,
                self._defect_initial_setting.structure,
                self._defect_initial_setting.irreducible_sites,
                self._defect_initial_setting.interstitial_coords).defect
            d.to_json_file(dir_name + "defect.json")
            d.initial_structure.to(filename=dir_name + "POSCAR-Initial")

            if not self._defect_initial_setting.distance == 0.0:
                perturbed_defect_structure, perturbed_sites = \
                    perturb_around_a_point(d.initial_structure,
                                           d.defect_coords,
                                           self._defect_initial_setting.cutoff,
                                           self._defect_initial_setting.distance)
                perturbed_defect_structure.to(
                    filename=dir_name + "POSCAR-DisplacedInitial")
                shutil.copyfile(
                    dir_name + "POSCAR-DisplacedInitial", dir_name + "POSCAR")
            else:
                shutil.copyfile(
                    dir_name + "POSCAR-Initial", dir_name + "POSCAR")

            # Construct POTCAR file
            elements = d.initial_structure.symbol_set
            make_potcar(dir_name, elements, potcar_dir())
            # Construct INCAR file
            shutil.copyfile(self._incar, dir_name + "INCAR")
            nions = get_nions(d.initial_structure)
            nelect = get_num_electrons_from_potcar(dir_name + "POTCAR",
                                                   nions, d.charge)

            with open(dir_name + 'INCAR', 'a') as fa:
                fa.write('NELECT = ' + str(nelect))

            # copy KPOINTS file
            shutil.copyfile(self._kpoints, dir_name + "KPOINTS")


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
    parser.add_argument("--add", dest="add", type=str, default=None,
                        help="Particular defect name added.")

    opts = parser.parse_args()
    defect_initial_setting = DefectInitialSetting.\
        from_defect_in(poscar=opts.dposcar, defect_in_file=opts.defect_in)

    if opts.add:
        for d in opts.add:
            a = VaspDefectInputSetMaker(d, defect_initial_setting, opts.incar,
                                        opts.kpoints)
            if a.is_directory:
                print_already_exist(d)
            else:
                print_is_being_constructed(d)
                a.constructor()
    else:
        VaspDefectInputSetMaker(defect_initial_setting,
                                particular_defects=opts.add,
                                incar=opts.incar,
                                kpoints=opts.kpoints)


if __name__ == "__main__":
    main()
