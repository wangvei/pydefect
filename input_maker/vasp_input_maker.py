import os
import shutil
import ruamel.yaml as yaml
from pymatgen.io.vasp.inputs import Potcar
from defect_in import DefectSetting
from pydefect.input_maker.defect import get_nions
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
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except IOError:
        print('.pydefect.yaml cannot be opened.')

    for k, v in d.items():
        if k == "DEFAULT_POTCAR":
            potcar_dir = v

    if not potcar_dir:
        raise ValueError('DEFAULT_POTCAR is not set in .pydefect.yaml')

    return potcar_dir


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


def get_charge(potcar, nions, charge):
    """
    Return total charge from POTCAR file, number of ions, and charge state.
    """
    p = Potcar.from_file(potcar)
    # check only the length of potcar and nions.
    if not len(p) == len(nions):
        raise ValueError("Size of elements in POTCAR file is different")
    return sum([v.nelectrons * nions[i] for i, v in enumerate(p)]) - charge


class VaspDefectInputSetMaker(DefectInputSetMaker):

    def __init__(self, defect_setting, particular_defects="",
                 incar="INCAR", kpoints="KPOINTS"):

        # construct self._defect_setting and self._defect_set
        super().__init__(defect_setting, particular_defects)

        self._incar = incar
        self._kpoints = kpoints
        # from defect_setting
        self._perfect_structure = defect_setting.structure
        self._irreducible_sites = defect_setting.irreducible_sites
        self._interstitial_coords = defect_setting.interstitial_coords
        self._distance = defect_setting.distance
        self._cutoff = defect_setting.cutoff
        self._symprec = defect_setting.cutoff

        # Construct defect input files.
        self._make_perfect_input()
        for d in self._defect_set:
            self._make_defect_input(d)

    def _make_perfect_input(self):
        dir_name = "perfect/"
        if os.path.exists(dir_name):
            print_already_exist(dir_name)
        else:
            print_is_being_constructed(dir_name)
            os.makedirs(dir_name)
            self._defect_setting.structure.to(filename=dir_name + "POSCAR")
            shutil.copyfile(self._incar, dir_name + "INCAR")
            shutil.copyfile(self._kpoints, dir_name + "KPOINTS")
            elements = self._defect_setting.structure.symbol_set
            make_potcar(dir_name, elements, potcar_dir())


    def _make_defect_input(self, defect_name):
        # Construct: defect_structure, defect_coords, defect_index
        dir_name = defect_name + "/"

        if os.path.exists(defect_name):
            print_already_exist(dir_name)
        else:
            print_is_being_constructed(dir_name)
            os.makedirs(dir_name)

            # Constructs three POSCAR-type files
            # POSCAR-Initial: POSCAR with a defect
            # POSCAR-DispInitial: perturbed defect POSCAR near the defect
            # POSCAR: = POSCAR-DispInitial if exists, otherwise POSCAR-Initial
            a = DefectMaker(defect_name, self._perfect_structure,
                            self._irreducible_sites, self._interstitial_coords)
            d = a.defect
            d.to_json_file(dir_name + "defect.json")
            d.initial_structure.to(filename=dir_name + "POSCAR-Initial")

            if not self._distance == 0.0:
                perturbed_defect_structure, perturbed_sites = \
                    perturb_around_a_point(d.initial_structure, d.defect_coords,
                                           self._cutoff, self._distance)
                perturbed_defect_structure.to(
                    filename=dir_name + "POSCAR-DispInitial")
                shutil.copyfile(
                    dir_name + "POSCAR-DispInitial", dir_name + "POSCAR")
            else:
                shutil.copyfile(
                    dir_name + "POSCAR-Initial", dir_name + "POSCAR")

            # Construct POTCAR file
            elements = d.initial_structure.symbol_set
            make_potcar(dir_name, elements, potcar_dir())
            # Construct INCAR file
            shutil.copyfile(self._incar, dir_name + "INCAR")
            nions = get_nions(d.initial_structure)
            nelect = get_charge(dir_name + "POTCAR", nions, d.charge)
            with open(dir_name + 'INCAR', 'a') as fa:
                fa.write('NELECT = ' + str(nelect))
            # Construct KPOINTS file
            shutil.copyfile(self._kpoints, dir_name + "KPOINTS")


def main():
    import argparse
    parser = argparse.ArgumentParser(
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--defectin", dest="defectin",
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
    defect_setting = DefectSetting.from_defect_in(poscar=opts.dposcar, 
                                                  defect_in_file=opts.defectin)
    if opts.add:
        for d in opts.add:
            a = VaspDefectInputMaker(d, defect_setting, opts.incar, opts.kpoints)
            if a.is_directory:
                print_already_exist(d)
            else:
                print_is_being_constructed(d)
                a.constructor()
    else:
        VaspDefectInputSetMaker(defect_setting, particular_defects=opts.add,
                                incar=opts.incar, kpoints=opts.kpoints)


if __name__ == "__main__":
    main()
