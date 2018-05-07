# -*- coding: utf-8 -*-

from enum import Enum, auto, unique
import numpy as np
import os
import shutil

from monty.serialization import loadfn
import ruamel.yaml as yaml

from pymatgen.io.vasp import Potcar, Kpoints, Incar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.defect_entry import get_num_atoms_for_elements, \
    get_num_electrons_from_potcar
from pydefect.core.supercell_dft_results import defect_center
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting, \
    element_set
from pydefect.input_maker.kpt_centering import kpt_centering
from pydefect.input_maker.input_maker import \
    DefectMaker, DefectInputSetMaker, print_is_being_removed, \
    print_already_exist, print_is_being_constructed, perturb_neighbors

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pydefect.yaml")
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
INCAR_SETTINGS_FILE = os.path.join(MODULE_DIR, "default_INCAR_setting.yaml")
DEFAULT_INCAR = os.path.join(MODULE_DIR, "default_INCAR")


def potcar_dir():
    """    
    Returns the name of POTCAR file directory.
    SETTINGS_FILE needs to be defined in the same module.
    """
    pydefect_yaml = None
    potcar_director_path = None
    try:
        with open(SETTINGS_FILE, "r") as f:
            pydefect_yaml = yaml.safe_load(f)
#    try:
#        with open(SETTINGS_FILE, "r") as fr:
#            pydefect_yaml = loadfn(fr)
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
    Writes POTCAR with a list of the given elements names.
    Now, only default POTCAR files are supported.
    """    
    with open(os.path.join(dirname, "POTCAR"), 'w') as potcar:
        for e in elements:
            potcar_e = "POTCAR_" + e
            potcar_file_name = os.path.join(default_potcar_dir, potcar_e)

            with open(potcar_file_name) as pot:
                potcar.write(pot.read())

@unique
class SupportedTask(Enum):
    structure_opt = auto()
    band = auto()
    dos = auto()
    dielectric = auto()
    dielectric_function = auto()
    competing_phase = auto()
    competing_phase_molecule = auto()
    defect = auto()

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s):
        for m in SupportedTask:
            print(m.name)
            if m.name == s:
                return m
        raise ValueError("Can't interpret supported task %s" % s)

    @staticmethod
    def hasattr(s):
        for m in SupportedTask:
            if m.name == s:
                return True
        return False


@unique
class SupportedFunctional(Enum):
    pbe = auto()
    hse06 = auto()
    hse06_pre_gga = auto()
    pbesol = auto()
    pbe_d3 = auto()

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s):
        for m in SupportedFunctional:
            if m.name == s:
                return m
        raise ValueError("Can't interpret supported functional %s" % s)

    @staticmethod
    def hasattr(s):
        for m in SupportedFunctional:
            if m.name == s:
                return True
        return False


def make_kpoints(task, dirname='.', poscar="DPOSCAR", ibzkpt=None,
                 is_metal=False, kpts_shift=None, kpt_density_opt=4,
                 kpt_density_defect=2, kpt_multiplier_factor=2):
    """
    Constructs a KPOINTS file based on default settings depending on the task.
    Args:
        task (SupportedTask Enum):
        dirname (str): Director name at which INCAR is written.
        poscar (str): POSCAR-type file name
        ibzkpt (str): IBZKPT-type file name
        is_metal (bool):
        kpt_shift (3 list): k-point shift vector

    """

    s = Structure.from_file(os.path.join(poscar))
    reciprocal_lattice = s.lattice.reciprocal_lattice

    task = SupportedTask.from_string(task)

    # check if



    if task == SupportedTask.band:
        if ibzkpt is None:
            # TODO: rewrite
            raise ValueError

    if task == SupportedTask.defect:
        kpt_mesh = [int(np.ceil(kpt_density_defect * r))
                    for r in reciprocal_lattice.abc]
    else:
        if is_metal is True and task == SupportedTask.competing_phase:
            kpt_mesh = [int(np.ceil(kpt_density_opt * r * 2))
                        for r in reciprocal_lattice.abc]
        else:
            kpt_mesh = [int(np.ceil(kpt_density_opt * r))
                        for r in reciprocal_lattice.abc]

        if task in (SupportedTask.dos,
                    SupportedTask.dielectric,
                    SupportedTask.dielectric_function):
            kpt_mesh = [k * kpt_multiplier_factor for k in kpt_mesh]

    if kpts_shift is None:
        if task == SupportedTask.defect:
            angles = reciprocal_lattice.angles

            if not kpts_shift:
                kpts_shift = []
                for i in range(3):
                    if angles[i - 2] == 90 and angles[i - 1] == 90:
                        kpts_shift.append(0.5)
                    else:
                        kpts_shift.append(0.0)

        else:
            sg = SpacegroupAnalyzer(structure=s).get_space_group_number()

            if 168 <= sg <= 194:
                angles = reciprocal_lattice.angles
                if angles[2] != 90:
                    raise ValueError("Axial in the hexagonal structure is not "
                                     "set properly.")

            kpts_shift = kpt_centering[sg]

    kpts = (tuple(kpt_mesh), )

    comment = "K-point mesh generated by pydefect."
    Kpoints(comment=comment, kpts=kpts, kpts_shift=kpts_shift).\
        write_file(os.path.join(dirname, 'KPOINTS'))


def get_max_enmax(element_names, dirname=potcar_dir()):

    enmax = []
    for e in element_names:
        potcar_file_name = os.path.join(dirname, "POTCAR_" + e)
        enmax.append(Potcar.from_file(potcar_file_name)[0].enmax)

    return max(enmax)


def make_incar(task, functional, hfscreen=None, aexx=None,
               is_magnetization=False, dirname='.', defect_in=None,
               poscar="DPOSCAR", set_pars=True):
    """
    Constructs an INCAR file based on default settings depending on the task.
    ENCUT can be determined from max(ENMAX).
    Args:
        task (SupportedTask Enum):
        functional (SupportedFunctional Enum):
        is_magnetization (bool):
        dirname (str): Director name at which INCAR is written.
        defect_in (str): defect.in file name parsed for checking elements
        poscar (str): POSCAR-type file name parsed for checking elements
        set_pars (book): if set the NPAR and KPAR at the command line

    defect.in file is parsed when exists, otherwise DPOSCAR file.
    """

    setting_file = loadfn(INCAR_SETTINGS_FILE)

    if SupportedTask.hasattr(task) is False:
        print("SupportedTask:", task, "is not a proper flag.")
    if SupportedFunctional.hasattr(functional) is False:
        print("SupportedFunctional:", functional, "is not a proper flag.")

    setting = {}
    for i in (task, functional):
        try:
            setting.update(setting_file[i])
        except IOError:
            print('default_INCAR_setting.yaml cannot be opened')

    if is_magnetization:
        setting["ISPIN"] = 2

    # TODO: read user setting yaml file.
    if defect_in:
        defect_initial_setting = \
            DefectInitialSetting.from_defect_in(poscar, defect_in)
        setting["ENCUT"] = \
            str(get_max_enmax(element_set(defect_initial_setting)))
    else:
        setting["ENCUT"] = input("Input ENCUT:")

    if set_pars:
        setting["NPAR"] = input("Input NPAR:")

    if set_pars:
        setting["KPAR"] = input("Input KPAR:")

    if hfscreen:
        setting["HFSCREEN"] = hfscreen

    if aexx:
        setting["AEXX"] = aexx

    incar = os.path.join(dirname, 'INCAR')
    Incar(setting).write_file(filename=incar)

    if SupportedFunctional.hasattr(functional):
        incar_pre_gga = os.path.join(dirname, 'INCAR-pre_gga')

        setting.pop("LHFCALC")
        setting.pop("HFSCREEN")
        setting.pop("PRECFOCK")
        setting.pop("ALGO")
        setting.pop("TIME")

        Incar(setting).write_file(filename=incar_pre_gga)


class VaspDefectInputSetMaker(DefectInputSetMaker):

    def __init__(self, defect_initial_setting, filtering_words=None,
                 particular_defects=None, incar="INCAR", kpoints="KPOINTS",
                 force_overwrite=False):

        if not os.path.exists("INCAR") or not os.path.exists("KPOINTS"):
            raise VaspInputFileError("INCAR and/or KPOINTS is absent.")

        # make self._defect_initial_setting and self._defect_name_set
        super().__init__(defect_initial_setting, filtering_words,
                         particular_defects, force_overwrite)

        self._incar = incar
        self._kpoints = kpoints

        self.make_input()

    def _make_perfect_input(self):

        if self._force_overwrite:
            if os.path.exists("perfect"):
                print_is_being_removed("perfect")
                shutil.rmtree("perfect")

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

        # TODO: check if the defect_name is proper or not.
        if self._force_overwrite:
            if os.path.exists(defect_name):
                print_is_being_removed(defect_name)
                shutil.rmtree(defect_name)

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
            d.to_json_file(os.path.join(defect_name, "defect_entry.json"))
            d.initial_structure.to(
                filename=os.path.join(defect_name, "POSCAR-Initial"))

            if not self._defect_initial_setting.distance == 0.0:
                center = defect_center(d)
#                center = defect_center(d.initial_structure, d)
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


class VaspInputFileError(Exception):
    pass


