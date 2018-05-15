# -*- coding: utf-8 -*-

from collections import OrderedDict
from copy import deepcopy
from enum import Enum, auto, unique
from math import ceil
import numpy as np
import os

import re
import ruamel.yaml as yaml

from monty.serialization import loadfn
from monty.io import zopen

from pymatgen.io.vasp import Potcar, Kpoints, Incar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.io_utils import clean_lines

from pydefect.input_maker.defect_initial_setting import DefectInitialSetting, \
    element_set
from pydefect.database.kpt_centering import kpt_centering
from pydefect.database.atom import symbols_to_atom
from pydefect.util.structure import structure_to_seekpath, find_hpkot_primitive

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

incar_flags = OrderedDict()
incar_flags["algo"] = ["ALGO"]
incar_flags["accuracy"] = ["PREC", "LREAL", "LASPH", "ENCUT", "NELM", "NELMIN", "EDIFF"]
incar_flags["symmetry"] = ["ISYM", "SYMPREC"]
incar_flags["ionic"] = ["ISIF", "EDIFFG", "IBRION", "NSW"]
incar_flags["occupation"] = ["ISMEAR", "SIGMA", "FERDO", "FERWE", "SMEARINGS"]
incar_flags["magnetization"] = ["LNONCOLLINEAR", "MAGMOM", "LSORBIT", "GGA_COMPAT"]
incar_flags["mixing"] = ["MIX", "INIMIX", "MAXMIX", "AMIX", "BMIX", "AMIX_MAG", "BMIX_MAG", "AMIN", "MIXPRE", "WC"]
incar_flags["io_control"] = ["ISTART", "ICHARG", "LWAVE", "LCHARG"]
incar_flags["analysis"] = ["NBANDS", "LORBIT", "NEDOS", "EMAX", "EMIN", "RWIGS"]
incar_flags["dielectric"] = ["LEPSILON", "LRPA", "LOPTICS", "LCALCEPS", "EFIELD", "CSHIFT", "OMEGAMIN", "OMEGAMAX", "LNABLA"]
incar_flags["dft"] = ["GGA", "IVDW", "VDW_S6", "VDW_A1", "VDW_S8", "VDW_A2", "METAGGA"]
incar_flags["hybrid"] = ["LHFCALC", "AEXX", "HFSCREEN", "NKRED", "PRECFOCK", "TIME"]
incar_flags["gw"] = ["NBANDSGW", "NOMEGA", "OMEGAMAX", "ANTIRES", "MAXMEM", "NMAXFOCKAE", "NBANDSLF", "ENCUTGW"]
incar_flags["ldau"] = ["LDAU", "LDAUTYPE", "LMAXMIX", "LDAUPRINT", "LDAUL", "LDAUU", "LDAUJ"]
incar_flags["electrostatic"] = ["NELECT", "EPSILON", "DIPOL", "IDIPOL", "LMONO", "LDIPOL"]
incar_flags["parallel"] = ["NPAR", "KPAR"]


class ModIncar(Incar):
    """
    Incar class modified for pretty writing of INCAR file.
    Since from_file and from_string methods in Incar class use Incar class
    constructor, we need to override them.
    """
    @staticmethod
    def from_file(filename):
        """
        Reads an Incar object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            ModIncar object
        """
        with zopen(filename, "rt") as f:
            return ModIncar.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an Incar object from a string.

        Args:
            string (str): Incar string

        Returns:
            ModIncar object
        """
        lines = list(clean_lines(string.splitlines()))
        params = {}
        for line in lines:
            m = re.match(r'(\w+)\s*=\s*(.*)', line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip()
                val = Incar.proc_val(key, val)
                params[key] = val
        return ModIncar(params)

    def pretty_write_file(self, filename):
        """
        Write Incar to a file in a prettier manner.

        Args:
            filename (str): filename to write to.
        """

        with zopen(filename, "wt") as fw:
            for key, val in incar_flags.items():
                blank_line = False
                for v in val:
                    if v in self.keys():
                        fw.write(v + " = " + str(self[v]) + "\n")
                        blank_line = True
                if blank_line:
                    fw.write("\n")
        return None


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


def make_hpkot_primitive_poscar(poscar="POSCAR", pposcar="PPOSCAR"):
    s = Structure.from_file(os.path.join(poscar))
    primitive_structure = find_hpkot_primitive(s)
    primitive_structure.to(filename=pposcar)


def make_supercell_poscar(scaling_matrix, poscar="POSCAR", sposcar="SPOSCAR"):
    s = Structure.from_file(os.path.join(poscar))
    supercell_structure = s.make_supercell(scaling_matrix)
    supercell_structure.to(filename=sposcar)


def make_band_kpoints(ibzkpt, dirname='.', poscar="POSCAR",
                      num_split_kpoints=1):
    """
    Write the KPOINTS file for the band structure calculation.
    Args:
        ibzkpt (str): Name of theIBZKPT-type file, which is used to set the
                      weighted k-points.
        dirname (str): Director name at which INCAR is written.
        poscar (str): input POSCAR-type file name
        num_split_kpoints (int): number of KPOINTS files used for a band
                                 structure calculation.
    """

    s = Structure.from_file(os.path.join(poscar))

    seekpath_full_info = structure_to_seekpath(s)

    # primitive structure
    lattice = seekpath_full_info["primitive_lattice"]
    element_types = seekpath_full_info["primitive_types"]
    species = [symbols_to_atom[i] for i in element_types]
    positions = seekpath_full_info["primitive_positions"]
    primitive = Structure(lattice, species, positions)

    # It would be great if sg is obtained from seekpath.
    # Note: Parameters used for the symmetry search are different between
    #       seekpath and pymatgen.
    sg = SpacegroupAnalyzer(structure=primitive).get_space_group_number()

    # k-path
    kpath = seekpath_full_info["explicit_kpoints_rel"]
    kpath_label = seekpath_full_info["explicit_kpoints_labels"]
    ibzkpt = Kpoints.from_file(ibzkpt)

    num_kpoints = ceil(len(kpath) / num_split_kpoints)
    for x in range(num_split_kpoints):
        kpoints = deepcopy(ibzkpt)

        kpoints.comment = \
            "Generated by pydefect. Formula: " + \
            primitive.composition.reduced_formula + " SG: " + str(sg)

        divided_kpath = kpath[num_kpoints * x: num_kpoints * (x + 1)]
        divided_kpath_label = \
            kpath_label[num_kpoints * x: num_kpoints * (x + 1)]

        for d, dl in zip(divided_kpath, divided_kpath_label):
            kpoints.num_kpts += 1
            kpoints.kpts.append(d)
            # weight zero is set here.
            kpoints.kpts_weights.append(0)
            kpoints.labels.append(dl)

        if num_split_kpoints > 1:
            output_filename = os.path.join(dirname, "KPOINTS-" + str(x))
        else:
            output_filename = os.path.join(dirname, "KPOINTS")

        kpoints.write_file(output_filename)


@unique
class SupportedTask(Enum):
    """
    Supported tasks in pydefect
    """
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
            if m.name == s:
                return m
        raise ValueError("Can't interpret supported task %s" % s)


@unique
class SupportedFunctional(Enum):
    """
    Supported functionals in pydefect
    """
    pbe = auto()
    hse06 = auto()
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


def make_kpoints(task, dirname='.', poscar="POSCAR", num_split_kpoints=1,
                 ibzkpt="IBZKPT", is_metal=False, kpts_shift=None,
                 kpts_density_opt=3, kpts_density_defect=1.5,
                 multiplier_factor=2, multiplier_factor_metal=2):
    """
    Constructs a KPOINTS file based on default settings depending on the task.
    Args:
        task (SupportedTask Enum):
        dirname (str): Director name at which INCAR is written.
        poscar (str): input POSCAR-type file name
        num_split_kpoints (int): number of KPOINTS files used for a band
                                 structure calculation.
        ibzkpt (str): Name of the IBZKPT-type file, which is used to set the
                      weighted k-points.
        is_metal (bool): If the system to be calculated is metal or not. This
                         has a meaning for the calculation of a competing phase
        kpts_shift (1x3 list): Shift of k-point set based on vasp setting
        kpts_density_opt (float): Density of k-point = (N_k / k_length) for
                                  structure optimization
        kpts_density_defect (float): Density of k-point = (N_k / k_length) for
                                     point-defect calculations
        multiplier_factor (float): Multiplier factor for dos and dielectric
                                   constants
        multiplier_factor_metal (float): Multiplier factor for metallic systems
    """

    s = Structure.from_file(os.path.join(poscar))
    reciprocal_lattice = s.lattice.reciprocal_lattice

    task = SupportedTask.from_string(task)
    # band
    if task == SupportedTask.band:
        make_band_kpoints(ibzkpt, dirname, poscar, num_split_kpoints)
        return True
    # defect
    elif task == SupportedTask.defect:
        kpt_mesh = [int(np.ceil(kpts_density_defect * r))
                    for r in reciprocal_lattice.abc]
    else:
        # metallic competing phase
        if task == SupportedTask.competing_phase and is_metal is True:
            kpt_mesh = \
                [int(np.ceil(kpts_density_opt * r * multiplier_factor_metal))
                 for r in reciprocal_lattice.abc]
        # dos, dielectric const, dielectric function
        elif task in (SupportedTask.dos,
                      SupportedTask.dielectric,
                      SupportedTask.dielectric_function):
            kpt_mesh = \
                [int(np.ceil(kpts_density_opt * r) * multiplier_factor)
                 for r in reciprocal_lattice.abc]
        # molecule
        elif task == SupportedTask.competing_phase_molecule:
            kpt_mesh = [1, 1, 1]
        # insulating competing phase
        else:
            kpt_mesh = [int(np.ceil(kpts_density_opt * r))
                        for r in reciprocal_lattice.abc]

    if kpts_shift is None:
        if task in (SupportedTask.structure_opt, SupportedTask.defect):
            angles = s.lattice.angles

            if not kpts_shift:
                kpts_shift = []
                for i in range(3):
                    if angles[i - 2] == 90 and angles[i - 1] == 90:
                        kpts_shift.append(0.5)
                    else:
                        kpts_shift.append(0.0)

        elif task == SupportedTask.competing_phase_molecule:
            kpts_shift = [0, 0, 0]

        else:
            # For primitive cell calculations
            sg = SpacegroupAnalyzer(structure=s).get_space_group_number()
            # Check the c-axis for hexagonal phases.
            if 168 <= sg <= 194:
                angles = s.lattice.angles
                if angles[0] != 90 or angles[1] != 90:
                    raise ValueError("Axial in the hexagonal structure is not "
                                     "set properly.")
            kpts_shift = kpt_centering[sg]

    kpts = (tuple(kpt_mesh),)

    comment = "K-point mesh generated by pydefect."
    Kpoints(comment=comment, kpts=kpts, kpts_shift=kpts_shift). \
        write_file(os.path.join(dirname, 'KPOINTS'))


def get_max_enmax(element_names, dirname=potcar_dir()):
    enmax = []
    for e in element_names:
        potcar_file_name = os.path.join(dirname, "POTCAR_" + e)
        enmax.append(Potcar.from_file(potcar_file_name)[0].enmax)

    return max(enmax)


def make_incar(task, functional, hfscreen=None, aexx=None,
               is_magnetization=False, dirname='.', defect_in=None,
               poscar="POSCAR", my_incar_setting=None):
    """
    Constructs an INCAR file based on default settings depending on the task.
    ENCUT can be determined from max(ENMAX).
    Args:
        task (SupportedTask Enum):
        functional (SupportedFunctional Enum):
        hfscreen (float): VASP parameter of screening distance HFSCREEN
        aexx (float): VASP parameter of the Fock exchange AEXX
        is_magnetization (bool):
        dirname (str): Director name at which INCAR is written.
        defect_in (str): defect.in file name parsed for checking elements
        poscar (str): POSCAR-type file name parsed for checking elements
        my_incar_setting (str): user setting yaml file name

    defect.in file is parsed when exists, otherwise DPOSCAR file is.
    """

    setting_file = loadfn(INCAR_SETTINGS_FILE)

    try:
        task = SupportedTask.from_string(task)
    except AttributeError:
        print("SupportedTask:", task, "is not a proper flag.")

    try:
        functional = SupportedFunctional.from_string(functional)
    except AttributeError:
        print("SupportedFunctional:", functional, "is not a proper flag.")

    setting = {}
    for i in (task, functional):
        try:
            setting.update(setting_file[i.name])
        except IOError:
            print('default_INCAR_setting.yaml cannot be opened')

    if is_magnetization:
        setting["ISPIN"] = 2

    if defect_in:
        defect_initial_setting = \
            DefectInitialSetting.from_defect_in(poscar, defect_in)
        setting["ENCUT"] = \
            str(get_max_enmax(element_set(defect_initial_setting)))
    else:
        setting["ENCUT"] = input("Input ENCUT:")

    if functional is SupportedFunctional.hse06:
        if hfscreen:
            setting["HFSCREEN"] = hfscreen
        if aexx:
            setting["AEXX"] = aexx

    if my_incar_setting:
        setting.update(loadfn(my_incar_setting))

    incar = os.path.join(dirname, 'INCAR')
    ModIncar(setting).pretty_write_file(filename=incar)

    if functional is SupportedFunctional.hse06:
        incar_pre_gga = os.path.join(dirname, 'INCAR-pre_gga')

        # pop out hybrid related flags.
        for i in incar_flags["hybrid"]:
            if i in setting.keys():
                setting.pop(i)

        ModIncar(setting).pretty_write_file(filename=incar_pre_gga)

