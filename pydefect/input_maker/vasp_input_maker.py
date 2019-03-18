# -*- coding: utf-8 -*-

from collections import OrderedDict
from copy import deepcopy
from enum import Enum, unique
from math import ceil
import numpy as np
import os
import warnings

import re
import ruamel.yaml as yaml

from monty.serialization import loadfn
from monty.io import zopen

from obadb.util.structure_handler import structure_to_seekpath, find_hpkot_primitive

from pymatgen.io.vasp import PotcarSingle, Potcar, Kpoints, Incar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.io_utils import clean_lines

from pydefect.core.prior_info import PriorInfo
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting, \
    element_set
from pydefect.database.kpt_centering import kpt_centering
from pydefect.database.atom import symbols_to_atom, u_parameter, \
    unoccupied_bands

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
POTCAR_LIST_FILE = os.path.join(MODULE_DIR, "default_POTCAR_list.yaml")

incar_flags = OrderedDict()
incar_flags["algo"] = ["ALGO"]
incar_flags["accuracy"] = ["PREC", "LREAL", "LASPH", "ENCUT", "NELM", "NELMIN", "EDIFF"]
incar_flags["symmetry"] = ["ISYM", "SYMPREC"]
incar_flags["ionic"] = ["ISIF", "EDIFFG", "IBRION", "NSW"]
incar_flags["occupation"] = ["ISMEAR", "SIGMA", "FERDO", "FERWE", "SMEARINGS"]
incar_flags["magnetization"] = ["ISPIN", "LNONCOLLINEAR", "MAGMOM", "LSORBIT", "GGA_COMPAT"]
incar_flags["mixing"] = ["MIX", "INIMIX", "MAXMIX", "AMIX", "BMIX", "AMIX_MAG", "BMIX_MAG", "AMIN", "MIXPRE", "WC"]
incar_flags["io_control"] = ["ISTART", "ICHARG", "LWAVE", "LCHARG"]
incar_flags["analysis"] = ["NBANDS", "LORBIT", "NEDOS", "EMAX", "EMIN", "RWIGS"]
incar_flags["dielectric"] = ["LEPSILON", "LRPA", "LOPTICS", "LCALCEPS", "EFIELD", "CSHIFT", "OMEGAMIN", "OMEGAMAX", "LNABLA"]
incar_flags["dft"] = ["GGA", "IVDW", "VDW_S6", "VDW_A1", "VDW_S8", "VDW_A2", "METAGGA"]
incar_flags["hybrid"] = ["LHFCALC", "PRECFOCK", "AEXX", "HFSCREEN", "NKRED", "TIME"]
incar_flags["gw"] = ["NBANDSGW", "NOMEGA", "OMEGAMAX", "ANTIRES", "MAXMEM", "NMAXFOCKAE", "NBANDSLF", "ENCUTGW"]
incar_flags["ldau"] = ["LDAU", "LDAUTYPE", "LMAXMIX", "LDAUPRINT", "LDAUL", "LDAUU", "LDAUJ"]
incar_flags["electrostatic"] = ["NELECT", "EPSILON", "DIPOL", "IDIPOL", "LMONO", "LDIPOL"]
incar_flags["parallel"] = ["NPAR", "NCORE", "KPAR"]


@unique
class Task(Enum):
    """
    Supported tasks in pydefect
    """
    structure_opt = "structure_opt"
    band = "band"
    dos = "dos"
    dielectric = "dielectric"
    dielectric_function = "dielectric_function"
    competing_phase = "competing_phase"
    competing_phase_molecule = "competing_phase_molecule"
    defect = "defect"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):
        for m in Task:
            if m.name == s:
                return m
        raise AttributeError("Task: " + str(s) + " is not proper.\n" +
                             "Supported Task:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])


@unique
class Functional(Enum):
    """
    Supported functionals in pydefect
    """
    pbe = "pbe"
    pbesol = "pbesol"
    pbe_d3 = "pbe_d3"
    scan = "scan"
    hse = "hse"
    nhse = "nhse"
    pbe0 = "pbe0"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):
        for m in Functional:
            if m.name == s:
                return m
        raise AttributeError("Functional: " + str(s) + " is not proper.\n" +
                             "Supported Functional:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])


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
                comment = False
                blank_line = False
                for v in val:
                    if v in self.keys():
                        if comment is False:
                            fw.write("# " + str(key) + "\n")
                            comment = True
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
        if k == "DEFAULT_POTCAR_DIR":
            potcar_director_path = v

    try:
        return potcar_director_path
    except ValueError:
        print('DEFAULT_POTCAR_DIR must be set in .pydefect.yaml')


def make_potcar(elements, potcar="POTCAR"):
    """
    Writes POTCAR with a list of the given element names.
    """

    list_file = loadfn(POTCAR_LIST_FILE)

    with open(potcar, 'w') as pot:
        for e in elements:
            potcar_file_name = \
                os.path.join(potcar_dir(), list_file[e], "POTCAR")

            with open(potcar_file_name) as each_pot:
                pot.write(each_pot.read())


def make_hpkot_primitive_poscar(poscar="POSCAR", pposcar="PPOSCAR",
                                symprec=1e-05, angle_tolerance=-1.0):
    s = Structure.from_file(os.path.join(poscar))
    primitive_structure = find_hpkot_primitive(s, symprec=symprec,
                                               angle_tolerance=angle_tolerance)
    primitive_structure.to(filename=pposcar)


def make_band_kpoints(ibzkpt, dirname='.', poscar="POSCAR",
                      num_split_kpoints=1, time_reversal=True):
    """
    Write the KPOINTS file for the band structure calculation.
    Args:
        ibzkpt (str): Name of theIBZKPT-type file, which is used to set the
                      weighted k-points.
        dirname (str): Director name at which INCAR is written.
        poscar (str): input POSCAR-type file name
        num_split_kpoints (int): number of KPOINTS files used for a band
                                 structure calculation.
        time_reversal:
    """

    s = Structure.from_file(os.path.join(poscar))

    seekpath_full_info = structure_to_seekpath(structure=s,
                                               time_reversal=time_reversal)

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


def make_kpoints(task, dirname='.', poscar="POSCAR", num_split_kpoints=1,
                 ibzkpt="IBZKPT", is_metal=False, kpts_shift=None,
                 kpts_density_opt=3, kpts_density_defect=1.5, factor_dos=2,
                 factor_metal=2, prior_info=None, is_magnetization=False):
    """
    Constructs a KPOINTS file based on default settings depending on the task.
    Args:
        task (Task Enum):
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
        factor_dos (int): Multiplier factor for dos and dielectric
            constants.
            This is integer to make k-points even numbers. Then, we can use
            NKRED = 2 for hybrid functional calculations.
        factor_metal (int): Multiplier factor for metallic systems
            This is integer to make k-points even numbers. Then, we can use
            NKRED = 2 for hybrid functional calculations.
        prior_info (str): mp.json file name
        is_magnetization (bool): If the magnetization is considered or not.
    """

    s = Structure.from_file(os.path.join(poscar))
    reciprocal_lattice = s.lattice.reciprocal_lattice

    task = Task.from_string(task)
    comment = "Generated by pydefect. task: " + str(task) + ", "

    if prior_info:
        pi = PriorInfo.load_json(filename=prior_info)
        if pi.is_metal:
            is_metal = True
        if pi.is_molecule:
            task = Task.competing_phase_molecule

    # band
    if task == Task.band:
        if is_magnetization:
            time_reversal = False
        else:
            time_reversal = True

        make_band_kpoints(ibzkpt, dirname, poscar, num_split_kpoints,
                          time_reversal)
        return True
    # molecule: task is competing_phase_molecule or is_molecule = .True.
    elif task == Task.competing_phase_molecule:
        kpt_mesh = [1, 1, 1]
        comment += "Gamma-only for a cluster."
    # defect
    elif task == Task.defect:
        kpt_mesh = [int(np.ceil(kpts_density_defect * r))
                    for r in reciprocal_lattice.abc]
        comment += "kpt density: " + str(kpts_density_defect)
    else:

        # For primitive cell calculations
        sg = SpacegroupAnalyzer(structure=s).get_space_group_number()

        # Note that the numbers of k-points along all the directions must be the
        # same for body centered orthorhombic and body centered tetragonal to
        # keep the symmetry.
        body_centered_orthorhombic = (23, 24, 44, 45, 46, 71, 72, 73, 74)
        body_centered_tetragonal = (79, 80, 81, 82, 87, 88, 97, 98, 107, 108,
                                    109, 110, 119, 120, 121, 122, 139, 140, 141,
                                    142)
        if sg in body_centered_orthorhombic or body_centered_tetragonal:
            average_abc = np.prod(reciprocal_lattice.abc) ** (1.0 / 3.0)
            revised_reciprocal_abc = (average_abc, average_abc, average_abc)
        else:
            revised_reciprocal_abc = reciprocal_lattice.abc

        # insulating phase
        if task in (Task.competing_phase, Task.structure_opt) \
                and is_metal is False:
            kpt_mesh = [int(np.ceil(kpts_density_opt * r))
                        for r in revised_reciprocal_abc]
            comment += "kpt density: " + str(kpts_density_defect)
        # metallic phase
        elif task == Task.competing_phase and is_metal is True:
            kpt_mesh = [int(np.ceil(kpts_density_opt * r)) * factor_metal
                        for r in revised_reciprocal_abc]
            comment += "kpt density: " + str(kpts_density_defect) + ", " + \
                       "factor: " + str(factor_metal)
        # dos, dielectric const, dielectric function
        elif task in (Task.dos, Task.dielectric, Task.dielectric_function):
            kpt_mesh = [int(np.ceil(kpts_density_opt * r)) * factor_dos
                        for r in revised_reciprocal_abc]
            comment += "kpt density: " + str(kpts_density_defect) + ", " + \
                       "factor: " + str(factor_dos)
        else:
            raise AttributeError("Task: " + str(task) + " is not supported.\n")

    kpts = (tuple(kpt_mesh),)

    if kpts_shift is None:
        if task in (Task.structure_opt, Task.defect):
            angles = s.lattice.angles

            # if not kpts_shift:
            kpts_shift = []
            for i in range(3):
                # shift kpt mesh center only for the lattice vector being normal
                # to a lattice plane and even number of k-points.
                if angles[i - 2] == 90 and angles[i - 1] == 90 \
                        and kpts[0][i] % 2 == 0:
                    kpts_shift.append(0.5)
                else:
                    kpts_shift.append(0.0)

        elif task in (Task.competing_phase_molecule, Task.dos,
                      Task.dielectric, Task.dielectric_function):
            kpts_shift = [0, 0, 0]

        elif task is Task.competing_phase:
            # Check the c-axis for hexagonal phases.
            if 168 <= sg <= 194:
                angles = s.lattice.angles
                if angles[0] != 90 or angles[1] != 90:
                    raise ValueError("Axial in the hexagonal structure is not "
                                     "set properly.")
            kpts_shift = kpt_centering[sg]
            for i in range(3):
                if kpts[0][i] % 2 == 1:
                    kpts_shift[i] = 0

    Kpoints(comment=comment, kpts=kpts, kpts_shift=kpts_shift). \
        write_file(os.path.join(dirname, 'KPOINTS'))


class MakeIncar:

    def __init__(self, task, functional, poscar="POSCAR", potcar=None,
                 hfscreen=None, aexx=None, is_metal=False,
                 is_magnetization=False, ldau=False, dirname='.',
                 defect_in=None, prior_info=None, my_incar_setting=None):
        """
        Constructs an INCAR file based on default settings depending on the task
        ENCUT can be determined from max(ENMAX).
        Args:
            task (Task Enum):
            functional (Functional Enum):
            poscar (str): POSCAR-type file name parsed for checking elements
            potcar (str): POTCAR-type file
            hfscreen (float): VASP parameter of screening distance HFSCREEN
            aexx (float): VASP parameter of the Fock exchange AEXX
            is_metal (bool): If the system to be calculated is metal or not.
                             This has a meaning for the calculation of a
                             competing phase
            is_magnetization (bool): If the magnetization is considered or not.
            ldau (bool): If the on-site Coulomb potential is considered.
            dirname (str): Director name at which INCAR is written.
            defect_in (str): defect.in file name parsed for checking elements
            prior_info (str): prior_info.json file name
            my_incar_setting (str): user setting yaml file name
        """

        self.task = Task.from_string(task)
        self.functional = Functional.from_string(functional)
        self.potcar = potcar
        self.setting_file = loadfn(INCAR_SETTINGS_FILE)
        self.defect_in = defect_in

        # raise error when incompatible combinations are given.
        if self.functional == Functional.nhse and \
                self.task not in (Task.band, Task.dos):
            raise IncarIncompatibleError(
                "nHSE is not compatible with Task: {}.".format(str(self.task)))

        self.setting = {}
        for i in (self.task, self.functional):
            try:
                s = self.setting_file[i.name]
                if s:
                    self.setting.update(s)
            except IOError:
                print('default_INCAR_self.setting.yaml cannot be opened')

        # Parse POSCAR file.
        try:
            self.structure = Structure.from_file(poscar)
        except:
            raise IOError("POSCAR type file is imperative.")
        self.elements = self.structure.symbol_set

        max_enmax = self._get_max_enmax()

        # In case the lattice shape is relaxed, ENCUT is multiplied by 1.3 to
        # accurately calculate stress.
        if "ISIF" in self.setting and self.setting["ISIF"] > 2:
            self.setting["ENCUT"] = str(max_enmax * 1.3)
        else:
            self.setting["ENCUT"] = str(max_enmax)

        if self.task in (Task.band, Task.dos, Task.dielectric_function):
            self._add_nbands_flag()

        if self.task in (Task.dos, Task.dielectric_function):
            self._add_emin_emax_nedos()

        if self.functional in (Functional.hse, Functional.nhse):
            if hfscreen:
                self.setting["HFSCREEN"] = hfscreen
            if aexx:
                self.setting["AEXX"] = aexx

        # read prior_info here.
        if prior_info:
            pi = PriorInfo.load_json(filename=prior_info)
            # check magnetization
            if pi.is_magnetic:
                self.setting["ISPIN"] = 2
            if pi.is_metal:
                is_metal = pi.is_metal

        if is_magnetization:
            self.setting["ISPIN"] = 2

        if my_incar_setting:
            self.setting.update(loadfn(my_incar_setting))

        # From here, we deal with particular issues related to INCAR
        if ldau:
            self._add_plus_u_flags()

        # remove NPAR for the calculations of dielectrics as it is not allowed.
        if self.task in (Task.dielectric, Task.dielectric_function):
            if "NPAR" in self.setting:
                self.setting.pop("NPAR")

        # remove KPAR for band structure calculations
        if self.task in (Task.band, Task.competing_phase_molecule):
            if "KPAR" in self.setting:
                self.setting.pop("KPAR")

        # add NKRED for dos, dielectric function or metal.
        if self.functional in (Functional.hse, Functional.nhse) \
                and (self.task in (Task.dos, Task.dielectric_function)
                     or is_metal):
            # the following number should be the same as factor_metal in
            # make_kpoints
            self.setting["NKRED"] = 2

        # split the function to write_file
        incar = os.path.join(dirname, 'INCAR')
        ModIncar(self.setting).pretty_write_file(filename=incar)

        if functional in (Functional.hse, Functional.nhse):
            incar_pre_gga = os.path.join(dirname, 'INCAR-pre_gga')

            # pop out hybrid related flags.
            for i in incar_flags["hybrid"]:
                if i in self.setting:
                    self.setting.pop(i)

            # insert some flags.
            self.setting.update(self.setting_file["hse_pre"])
            ModIncar(self.setting).pretty_write_file(filename=incar_pre_gga)

    def _get_max_enmax(self, ):

        # ENCUT is determined by one of three ways.
        # 1. READ my_incar_self.setting.yaml file. (overwrite later)
        # 2. Parse POTCAR file.
        # 3. Parse defect.in file, check relevant elements including foreign
        #    atoms and determine the max ENMAX value.
        # 4. Parse POSCAR file, and check the relevant elements. Note that
        #    dopants are not considered in this case.
        if self.potcar:
            potcar = Potcar.from_file(self.potcar)
            enmax = [p.enmax for p in potcar]
        else:
            if self.defect_in:
                defect_initial_setting = \
                    DefectInitialSetting.from_defect_in(self.poscar,
                                                        self.defect_in)
                name_set = element_set(defect_initial_setting)
            else:
                name_set = self.elements
            list_file = loadfn(POTCAR_LIST_FILE)
            enmax = []
            for e in name_set:
                potcar_file_name = \
                    os.path.join(potcar_dir(), list_file[e], "POTCAR")
            enmax.append(PotcarSingle.from_file(potcar_file_name).enmax)

        return max(enmax)

    def _add_nbands_flag(self):

        composition = self.structure.composition

        num_electrons = {}
        if self.potcar:
            potcar = Potcar.from_file(self.potcar)
            for p in potcar:
                num_electrons[p.element] = p.nelectrons
        else:
            list_file = loadfn(POTCAR_LIST_FILE)
            for e in self.elements:
                potcar = os.path.join(potcar_dir(), list_file[e], "POTCAR")
                num_electrons[e] = PotcarSingle.from_file(potcar).nelectrons

        if len(num_electrons) != len(composition):
            raise ValueError("The number of elements is not consistent.")

        num_bands = 0.0
        for e in self.elements:
            num_bands += \
                composition[e] * (num_electrons[e] / 2 + unoccupied_bands[e])

        self.setting["NBANDS"] = ceil(num_bands)

    def _add_emin_emax_nedos(self):
        self.setting["EMIN"] = -20
        self.setting["EMAX"] = 10
        self.setting["NEDOS"] = 3001

    def _add_plus_u_flags(self):
        # Only LDAUTYPE = 2 is supported
        ldaul = []
        ldauu = []
        ldauj = []
        for e in self.elements:
            if e in u_parameter.keys():
                a = u_parameter[e]
                ldaul.append(a[0])
                ldauu.append(a[1])
            else:
                ldaul.append(-1)
                ldauu.append(0)
            # LDAUJ is always 0.
            ldauj.append(0)

        if max(ldaul) == 3:
            self.setting["LMAXMIX"] = 6
        elif max(ldaul) == 2:
            self.setting["LMAXMIX"] = 4
        elif max(ldaul) == -1:
            return False

        # These setting must be after checking the LMAXMIX.
        self.setting["LDAU"] = True
        self.setting["LDAUTYPE"] = 2
        self.setting["LDAUPRINT"] = 2

        self.setting["LDAUL"] = \
            " ".join(['{0:3d}'.format(l, end="") for l in ldaul])
        self.setting["LDAUU"] = \
            " ".join(['{0:3d}'.format(u, end="") for u in ldauu])
        self.setting["LDAUJ"] = \
            " ".join(['{0:3.1f}'.format(j, end="") for j in ldauj])


class IncarIncompatibleError(Exception):
    pass
