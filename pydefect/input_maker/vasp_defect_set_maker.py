# -*- coding: utf-8 -*-

import os
import shutil

from pydefect.core.defect_entry import get_num_atoms_for_elements
from pydefect.util import get_num_electrons_from_potcar
from pydefect.core.supercell_dft_results import defect_center
from pydefect.input_maker.defect_set_maker import DefectMaker, \
    DefectInputSetMaker, print_is_being_removed, print_already_exist, \
    print_is_being_constructed
from pydefect.input_maker.vasp_input_maker import potcar_dir, make_potcar
from pydefect.util.structure import perturb_neighbors

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


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
                perturbed_defect_structure, perturbed_sites = \
                    perturb_neighbors(d.initial_structure,
                                      center,
                                      self._defect_initial_setting.cutoff,
                                      self._defect_initial_setting.distance)
                perturbed_defect_structure. \
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
