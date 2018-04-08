# -*- coding: utf-8 -*-

from abc import ABCMeta
from itertools import product
import json
import numpy as np
import warnings

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.electronic_structure.core import Spin

from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def min_distance_and_its_v2coord(v1, v2, lattice_vector_matrix):
    """
    """

    candidate = []
    v = np.dot(lattice_vector_matrix, v2 - v1)

    for index in product((-1, 0, 1), repeat=3):
        index = np.array(index)
        delta_v = np.dot(lattice_vector_matrix, index)
        distance = np.linalg.norm(delta_v + v)
        candidate.append(distance)

    return min(candidate)


def distance_list(structure, defect_coords):
    """
    """

    lattice_vector_matrix = structure.lattice.matrix

    return [min_distance_and_its_v2coord(host_atom_coords,
                                         defect_coords,
                                         lattice_vector_matrix)
            for host_atom_coords in structure.frac_coords]


class DftResults(metaclass=ABCMeta):
    """
    Abstract class holding some DFT results used for defect analysis.
    Args:
        final_structure (Structure): Optimized Structure
        total_energy (float): total energy
        eigenvalues (N_k x N_band numpy array): eigenvalues
        electrostatic_potential (list): electrostatic potential at atomic sites
    """
    def __init__(self, final_structure, total_energy, eigenvalues,
                 electrostatic_potential):
        self._final_structure = final_structure
        self._total_energy = total_energy
        self._eigenvalues = eigenvalues
        self._electrostatic_potential = electrostatic_potential

    @classmethod
    def from_vasp_files(cls, directory_path, contcar_name="CONTCAR",
                        outcar_name="OUTCAR", vasprun_name="vasprun.xml"):
        """
        Although electrostatic_potential is not used for UnitcellDftResults,
        this method is implemented in DftResults class because constructor is
        easily written.

        Args:
            directory_path (str): path of directory.
            contcar_name (str): Name of converged CONTCAR file.
                                Defaults to CONTCAR.
            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
            vasprun_name (str): Name of vasprun.xml file.
                                Defaults to vasprun.xml.
        """
        contcar = Poscar.from_file(directory_path + "/" + contcar_name)
        outcar = Outcar(directory_path + "/" + outcar_name)
        vasprun = Vasprun(directory_path + "/" + vasprun_name)

        final_structure = contcar.structure
        total_energy = outcar.final_energy
        eigenvalues = vasprun.eigenvalues
        electrostatic_potential = outcar.electrostatic_potential

        return cls(final_structure, total_energy, eigenvalues,
                   electrostatic_potential)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        return cls(d["final_structure"], d["total_energy"], eigenvalues,
                   d["electrostatic_potential"])

    @classmethod
    def json_load(cls, filename):
        """
        Constructs a class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    def as_dict(self):
        """
        Dict representation of DefectInitialSetting class object.
        Json-serializable dict representation.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self._eigenvalues.items()}

        d = {"final_structure":         self._final_structure,
             "total_energy":            self._total_energy,
             "eigenvalues":             eigenvalues,
             "electrostatic_potential": self._electrostatic_potential}

        return d

    def to_json_file(self, filename):
        """
        Returns a json file.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @property
    def final_structure(self):
        return self._final_structure

    @property
    def eigenvalues(self):
        return self._eigenvalues

    @property
    def total_energy(self):
        return self._total_energy

    @property
    def electrostatic_potential(self):
        return self._electrostatic_potential


class SupercellDftResults(DftResults):
    """
    Subclass holding DFT results for supercell systems both w/ and w/o a defect.
    """
    def relative_total_energy(self, perfect_dft_results):
        """
        Return relative total energy w.r.t. the perfect supercell.

        Args:
            perfect_dft_results (SupercellDftResults):
                SupercellDftResults class object for the perfect supercell.
        """
        return self._total_energy - perfect_dft_results.total_energy

    def relative_potential(self, perfect_dft_results, defect_entry):
        """
        Return a list of relative site potential w.r.t. the perfect supercell.
        Note that None is inserted for interstitial sites.

        Args:
            perfect_dft_results (SupercellDftResults):
                SupercellDftResults class object for the perfect supercell.
            defect_entry (DefectEntry):
                DefectEntry class object.
        """
        mapping = defect_entry.atom_mapping_to_perfect

        relative_potential = []

        for d_atom, p_atom in enumerate(mapping):

            if p_atom is None:
                relative_potential.append(None)
            else:
                ep_defect = self.electrostatic_potential[d_atom]
                ep_perfect = perfect_dft_results.electrostatic_potential[p_atom]
                relative_potential.append(ep_defect - ep_perfect)

        return relative_potential

#    def inserted_atom_displacements(self, defect_entry):
#        """
#        Returns coordinates of defect center by calculating the averaged
#        coordinates. If len(defect_coords) == 1, returns defect_coords[0].
#        Args:
#            defect_entry (DefectEntry):
#                related DefectEntry class object
#        """
#        displacements = []
#
#        for k in defect_entry.inserted_atoms.keys:
#            before_relaxation = defect_entry.initial_structure.frac_coords[k]
#            after_relaxation = self.final_structure.frac_coords[k]
#            displacements.append(
#                min_distance_and_its_v2coord(before_relaxation,
#                                             after_relaxation,
#                                             self.final_structure.axis))
#        return displacements

    def defect_center(self, defect_entry):
        """
        Returns a fractional coordinates of the defect center which is
        calculated by averaging the coordinates of vacancies and interstitials.
        If len(defect_coords) == 1, returns defect_coords[0].

        Args:
            defect_entry (DefectEntry): related DefectEntry class object
        """
        inserted_atom_coords = \
            list([self.final_structure.frac_coords[k]
                  for k in defect_entry.inserted_atoms.keys()])
        removed_atom_coords = list(defect_entry.removed_atoms.values())

        defect_coords = inserted_atom_coords + removed_atom_coords

        # np.array([[0, 0.1, 0.2], [0.3, 0.4, 0.5]]).transpose() =
        # np.array([[0, 0.3], [0.1, 0.4], [0.2, 0.5]])
        return [np.mean(i) for i in np.array(defect_coords).transpose()]

    def distances_from_a_point(self, defect_entry):
        """
        Returns a list of distances at atomic sites from a defect center defined
        by defect_entry. Note that in the case of an interstitial-type defect,
        zero is also set.

        Args:
            defect_entry (DefectEntry): related DefectEntry class object
        """
        return distance_list(self.final_structure,
                             self.defect_center(defect_entry))


class UnitcellDftResults(DftResults):
    """
    DFT results for a unitcell
    Args:
        static_dielectric_tensor (3x3 numpy array):
        ionic_dielectric_tensor (3x3 numpy array):
    """

    def __init__(self, final_structure, total_energy, eigenvalues,
                 electrostatic_potential, static_dielectric_tensor=None,
                 ionic_dielectric_tensor=None):
        """ """
        super().__init__(final_structure, total_energy, eigenvalues,
                         electrostatic_potential)

        self._static_dielectric_tensor = static_dielectric_tensor
        self._ionic_dielectric_tensor = ionic_dielectric_tensor

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        return cls(d["final_structure"], d["total_energy"], eigenvalues,
                   d["electrostatic_potential"], d["static_dielectric_tensor"],
                   d["ionic_dielectric_tensor"])

    def set_dielectric_constants_from_outcar(self, directory_path,
                                             outcar_name="OUTCAR"):
        outcar = Outcar(directory_path + "/" + outcar_name)
        self._static_dielectric_tensor = np.array(outcar.dielectric_tensor)
        self._ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)

    @property
    def static_dielectric_tensor(self):
        if self._static_dielectric_tensor is not None:
            return self._static_dielectric_tensor
        else:
            warnings.warn(message="Static dielectric tensor is not set yet.")
            return None

    @property
    def ionic_dielectric_tensor(self):
        if self._ionic_dielectric_tensor is not None:
            return self._ionic_dielectric_tensor
        else:
            warnings.warn(message="Ionic dielectric tensor is not set yet.")
            return None

    @property
    def total_dielectric_tensor(self):
        return self._static_dielectric_tensor + self._ionic_dielectric_tensor

    @static_dielectric_tensor.setter
    def static_dielectric_tensor(self, static_dielectric_tensor):
        self._static_dielectric_tensor = static_dielectric_tensor

    @ionic_dielectric_tensor.setter
    def ionic_dielectric_tensor(self, ionic_dielectric_tensor):
        self._ionic_dielectric_tensor = ionic_dielectric_tensor

    def set_static_dielectric_tensor_from_outcar(self, directory_path,
                                                 outcar_name="OUTCAR"):
        outcar = Outcar(directory_path + "/" + outcar_name)
        self._static_dielectric_tensor = np.array(outcar.dielectric_tensor)

    def set_ionic_dielectric_tensor_from_outcar(self, directory_path,
                                                outcar_name="OUTCAR"):
        outcar = Outcar(directory_path + "/" + outcar_name)
        self._ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)

    def set_dielectric_tensor_from_outcar(self, directory_path,
                                          outcar_name="OUTCAR"):
        self.set_static_dielectric_tensor_from_outcar(directory_path,
                                                      outcar_name)
        self.set_ionic_dielectric_tensor_from_outcar(directory_path,
                                                     outcar_name)

    def as_dict(self):
        """
        Dict representation of DefectInitialSetting class object.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self._eigenvalues.items()}

        d = {"final_structure":          self._final_structure,
             "total_energy":             self._total_energy,
             "eigenvalues":              eigenvalues,
             "electrostatic_potential":  self._electrostatic_potential,
             "static_dielectric_tensor": self._static_dielectric_tensor,
             "ionic_dielectric_tensor":  self._ionic_dielectric_tensor}

        return d
