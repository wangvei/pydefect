# -*- coding: utf-8 -*-

from collections import defaultdict
from itertools import product
import json
import numpy as np
import os
import warnings

from monty.json import MontyEncoder
from monty.serialization import loadfn

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.electronic_structure.core import Spin

from pydefect.core.defect_entry import DefectEntry
from pydefect.vasp_util.script.vasp_process_analyzer import VaspNotConvergedError


__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"


def defect_center_from_coords(inserted_atom_coords, removed_atom_coords, structure):
    """
    """
    defect_coords = inserted_atom_coords + removed_atom_coords
    # First defect_coords is used as a base point under periodic boundary
    # condition. Here, we are aware of the case when two defect positions
    # are, e.g., [0.01, 0.01, 0.01] and [0.99, 0.99, 0.99].
    base = defect_coords[0]
    shortest_defect_coords = []

    for dc in defect_coords:
        diff = min_distance_under_pbc(
            np.array(base), np.array(dc), structure.lattice.matrix)[1]
        shortest_defect_coords.append(dc + diff)

    # np.array([[0, 0.1, 0.2], [0.3, 0.4, 0.5]]).transpose() =
    # np.array([[0, 0.3], [0.1, 0.4], [0.2, 0.5]])
    return [np.mean(i) for i in np.array(shortest_defect_coords).transpose()]


def defect_center(defect_entry, structure=None):
    """
    Return a fractional coordinates of the defect center that is calculated
    by averaging the coordinates of vacancies and interstitials.
    When len(defect_coords) == 1, simply returns defect_coords[0].
    First defect_coords is used as a base point when the periodic boundary
    condition is considered.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        defect_entry (DefectEntry):
            Related DefectEntry class object
    """
    # If structure is not given, initial_structure of defect_entry is used.
    if structure is None:
        structure = defect_entry.initial_structure

    inserted_atom_coords = list([structure.frac_coords[k]
                                 for k in defect_entry.inserted_atoms])
    removed_atom_coords = list(defect_entry.removed_atoms.values())

    return defect_center_from_coords(inserted_atom_coords, removed_atom_coords,
                                     structure)


def min_distance_under_pbc(frac1, frac2, lattice_parameters):
    """
    Return the shortest distance between two points in fractional coordinates
    under periodic boundary condition.

    Args:
       frac1 (1x3 np.array):
           1st fractional coordinates
       frac2 (1x3 np.array):
           2nd fractional coordinates
       lattice_parameters (3x3 numpy array):
           a, b, c lattice vectors
    """
    candidate = []
    diff = np.dot(lattice_parameters, frac2 - frac1)

    # (-1, -1, -1), (-1, -1, 0), ..., (1, 1, 1)
    for index in product((-1, 0, 1), repeat=3):
        index = np.array(index)
        delta_diff = np.dot(lattice_parameters, index)
        distance = np.linalg.norm(delta_diff + diff)
        candidate.append([distance, index])

    return min(candidate, key=lambda x: x[0])


def distance_list(structure, coords):
    """
    Return a list of the shortest distances between a point and its images
    under periodic boundary condition.
    Args:
       structure (Structure):
           pmg structure class object
       coords (1x3 numpy array):
           Fractional coordinates
    """
    lattice_parameters = structure.lattice.matrix

    return [min_distance_under_pbc(host, coords, lattice_parameters)[0]
            for host in structure.frac_coords]


def distances_from_point(structure, defect_entry):
    """
    Returns a list of distances at atomic sites from a defect center defined
    by defect_entry. Note that in the case of an interstitial-type defect,
    zero is also set to the interstitial site.

    Args:
        structure (Structure):
            pmg Structure class object for perfect supercell
        defect_entry (DefectEntry):
            DefectEntry class object considered
    """
    return distance_list(structure, defect_center(defect_entry, structure))


class SupercellDftResults:
    """
    A class holding DFT results for supercell systems both w/ and w/o a defect.
    """

    def __init__(self, final_structure, total_energy, magnetization,
                 eigenvalues, kpoints, kpoint_weights, electrostatic_potential,
                 eigenvalue_properties, volume, fermi_level):
        """
        Args:
            final_structure (Structure):
                pmg Structure class object. Usually relaxed structures
            total_energy (float):
                Final total energy in eV.
            magnetization (float):
                Total magnetization in \mu_B
            eigenvalues (N_spin x N_kpoint x N_band np.array):
                Numpy array for the electron eigenvalues in absolute scale.
            kpoints (list):
            kpoint_weights (list):
            electrostatic_potential (list):
                Atomic site electrostatic potential.
            fermi_level (float):
               Fermi level in the absolute scale.
        """
        self._final_structure = final_structure
        self._total_energy = total_energy
        self._magnetization = magnetization
        self._eigenvalues = eigenvalues
        self._kpoints = kpoints
        self._kpoint_weights = kpoint_weights
        self._electrostatic_potential = electrostatic_potential
        self._eigenvalue_properties = eigenvalue_properties
        self._volume = volume
        self._fermi_level = fermi_level

    def __str__(self):
        outs = ["total energy (eV): " + str(self._total_energy),
                "total magnetization (\mu_B): " + str(self._magnetization),
                "electrostatic potential: " + str(self._electrostatic_potential),
                "eigenvalues_properties: " + str(self._eigenvalue_properties),
                "final structure: \n" + str(self._final_structure),
                "volume: \n" + str(self._volume),
                "Fermi level (eV): \n" + str(self._fermi_level),
                "K points \n" + str(self._kpoints)]
        return "\n".join(outs)

    @classmethod
    def from_vasp_files(cls, directory_path, contcar_name="CONTCAR",
                        outcar_name="OUTCAR", vasprun_name="vasprun.xml",
                        continue_forcedly=False):
        """
        Construct class object from vasp output files.
        Args:
            directory_path (str):
                path to the directory storing calc results.
            contcar_name (str):
                Name of the converged CONTCAR file.
            outcar_name (str):
                Name of the OUTCAR file.
            vasprun_name (str):
                Name of the vasprun.xml file.
            continue_forcedly (bool):
                Whether to continue parsing the results even though the
                electronic and/or ionic steps are not converged.
        """

        filename = os.path.join(directory_path, contcar_name)
        try:
            contcar = Poscar.from_file(filename)
        except IOError:
            print("Parsing {] failed.".format(filename))

        filename = os.path.join(directory_path, outcar_name)
        try:
            outcar = Outcar(os.path.join(filename))
        except IOError:
            print("Parsing {] failed.".format(filename))

        filename = os.path.join(directory_path, vasprun_name)
        try:
            vasprun = Vasprun(os.path.join(directory_path, vasprun_name))
        except IOError:
            print("Parsing {] failed.".format(filename))

        # Check if the electronic and ionic steps are converged.
        def managing(continue_or_not, message):
            if continue_or_not:
                raise warnings.warn(message)
            else:
                raise VaspNotConvergedError(message)

        if vasprun.converged_electronic is False:
            managing(continue_forcedly, "Electronic step is not converged.")
        if vasprun.converged_ionic is False:
            managing(continue_forcedly, "Ionic step is not converged.")

        final_structure = contcar.structure
        total_energy = outcar.final_energy
        magnetization = outcar.total_mag
        eigenvalues = vasprun.eigenvalues
        kpoints = vasprun.actual_kpoints
        kpoint_weights = vasprun.actual_kpoints_weights
        electrostatic_potential = outcar.electrostatic_potential
        eigenvalue_properties = vasprun.eigenvalue_band_properties
        volume = contcar.structure.volume
        fermi_level = vasprun.efermi

        return cls(final_structure, total_energy, magnetization, eigenvalues,
                   kpoints, kpoint_weights, electrostatic_potential,
                   eigenvalue_properties, volume, fermi_level)

    @classmethod
    def from_dict(cls, d):
        """
        Construct a class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        return cls(d["final_structure"], d["total_energy"], d["magnetization"],
                   eigenvalues, d["kpoints"], d["kpoint_weights"],
                   d["electrostatic_potential"], d["eigenvalue_properties"],
                   d["volume"], d["fermi_level"])

    @classmethod
    def load_json(cls, filename):
        """
        Construct a class object from a json file.
        defaultdict is imperative to keep the backward compatibility.
        For instance, when one adds new attributes, they do not exist in old
        json files. Then, the corresponding values are set to None.
        """
        dd = defaultdict(lambda: None, loadfn(filename))
        return cls.from_dict(dd)

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
             "magnetization":           self._magnetization,
             "eigenvalues":             eigenvalues,
             "kpoints":                 self._kpoints,
             "kpoint_weights":          self._kpoint_weights,
             "electrostatic_potential": self._electrostatic_potential,
             "eigenvalue_properties":   self._eigenvalue_properties,
             "volume":                  self._volume,
             "fermi_level":             self._fermi_level}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file, named dft_results.json.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @property
    def eigenvalues(self):
        return self._eigenvalues

    @property
    def kpoints(self):
        return self._kpoints

    @property
    def kpoint_weights(self):
        return self._kpoint_weights

    @property
    def final_structure(self):
        return self._final_structure

    @property
    def total_energy(self):
        return self._total_energy

    @property
    def magnetization(self):
        return self._magnetization

    @property
    def electrostatic_potential(self):
        return self._electrostatic_potential

    @property
    def eigenvalue_properties(self):
        return self._eigenvalue_properties

    @property
    def volume(self):
        return self._volume

    @property
    def fermi_level(self):
        return self._fermi_level

    def relative_total_energy(self, referenced_dft_results):
        """
        Return a relative total energy w.r.t. referenced supercell.

        Args:
            referenced_dft_results (SupercellDftResults):
                SupercellDftResults object for referenced supercell dft results.
        """
        return self._total_energy - referenced_dft_results.total_energy

    def relative_potential(self, referenced_dft_results, defect_entry):
        """
        Return a list of relative site potential w.r.t. the perfect supercell.
        Note that None is inserted for interstitial sites.

        Args:
            referenced_dft_results (SupercellDftResults):
                SupercellDftResults object for referenced supercell dft results.
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
                ep_perfect = \
                    referenced_dft_results.electrostatic_potential[p_atom]
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
