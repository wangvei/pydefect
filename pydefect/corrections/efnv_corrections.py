# -*- coding: utf-8 -*-
import json
from copy import deepcopy
from functools import reduce
from itertools import product, groupby
from math import sqrt, pow, ceil
from operator import itemgetter
from typing import Optional, Union, Tuple, List

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn
from numpy import cos, sqrt, dot, cross, pi, exp, mean
from numpy.linalg import norm

from pydefect.core.config import COLOR
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.corrections.corrections import Correction
from pydefect.util.logger import get_logger

from pymatgen.core import Structure
from pymatgen.core.lattice import Lattice

from scipy.constants import elementary_charge, epsilon_0
from scipy.special import erfc
from scipy.stats import mstats

"""
This module provides functions used to correct defect formation energies
with finite supercell-size dependencies.
"""

logger = get_logger(__name__)


def calc_max_sphere_radius(lattice_vectors: np.ndarray) -> float:
    """Calculate Maximum radius of a sphere fitting inside the unit cell.

    Calculate three distances between two parallel planes using
    (a_i x a_j) . a_k / |a_i . a_j|

    Args:
        lattice_vectors (np.ndarray): 3x3 lattice vectors.

    Returns:
        Float of the radius
    """
    distances = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_a_j = cross(lattice_vectors[i - 2], lattice_vectors[i - 1])
        a_k = lattice_vectors[i]
        distances[i] = abs(dot(a_i_a_j, a_k)) / norm(a_i_a_j)
    return max(distances) / 2.0


def create_lattice_set(lattice_vectors: np.ndarray,
                       max_length: float) -> List[list]:
    """Return a set of lattice vectors within the max length.

    Note that angles between any two axes are assumed to be between 60 and
    120 degrees.

    Args:
        lattice_vectors (np.ndarray): 3x3 matrix.
        max_length (float): Max length to search lattice set.

    Returns:
        List of list of Cartesian vectors of lattice points.
    """
    max_int = [ceil(max_length / norm(lattice_vectors[i])) for i in range(3)]
    range_list = [range(-mi, mi + 1) for mi in max_int]

    lattice_set = []
    for index in product(range_list[0], range_list[1], range_list[2]):
        cart_vector = dot(lattice_vectors.transpose(), np.array(index))
        if norm(cart_vector) < max_length:
            lattice_set.append(cart_vector.tolist())

    return lattice_set


def calc_relative_potential(defect: SupercellCalcResults,
                            perfect: SupercellCalcResults,
                            defect_entry: DefectEntry) -> List[float]:
    """Return a list of relative site potential w.r.t. the perfect supercell.

    None is inserted for foreign atoms such as interstitials and substituted
    atoms.

    Args:
        defect (SupercellCalcResults):
            SupercellDftResults object for defect supercell dft results.
        perfect (SupercellCalcResults):
            SupercellDftResults object for referenced supercell dft results.
            Usually it is for perfect supercell.
        defect_entry (DefectEntry):
            DefectEntry class object corresponding.

    Returns:
        List of floats showing relative potential.
    """
    mapping = defect_entry.atom_mapping_to_perfect
    relative_potential = list()

    for d_atom, p_atom in enumerate(mapping):

        if p_atom is None:
            relative_potential.append(None)
        else:
            potential_defect = defect.electrostatic_potential[d_atom]
            potential_perfect = perfect.electrostatic_potential[p_atom]
            relative_potential.append(potential_defect - potential_perfect)

    return relative_potential


class Ewald(MSONable):
    """Container class for anisotropic Ewald sum."""

    def __init__(self,
                 lattice: Lattice,
                 dielectric_tensor: np.ndarray,
                 ewald_param: float,
                 prod_cutoff_fwhm: float,
                 real_neighbor_lattices: list,
                 reciprocal_neighbor_lattices: list):
        """ Store Ewald parameter and its related properties

        Args:
            lattice (Lattice):
                The given lattice.
            dielectric_tensor (3x3 np.array):
                Static dielectric tensor where the directions are compatible
                with the lattice.
            ewald_param (float):
                A parameter used for evaluating Ewald sum.
            prod_cutoff_fwhm (float):
                product of cutoff radius of G-vector and gaussian FWHM.
                Increasing this value, calculation will be more accurate, but
                slower.
            real_neighbor_lattices (list):
                List of lattice vectors in real space.
            reciprocal_neighbor_lattices (list):
                List of reciprocal lattice vectors.
        """

        self.lattice = lattice
        self.volume = lattice.volume
        self.reciprocal_lattice_matrix = self.lattice.reciprocal_lattice.matrix
        self.dielectric_tensor = dielectric_tensor
        self.ewald_param = ewald_param
        self.prod_cutoff_fwhm = prod_cutoff_fwhm

        # including self lattice site
        self.real_neighbor_lattices = real_neighbor_lattices[:]
        self.reciprocal_neighbor_lattices = reciprocal_neighbor_lattices[:]

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    @property
    def max_r_vector_norm(self):
        root_det_dielectric = sqrt(np.linalg.det(self.dielectric_tensor))
        cube_root_vol = pow(self.volume, 1 / 3)
        ewald = self.ewald_param / cube_root_vol * root_det_dielectric

        return self.prod_cutoff_fwhm / ewald

    @property
    def max_g_vector_norm(self):
        root_det_dielectric = sqrt(np.linalg.det(self.dielectric_tensor))
        cube_root_vol = pow(self.volume, 1 / 3)
        ewald = self.ewald_param / cube_root_vol * root_det_dielectric

        return 2 * ewald * self.prod_cutoff_fwhm

    def to_json_file(self, filename: str) -> None:
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def from_optimization(cls,
                          structure: Structure,
                          dielectric_tensor: np.ndarray,
                          initial_ewald_param: float = None,
                          convergence: float = 1.05,
                          prod_cutoff_fwhm: float = 25.0):
        """ Get optimized ewald parameter.

        Args:
            structure (Structure):
                Structure
            dielectric_tensor (numpy 3x3 array):
                dielectric tensor
            initial_ewald_param (float):
                Initial guess of parameter.
            convergence (float):
                If 1/convergence < n_(real)/n_(reciprocal) < convergence,
                where n_(real) and n_(reciprocal) is number of real lattices
                and reciprocal lattices, finishes optimization and
                returns ewald_param.
            prod_cutoff_fwhm (float):
                product of cutoff radius of G-vector and gaussian FWHM.
                Increasing this value, calculation will be more accurate, but
                slower.

        Returns:
            Optimized Ewald parameter (float).
        """
        root_det_dielectric = sqrt(np.linalg.det(dielectric_tensor))
        real_lattice_matrix = structure.lattice.matrix
        reciprocal_lattice_matrix = structure.lattice.reciprocal_lattice.matrix
        cube_root_vol = pow(structure.lattice.volume, 1 / 3)

        if initial_ewald_param is not None:
            ewald_param = initial_ewald_param
        else:
            # determine initial ewald parameter to satisfy following:
            # max_int(Real) = max_int(Reciprocal)
            # in neighbor_lattices function.

            # Left term:
            # max_int(Real) = 2 * x * Y  / l_r, where x, Y, and l_r are ewald,
            # prod_cutoff_fwhm, and axis length of real lattice, respectively.

            # Right term:
            # max_int(reciprocal) = Y  / (x * l_g)
            # where l_g is axis length of reciprocal lattice, respectively.
            # Then, x = sqrt(l_g / l_r / 2)
            # gmean : geometric mean,  (a1 * a2 * a3)^(1/3)

            l_r = mstats.gmean([norm(v) for v in real_lattice_matrix])
            l_g = mstats.gmean([norm(v) for v in reciprocal_lattice_matrix])
            ewald_param \
                = sqrt(l_g / l_r / 2) * cube_root_vol / root_det_dielectric

        for i in range(10):
            ewald = ewald_param / cube_root_vol * root_det_dielectric
            # count number of real lattice
            max_r_vector_norm = prod_cutoff_fwhm / ewald
            real_neighbor_lattices = \
                create_lattice_set(real_lattice_matrix, max_r_vector_norm)

            # count number of reciprocal lattice
            max_g_vector_norm = 2 * ewald * prod_cutoff_fwhm

            reciprocal_neighbor_lattices = \
                create_lattice_set(reciprocal_lattice_matrix,
                                   max_g_vector_norm)

            real_reciprocal_ratio = \
                len(real_neighbor_lattices) / len(reciprocal_neighbor_lattices)

            if 1 / convergence < real_reciprocal_ratio < convergence:
                return cls(structure.lattice, dielectric_tensor, ewald_param,
                           prod_cutoff_fwhm, real_neighbor_lattices,
                           reciprocal_neighbor_lattices)
            else:
                ewald_param *= real_reciprocal_ratio ** 0.17
        else:
            raise ValueError("The initial ewald param may not be adequate.")

    def reciprocal_lattice_set(self) -> np.ndarray:
        """Return reciprocal lattice set."""
        # remove [0, 0, 0] G vector.
        lattices = np.array(self.reciprocal_neighbor_lattices)
        # Note: numpy.delete returns new object.
        return np.delete(lattices, int((len(lattices) - 1) / 2), 0)

    def real_lattice_set(self,
                         include_self: bool = True,
                         shift: np.array = np.zeros(3)) -> np.ndarray:
        """Return real lattice set.
        Args:
             include_self (bool):
                Whether to include the central lattice point itself.
             shift (np.array):
        """
        lattices = np.array(deepcopy(self.real_neighbor_lattices))

        if not include_self:
            # Remove the middle constituent.
            # Note: numpy.delete returns new object.
            lattices = np.delete(lattices, int((len(lattices) - 1) / 2), 0)

        return lattices + shift


class ExtendedFnvCorrection(Correction, MSONable):
    """Calculate and manage the plot for the extended FNV corrections

    See [Kumagai and Oba, PRB 89, 195205 (2014)] for details.
    """
    method = "extended_FNV"

    def __init__(self,
                 ewald_json: str,
                 charge: int,
                 lattice_matrix: np.ndarray,
                 lattice_energy: float,
                 ave_pot_diff: float,
                 symbols_without_defect: List[str],
                 defect_center_coords: List[float],
                 atomic_coords_without_defect: List[np.ndarray],
                 distances_from_defect: List[float],
                 electrostatic_pot: List[float],
                 pc_pot: List[float],
                 defect_region_radius: float):
        """
        Args:
            ewald_json (str):
                Since size of Ewald class instance object is large, filename is
                stored instead.
            charge (int):
                Defect charge state.
            lattice_matrix (np.ndarray):
                3x3 lattice matrix [[a_x, a_y, a_z], [b_x,..], [c_x,..]].
            lattice_energy (float):
                Lattice energy at a given charge and lattice_matrix.
            ave_pot_diff (float):
                Averaged potential difference at atoms outside of
                defect_region_radius.
            symbols_without_defect (list of str):
                List of element symbols for all the atomic sites except for the
                defect.
                e.g., ["Mg", "Mg", ..., "O", "O", ..]
            defect_center_coords (list):
                Defect center position in fractional coordinates.
            atomic_coords_without_defect (list):
                List of fractional coordinates for all the atomic sites except
                for the defect.
                e.g., [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], ...]
            distances_from_defect (list of float):
                List of distances from the defect for all the atomic sites
            electrostatic_pot (list of float):
                List of electrostatic potential w.r.t the perfect one from
                the defect for all the atomic sites.
            pc_pot (list of float):
                List of point-charge potential from the defect for all the
                atomic sites.
            defect_region_radius (float):
                Maximum radius of a sphere touching to the lattice plane, used
                for defining the outside region of the defect.
        """
        # error check
        symbol_len = len(symbols_without_defect)
        dist_len = len(distances_from_defect)
        ele_pot_len = len(electrostatic_pot)
        pc_pot_len = len(pc_pot)
        if not (symbol_len == dist_len == ele_pot_len == pc_pot_len):
            raise IndexError(
                f"Lengths of symbols({symbol_len}), distances({dist_len}),  "
                f"electrostatic_pot({ele_pot_len}), model_pot({pc_pot_len}) "
                f"differ.")

        self.ewald_json = ewald_json
        self.charge = charge
        self.lattice_matrix = lattice_matrix
        self.lattice_energy = lattice_energy
        self.ave_pot_diff = ave_pot_diff
        self.defect_center_coords = defect_center_coords
        self.symbols_without_defect = symbols_without_defect
        self.atomic_coords_without_defect = atomic_coords_without_defect
        self.distances_from_defect = distances_from_defect[:]
        self.electrostatic_pot = electrostatic_pot[:]
        self.pc_pot = pc_pot[:]
        self.defect_region_radius = defect_region_radius
        self._manually_added_correction_energy = 0.0

    def __repr__(self):
        outs = \
            [f"Point-charge correction (eV): "
             f"{self.point_charge_correction_energy}",
             f"Alignment-like correction (eV): "
             f"{self.alignment_correction_energy}",
             f"Manually added correction (eV): "
             f"{self.manually_added_correction_energy}",
             f"Total correction energy (eV): {self.correction_energy}"]
        return "\n".join(outs)

    @property
    def point_charge_correction_energy(self) -> float:
        return - self.lattice_energy

    @property
    def alignment_correction_energy(self) -> float:
        return - self.ave_pot_diff * self.charge

    @property
    def manually_added_correction_energy(self) -> float:
        return self._manually_added_correction_energy

    @manually_added_correction_energy.setter
    def manually_added_correction_energy(self, corr_energy: float) -> None:
        self._manually_added_correction_energy = corr_energy

    @property
    def correction_energy(self) -> float:
        return (self.point_charge_correction_energy
                + self.alignment_correction_energy
                + self.manually_added_correction_energy)

    @property
    def max_sphere_radius(self) -> float:
        return calc_max_sphere_radius(self.lattice_matrix)

    def plot_potential(self,
                       file_name: str,
                       yrange: Optional[list] = None) -> None:
        """Plotter for the potential as a function of distance."""
        plot_distance_vs_potential(self.symbols_without_defect,
                                   self.distances_from_defect,
                                   self.electrostatic_pot,
                                   self.pc_pot,
                                   self.max_sphere_radius,
                                   self.ave_pot_diff,
                                   yrange)
        plt.tight_layout()
        plt.savefig(file_name, format="pdf", transparent=True)

    @classmethod
    def compute_correction(cls,
                           defect_entry: DefectEntry,
                           defect_dft: SupercellCalcResults,
                           perfect_dft: SupercellCalcResults,
                           dielectric_tensor: np.ndarray,
                           defect_center: list = None,
                           ewald: Union[str, Ewald] = None,
                           to_filename: str = "ewald.json"
                           ) -> "ExtendedFnvCorrection":
        """ Estimate correction energy for point defect formation energy.

        Args:
            defect_entry (DefectEntry):
                DefectEntry object of the considering defect.
            defect_dft (SupercellCalcResults):
                Calculated defect DFT results of the considering defect.
            perfect_dft (SupercellCalcResults):
                Calculated defect DFT results for perfect supercell.
            dielectric_tensor (np.ndarray):
                Dielectric tensor.
            defect_center (list):
                Defect center in fractional coordinates.
            ewald (str / Ewald):
                Ewald object or ewald.json filename.
            to_filename (str):
                Filename to jump json data.

        Returns
            ExtendedFnvCorrection class instance object.
        """
        if isinstance(ewald, str):
            ewald = Ewald.load_json(ewald)
        elif ewald is None:
            ewald = \
                Ewald.from_optimization(defect_dft.final_structure,
                                        dielectric_tensor)
            ewald.to_json_file(to_filename)

        relative_potential = calc_relative_potential(defect=defect_dft,
                                                     perfect=perfect_dft,
                                                     defect_entry=defect_entry)

        # not None is a must as 0 is also judged as False.
        diff_potential = [-ep for ep in relative_potential if ep is not None]

        defective_structure = defect_dft.final_structure

        atomic_coords_without_defect = []
        symbols_without_defect = []
        for i, j in enumerate(defect_entry.atom_mapping_to_perfect):
            if j is not None:
                atomic_coords_without_defect.append(
                    defective_structure.frac_coords[i])
                symbols_without_defect.append(
                    defective_structure.sites[i].specie.symbol)

        charge = defect_entry.charge
        lattice = defect_dft.final_structure.lattice
        if defect_center:
            logger.warning(f"Defect center {defect_center} is explicitly set."
                           f"User needs to know what's going on.")
        else:
            defect_center = defect_entry.defect_center_coords

        distances_from_defect = \
            deepcopy(defect_dft.displacements["final_distances"])

        inserted_atom_indices = \
            [i["index"] for i in defect_entry.inserted_atoms]
        for i in sorted(inserted_atom_indices, reverse=True):
            del distances_from_defect[i]

        lattice_energy, model_pot = \
            calc_lattice_energy_and_pot(atomic_coords_without_defect,
                                        charge,
                                        defect_center,
                                        ewald,
                                        lattice)
        # calc ave_pot_diff
        defect_region_radius = calc_max_sphere_radius(np.array(lattice.matrix))
        pot_diff = []

        # error check just in case (should be removed in the future)
        if not (len(distances_from_defect) == len(diff_potential)
                == len(model_pot)):
            raise IndexError(
                f"Size of distances_from_defect({len(distances_from_defect)}),"
                f" diff_ep({len(diff_potential)}), "
                f"model_pot({len(model_pot)}) differ.")

        for (d, a, m) in zip(distances_from_defect, diff_potential, model_pot):
            if d > defect_region_radius:
                pot_diff.append(a - m)
        ave_pot_diff = float(mean(pot_diff))

        return cls(ewald_json=ewald,
                   charge=charge,
                   lattice_matrix=lattice.matrix,
                   lattice_energy=lattice_energy,
                   ave_pot_diff=ave_pot_diff,
                   symbols_without_defect=symbols_without_defect,
                   defect_center_coords=defect_center,
                   atomic_coords_without_defect=atomic_coords_without_defect,
                   distances_from_defect=distances_from_defect,
                   electrostatic_pot=diff_potential,
                   pc_pot=model_pot,
                   defect_region_radius=defect_region_radius)


def calc_lattice_energy_and_pot(atomic_coords_without_defect: List[np.ndarray],
                                charge: int,
                                defect_center_coords: List[float],
                                ewald: "Ewald",
                                lattice: Lattice) -> Tuple[float, List[float]]:
    """

    Args:
        atomic_coords_without_defect (List):
            List of fractional coordinates for all the atomic sites except
            for the defect.
            e.g., [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], ...]
        charge:
            Defect charge state.
        defect_center_coords:
            Defect center position in fractional coordinates.
        ewald:
            Ewald object.
        lattice:
            Lattice object for the supercell.

    Returns:
        Tuple of the lattice energy and list of the model potential.
    """

    volume = lattice.volume
    coeff, diff_pot, mod_ewald_param, root_det_epsilon = \
        constants_for_anisotropic_ewald_sum(charge, ewald, volume)
    # model potential and lattice energy
    model_pot = []
    for r in atomic_coords_without_defect:
        # Ewald real part
        # \sum erfc(ewald*\sqrt(R*\epsilon_inv*R))
        # / \sqrt(det(\epsilon)) / \sqrt(R*\epsilon_inv*R) [1/A]
        shift = lattice.get_cartesian_coords(r - defect_center_coords)

        real_part, reciprocal_part = \
            calc_ewald_sum(ewald=ewald,
                           mod_ewald_param=mod_ewald_param,
                           root_det_epsilon=root_det_epsilon,
                           volume=volume,
                           include_self=True,
                           shift=shift)

        model_pot.append((real_part + reciprocal_part + diff_pot) * coeff)
    lattice_energy = point_charge_energy(charge, ewald, volume)
    return lattice_energy, model_pot


def point_charge_energy(charge: int,
                        ewald: "Ewald",
                        volume: float) -> float:
    """Return point-charge energy under periodic boundary condition."""
    if charge == 0:
        return 0.0

    coeff, diff_pot, mod_ewald_param, root_det_epsilon = \
        constants_for_anisotropic_ewald_sum(charge, ewald, volume)
    # Real part: sum erfc(ewald * sqrt(R * epsilon_inv * R))
    #                    / sqrt(det(epsilon)) / sqrt(R * epsilon_inv * R) [1/A]
    real_part, reciprocal_part = \
        calc_ewald_sum(ewald, mod_ewald_param, root_det_epsilon, volume)

    det_epsilon = np.linalg.det(ewald.dielectric_tensor)
    self_pot = - mod_ewald_param / (2.0 * pi * sqrt(pi * det_epsilon))
    lattice_energy = ((real_part + reciprocal_part + diff_pot + self_pot)
                      * coeff * charge / 2)

    return lattice_energy


def calc_ewald_sum(ewald: "Ewald",
                   mod_ewald_param: float,
                   root_det_epsilon: np.ndarray,
                   volume: float,
                   include_self: bool = False,
                   shift: np.ndarray = np.array([0, 0, 0])
                   ) -> Tuple[float, float]:
    """Return real and reciprocal Ewald summations at given parameters"""

    epsilon_inv = np.linalg.inv(ewald.dielectric_tensor)
    real_sum = 0
    # Skip the potential caused by the defect itself
    for v in ewald.real_lattice_set(include_self, shift):
        root_r_inv_epsilon_r = np.sqrt(reduce(dot, [v.T, epsilon_inv, v]))
        real_sum += \
            erfc(mod_ewald_param * root_r_inv_epsilon_r) / root_r_inv_epsilon_r
    real_part = real_sum / (4 * pi * root_det_epsilon)

    # Ewald reciprocal part
    # sum exp(-g * epsilon * g / (4 * ewald ** 2)) / g * epsilon * g [1/A]
    reciprocal_sum = 0
    for g in ewald.reciprocal_lattice_set():
        g_epsilon_g = reduce(dot, [g.T, ewald.dielectric_tensor, g])
        reciprocal_sum += \
            (exp(- g_epsilon_g / 4.0 / mod_ewald_param ** 2)
             / g_epsilon_g * cos(dot(g, np.zeros(3))))  # [A^2]
    reciprocal_part = reciprocal_sum / volume

    return real_part, reciprocal_part


def constants_for_anisotropic_ewald_sum(charge: int,
                                        ewald: "Ewald",
                                        volume: float) -> Tuple[float,
                                                                float,
                                                                float,
                                                                np.ndarray]:
    """Derive some constants used for anisotropic Ewald sum.

    YK2014: Kumagai and Oba, PRB 89, 195205 (2014)
    Note that in this program the formula written in YK2014 are divided by 4pi
    to keep the SI unit.

    Args:
        charge (int): Point charge
        ewald (Ewald): Ewald object with dielectric_tensor and Ewald parameter.
        volume (float): Volume in [A^3]

    Returns: Tuple of the following float values:
        coeff:
            Common coefficient in anisotropic Ewald sum.
        diff_pot:
            2nd term in Eq.(14) in YK2014, which comes from the potential
            difference caused by the finite size gaussian charge.
        mod_ewald_param:
            Modified Ewald parameter which is the ita in the Eqs in YK2014.
        root_det_epsilon:
            Square root of determinant of dielectric tensor
    """
    coeff = charge * elementary_charge * 1e10 / epsilon_0  # [V]
    cube_root_vol = pow(volume, 1 / 3)
    root_det_epsilon = sqrt(np.linalg.det(ewald.dielectric_tensor))
    mod_ewald_param = ewald.ewald_param / cube_root_vol * root_det_epsilon
    diff_pot = -0.25 / volume / mod_ewald_param ** 2  # [1/A]

    return coeff, diff_pot, mod_ewald_param, root_det_epsilon


def plot_distance_vs_potential(symbols_without_defect,
                               distances_from_defect,
                               difference_electrostatic_pot,
                               model_pot,
                               max_sphere_radius,
                               ave_pot_diff,
                               yrange: Optional[list] = None) -> None:
    """Plotter for the potential as a function of distance."""
    property_without_defect = list(zip(symbols_without_defect,
                                       distances_from_defect,
                                       difference_electrostatic_pot))
    # E.g. points_dictionary = {'Mg': [(3.67147, -0.7019),  ..], 'O': [..]}
    points_dictionary = {}
    for element, properties in \
            groupby(property_without_defect, key=itemgetter(0)):
        points_dictionary[element] = [(x, y) for _, x, y in properties]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # DFT electrostatic potential
    for i, (symbol, points) in enumerate(points_dictionary.items()):
        x_set = np.array([point[0] for point in points])
        y_set = np.array([point[1] for point in points])
        ax.scatter(x_set, y_set, c=COLOR[i], marker="x", label=symbol)

    # PC model potential
    ax.scatter(distances_from_defect, model_pot,
               marker=".", label="model potential", color=COLOR[-1])

    # difference between PC model potential and DFT electrostatic potential
    diff_model_electrostatic = (np.array(difference_electrostatic_pot)
                                - np.array(model_pot))

    ax.scatter(distances_from_defect, diff_model_electrostatic,
               marker="o", label="potential diff",
               facecolors='none', edgecolors=COLOR[-2])

    # potential difference
    point_x = [max_sphere_radius, max(distances_from_defect)]
    point_y = [ave_pot_diff, ave_pot_diff]

    ax.set_xlabel(r"Distance from a defect (${\rm \AA}$)", fontsize=15)
    ax.set_ylabel("Electrostatic potential (V)", fontsize=15)

    ax.plot(point_x, point_y, c=(0, 0, 0), label="potential difference")
    ax.legend(bbox_to_anchor=(1, 0.5), loc='center left')
    ax.tick_params(
        direction='in', bottom=True, top=True, left=True, right=True)

    # change 0.0 to 0
    from pydefect.util.matplotlib import formatter
    ax.xaxis.set_major_formatter(formatter)
    # ax.yaxis.set_major_formatter(formatter)

    if yrange:
        plt.ylim(yrange[0], yrange[1])

