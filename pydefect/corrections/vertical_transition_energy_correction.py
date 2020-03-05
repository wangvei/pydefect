# -*- coding: utf-8 -*-
from typing import Optional, List

import matplotlib.pyplot as plt
from monty.json import MSONable
import numpy as np
from numpy import mean

from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.defect_entry import DefectEntry

from pydefect.corrections.corrections import Correction
from pydefect.corrections.efnv_corrections import \
    (Ewald, ExtendedFnvCorrection, calc_lattice_energy_and_pot,
     plot_distance_vs_potential)


"""
This module provides classes used for correcting vertical transition levels
estimated from the total energies.
"""


class VerticalTransitionEnergyCorrection(Correction, MSONable):
    method = "gkfo"

    def __init__(self,
                 charge: int,
                 additional_charge: int,
                 pc_corr_energy: float,
                 ele_pc_corr_energy: float,
                 ave_pot_diff: float,
                 ave_pot_diff_by_addition: float,
                 dielectric_tensor: np.ndarray,
                 electronic_dielectric_tensor: np.ndarray,
                 symbols_without_defect: List[str],
                 distances_from_defect: List[float],
                 relative_potential: List[float],
                 pc_pot: List[float],
                 max_sphere_radius: float,
                 ):
        """Calculate the correction energy for the vertical transition.

        [Gake, Kumagai, Freysoldt and Oba, PRB 101, 020102 (2020)] for details.
        Here, "GKFO" means the correction energy of Eq. (13) in this paper.

        Args:
            charge (int):
                Charge state before the vertical transition (Q in GKFO).
            additional_charge (int):
                Additional charge, so must be 1 or -1 (Delta Q in GKFO).
            pc_corr_energy (float):
                Point-charge correction energy including ionic contribution
                (E^Q_{rm PC}(epsilon_0) in GKFO).
            ele_pc_corr_energy (float):
                Point-charge (PC) correction energy screened by only electronic
                dielectric constant for an additional charge
                (E^{Delta Q}_{\rm PC}(epsilon_infty) in GKFO).
            ave_pot_diff (float):
                Average potential difference before the transition
                (C^Q in GKFO correction).
            ave_pot_diff_by_addition (float):
                Average potential difference caused by adding charge Delta Q
                (C^{Delta Q} in GKFO correction).
            dielectric_tensor (3x3 np.array):
                Static dielectric tensor where the directions are compatible
                with the lattice (epsilon_0 in GKFO).
            electronic_dielectric_tensor (3x3 np.array):
                Dielectric tensor of electronic part (epsilon_infty in GKFO).
        """
        if abs(additional_charge) != 1:
            raise ValueError(f"{additional_charge} is invalid.")

        self.charge = charge
        self.additional_charge = additional_charge
        self.pc_corr_energy = pc_corr_energy
        self.ele_pc_corr_energy = ele_pc_corr_energy

        self.ave_pot_diff = ave_pot_diff
        self.ave_pot_diff_by_addition = ave_pot_diff_by_addition

        self.dielectric_tensor = dielectric_tensor
        self.electronic_dielectric_tensor = electronic_dielectric_tensor

        self.symbols_without_defect = symbols_without_defect
        self.distances_from_defect = distances_from_defect
        self.relative_potential = relative_potential
        self.pc_pot = pc_pot
        self.max_sphere_radius = max_sphere_radius

    @classmethod
    def from_files(cls,
                   dielectric_tensor: np.ndarray,
                   static_dielectric_tensor: np.ndarray,
                   initial_efnv: ExtendedFnvCorrection,
                   initial_calc_results: SupercellCalcResults,
                   final_defect_entry: DefectEntry,
                   final_calc_results: SupercellCalcResults,
                   prod_cutoff_fwhm: float = 25.0
                   ) -> "VerticalTransitionEnergyCorrection":
        """Create instance object by extracting data from several objects.

        Args:
            dielectric_tensor (np.ndarray):
                Sum of static and ionic dielectric tensor.
            static_dielectric_tensor (np.ndarray):
                static dielectric tensor.
            initial_efnv:
                ExtendedFnvCorrection instance object for initial charge state.
            initial_calc_results:
                SupercellCalcResults instance object for initial charge state.
            final_defect_entry:
                DefectEntry instance object for initial charge state.
            final_calc_results:
                SupercellCalcResults instance object for final charge state.
            prod_cutoff_fwhm:
                product of cutoff radius of G-vector and gaussian FWHM.
                Increasing this value, calculation will be more accurate, but
                slower.

        Returns:
            VerticalTransitionEnergyCorrection instance object.
        """
        added_charge = final_defect_entry.charge - initial_efnv.charge

        pc_corr_energy = initial_efnv.point_charge_correction_energy
        ave_pot_diff = initial_efnv.ave_pot_diff

        ewald = Ewald.from_optimization(final_calc_results.final_structure,
                                        static_dielectric_tensor,
                                        prod_cutoff_fwhm=prod_cutoff_fwhm)

        atom_coords = initial_efnv.atomic_coords_without_defect
        center_coords = initial_efnv.defect_center_coords
        lattice = final_calc_results.final_structure.lattice
        ele_pc_energy, model_pot = calc_lattice_energy_and_pot(
            atom_coords, added_charge, center_coords, ewald, lattice)

        initial_pot = initial_calc_results.electrostatic_potential
        final_pot = final_calc_results.electrostatic_potential
        # Need to reverse vasp potential.
        relative_potential = [-(f - i) for i, f in zip(initial_pot, final_pot)]

        distances = initial_efnv.distances_from_defect
        pot_diff = []
        for (d, a, m) in zip(distances, relative_potential, model_pot):
            if d > initial_efnv.defect_region_radius:
                pot_diff.append(a - m)
        ave_pot_diff_by_addition = float(mean(pot_diff))

        return cls(charge=initial_efnv.charge,
                   additional_charge=added_charge,
                   pc_corr_energy=pc_corr_energy,
                   ele_pc_corr_energy=-ele_pc_energy,
                   ave_pot_diff=ave_pot_diff,
                   ave_pot_diff_by_addition=ave_pot_diff_by_addition,
                   dielectric_tensor=dielectric_tensor,
                   electronic_dielectric_tensor=static_dielectric_tensor,
                   symbols_without_defect=initial_efnv.symbols_without_defect,
                   distances_from_defect=initial_efnv.distances_from_defect,
                   relative_potential=relative_potential,
                   pc_pot=model_pot,
                   max_sphere_radius=initial_efnv.max_sphere_radius)

    def plot_potential(self,
                       file_name: str = "vertical_correction.pdf",
                       yrange: Optional[list] = None) -> None:
        """Plotter for the potential as a function of distance."""
        plot_distance_vs_potential(self.symbols_without_defect,
                                   self.distances_from_defect,
                                   self.relative_potential,
                                   self.pc_pot,
                                   self.max_sphere_radius,
                                   self.ave_pot_diff,
                                   yrange)
        plt.savefig(file_name, format="pdf", transparent=True)

    def __repr__(self):
        outs = \
            [f"Charge : {self.charge}",
             f"Added charge: {self.additional_charge} ", "",
             f"PC 1st: {self.pc_corr_first}",
             f"PC 2nd: {self.pc_corr_second}",
             f"Point-charge contribution: {self.pc_correction_energy}", "",
             f"Alignment 1st: {self.align_first}",
             f"Alignment 2nd: {self.align_second}",
             f"Alignment 3rd: {self.align_third}",
             f"Alignment contribution: {self.alignment_correction_energy}", "",
             f"Total correction: {self.correction_energy}", "",
             ]

        return "\n".join(outs)

    @property
    def correction_energy(self) -> float:
        return self.pc_correction_energy + self.alignment_correction_energy

    @property
    def pc_corr_first(self) -> float:
        return 2 * self.additional_charge / self.charge * self.pc_corr_energy

    @property
    def pc_corr_second(self) -> float:
        return self.ele_pc_corr_energy

    @property
    def pc_correction_energy(self) -> float:
        return self.pc_corr_first + self.pc_corr_second

    @property
    def align_first(self) -> float:
        return - self.additional_charge * self.ave_pot_diff_by_addition

    @property
    def align_second(self) -> float:
        return - self.additional_charge * self.ave_pot_diff

    @property
    def align_third(self) -> float:
        ave_dielectric = np.trace(self.dielectric_tensor) / 3
        ave_ele_dielectric = np.trace(self.electronic_dielectric_tensor) / 3
        dielectric_ratio = ave_ele_dielectric / ave_dielectric
        return - dielectric_ratio * self.charge * self.ave_pot_diff_by_addition

    @property
    def alignment_correction_energy(self) -> float:
        return self.align_first + self.align_second + self.align_third
