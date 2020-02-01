from monty.json import MSONable
import numpy as np
from numpy import mean
from typing import Optional

from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.defect_entry import DefectEntry

from pydefect.corrections.corrections import Correction
from pydefect.corrections.efnv_corrections import \
    (Ewald, ExtendedFnvCorrection, calc_lattice_energy_and_pot)

from pymatgen.core.structure import Structure

"""
This module provides classes used for correcting vertical transition levels
estimated from the total energies.
"""

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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
                 electronic_dielectric_tensor: np.ndarray):
        """Calculate the correction energy for the vertical transition

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

    @classmethod
    def from_files(cls,
                   unitcell: UnitcellCalcResults,
                   initial_efnv_cor: ExtendedFnvCorrection,
                   initial_calc_results: SupercellCalcResults,
                   final_defect_entry: DefectEntry,
                   final_calc_results: SupercellCalcResults,
                   prod_cutoff_fwhm: float = 25.0):
        """

        """

        static_dielectric_tensor = unitcell.static_dielectric_tensor
        added_charge = final_defect_entry.charge - initial_efnv_cor.charge

        pc_corr_energy = initial_efnv_cor.point_charge_correction_energy
        ave_pot_diff = initial_efnv_cor.ave_pot_diff

        ewald = Ewald.from_optimization(final_calc_results.final_structure,
                                        unitcell.static_dielectric_tensor,
                                        prod_cutoff_fwhm=prod_cutoff_fwhm)

        atom_coords = initial_efnv_cor.atomic_coords_without_defect
        center_coords = initial_efnv_cor.defect_center_coords
        lattice = final_calc_results.final_structure.lattice
        ele_pc_energy, model_pot = calc_lattice_energy_and_pot(atom_coords,
                                                               added_charge,
                                                               center_coords,
                                                               ewald,
                                                               lattice)
        relative_potential = []
        initial_pot = initial_calc_results.electrostatic_potential
        final_pot = initial_calc_results.electrostatic_potential
        for i, j in zip(initial_pot, final_pot):
            # Need to reverse vasp potential.
            relative_potential.append(-(j - i))

        distances = initial_efnv_cor.distances_from_defect
        pot_diff = []
        for (d, a, m) in zip(distances, relative_potential, model_pot):
            if d > initial_efnv_cor.defect_region_radius:
                pot_diff.append(a - m)
        ave_pot_diff_by_addition = float(mean(pot_diff))

        return cls(charge=initial_efnv_cor.charge,
                   additional_charge=added_charge,
                   pc_corr_energy=pc_corr_energy,
                   ele_pc_corr_energy=-ele_pc_energy,
                   ave_pot_diff=ave_pot_diff,
                   ave_pot_diff_by_addition=ave_pot_diff_by_addition,
                   dielectric_tensor=unitcell.total_dielectric_tensor,
                   electronic_dielectric_tensor=static_dielectric_tensor)

    def __repr__(self):
        outs = [f"Charge : {self.charge}",
                f"Added charge: {self.additional_charge} ",
                f"Total correction: {self.correction_energy}",
                f"Point-charge contribution: {self.pc_correction_energy}",
                f"Alignment contribution: {self.alignment_correction_energy}"]

        return "\n".join(outs)

    @property
    def correction_energy(self) -> float:
        return self.pc_correction_energy + self.alignment_correction_energy

    @property
    def pc_correction_energy(self) -> float:

        first_term = \
            2 * self.additional_charge / self.charge * self.pc_corr_energy
        second_term = self.ele_pc_corr_energy
        print("first_term")
        print(first_term)
        print("second_term")
        print(second_term)
        return first_term + second_term

    @property
    def alignment_correction_energy(self) -> float:
        ave_dielectric = np.trace(self.dielectric_tensor) / 3
        ave_ele_dielectric = np.trace(self.electronic_dielectric_tensor) / 3

        return - (self.additional_charge * self.ave_pot_diff_by_addition
                  + self.additional_charge * self.ave_pot_diff
                  + ave_ele_dielectric / ave_dielectric
                  * self.charge * self.ave_pot_diff_by_addition)
