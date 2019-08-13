from monty.json import MontyEncoder, MSONable
from pydefect.corrections.corrections import Correction
from pydefect.corrections.efnv_corrections import Ewald, ExtendedFnvCorrection
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pymatgen.core.structure import Structure
from pydefect.corrections.efnv_corrections import point_charge_energy

"""
This module provides classes used for correcting vertical transition levels
estimated from the total energies.
"""

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class GkoCorrection(Correction, MSONable):
    method = "GKO"

    def __init__(self,
                 before_charge: int,
                 after_charge: int,
                 ewald_electronic_part: Ewald,
                 ewald_ionic_part: Ewald,
                 before_electronic_pc_correction: float,
                 after_electronic_pc_correction: float,
                 before_ionic_pc_correction: float,
                 ave_pot_diff: float):

        self.before_charge = before_charge
        self.after_charge = after_charge
        self.ewald_electronic_part = ewald_electronic_part
        self.ewald_ionic_part = ewald_ionic_part
        self.before_electronic_pc_correction = before_electronic_pc_correction
        self.after_electronic_pc_correction = after_electronic_pc_correction
        self.before_ionic_pc_correction = before_ionic_pc_correction
        self.ave_pot_diff = ave_pot_diff

    @classmethod
    def from_scratch(cls,
                     structure: Structure,
                     before_charge: int,
                     after_charge: int,
                     unitcell_dft: UnitcellCalcResults,
                     efnv_correction: ExtendedFnvCorrection):

        charge_diff = after_charge - before_charge
        if abs(charge_diff) != 1:
            raise ValueError(f"Charge transition {before_charge} -> "
                             f"{after_charge} is invalid.")

        volume = structure.lattice.volume
        ave_pot_diff = efnv_correction.ave_pot_diff
        dielectric_ele = unitcell_dft.static_dielectric_tensor
        dielectric_ion = unitcell_dft.ionic_dielectric_tensor
        ewald_ele = Ewald.from_optimization(structure, dielectric_ele)
        ewald_ion = Ewald.from_optimization(structure, dielectric_ion)

        before_electronic_pc_correction = \
            point_charge_energy(before_charge, ewald_ele, volume)
        after_electronic_pc_correction = \
            point_charge_energy(after_charge, ewald_ele, volume)
        before_ionic_pc_correction = \
            point_charge_energy(before_charge, ewald_ion, volume)

        return cls(before_charge, after_charge, ewald_ele, ewald_ion,
                   before_electronic_pc_correction,
                   after_electronic_pc_correction,
                   before_ionic_pc_correction, ave_pot_diff)

    @property
    def correction_energy(self):
        if self.after_charge - self.before_charge == 1:
            return (self.after_electronic_pc_correction
                    - self.before_electronic_pc_correction
                    - self.before_ionic_pc_correction
                    - self.ave_pot_diff)
        elif self.after_charge - self.before_charge == -1:
            return (self.after_electronic_pc_correction
                    - self.before_electronic_pc_correction
                    + self.before_ionic_pc_correction
                    + self.ave_pot_diff)
        else:
            raise ValueError(f"Charge transition {self.before_charge} -> "
                             f"{self.after_charge} is invalid.")

    @property
    def pc_correction_energy(self):
        if self.after_charge - self.before_charge == 1:
            return (self.after_electronic_pc_correction
                    - self.before_electronic_pc_correction
                    - self.before_ionic_pc_correction)
        elif self.after_charge - self.before_charge == -1:
            return (self.after_electronic_pc_correction
                    - self.before_electronic_pc_correction
                    + self.before_ionic_pc_correction)
        else:
            raise ValueError(f"Charge transition {self.before_charge} -> "
                             f"{self.after_charge} is invalid.")
