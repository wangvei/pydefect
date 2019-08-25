from monty.json import MSONable
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.corrections.corrections import Correction
from pydefect.corrections.efnv_corrections import Ewald, ExtendedFnvCorrection
from pydefect.corrections.efnv_corrections import point_charge_energy
from pymatgen.core.structure import Structure

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
                 before_ele_pc: float,
                 after_ele_pc: float,
                 before_ion_pc: float,
                 ave_pot_diff: float):
        """
        Args:
            before_charge (int):
                Charge state before the transition.
            after_charge (int):
                Charge state after the transition.
            ewald_electronic_part (Ewald):
                Ewald object for the electronic part.
            ewald_ionic_part (Ewald):
                Ewald object for the ionic part.
            before_ele_pc (float):
                Point-charge (PC) correction energy for electronic part before
                the transition.
            after_ele_pc (float):
                Point-charge (PC) correction energy for electronic part after
                the transition.
            before_ion_pc (float):
                Point-charge (PC) correction energy for ionic part before the
                transition.
            ave_pot_diff (float)):
                Average potential difference before the transition.

        """
        if abs(after_charge - before_charge) != 1:
            raise ValueError(f"Charge transition {before_charge} -> "
                             f"{after_charge} is invalid.")

        self.before_charge = before_charge
        self.after_charge = after_charge
        self.ewald_electronic_part = ewald_electronic_part
        self.ewald_ionic_part = ewald_ionic_part
        self.before_ele_pc = before_ele_pc
        self.after_ele_pc = after_ele_pc
        self.before_ion_pc = before_ion_pc
        self.ave_pot_diff = ave_pot_diff

    @classmethod
    def from_scratch(cls,
                     structure: Structure,
                     before_charge: int,
                     after_charge: int,
                     unitcell_dft: UnitcellCalcResults,
                     efnv_correction: ExtendedFnvCorrection):

        volume = structure.lattice.volume
        ave_pot_diff = efnv_correction.ave_pot_diff
        dielectric_ele = unitcell_dft.static_dielectric_tensor
        dielectric_ion = unitcell_dft.ionic_dielectric_tensor
        ewald_ele = Ewald.from_optimization(structure, dielectric_ele)
        ewald_ion = Ewald.from_optimization(structure, dielectric_ion)

        print(dielectric_ele, dielectric_ion)

        before_ele_pc = point_charge_energy(before_charge, ewald_ele, volume)
        after_ele_pc = point_charge_energy(after_charge, ewald_ele, volume)
        before_ion_pc = -2 * point_charge_energy(before_charge, ewald_ion, volume) / before_charge

        return cls(before_charge, after_charge, ewald_ele, ewald_ion,
                   before_ele_pc, after_ele_pc, before_ion_pc, ave_pot_diff)

    def __repr__(self):
        outs = [f"Total correction: {self.correction_energy}",
                f"After electronic PC: {self.after_ele_pc}",
                f"Before electronic PC: {self.before_ele_pc}",
                f"Before ionic PC: {self.before_ion_pc}",
                f"Finite-size contribution: {self.ave_pot_diff}"]

        return "\n".join(outs)

    @property
    def correction_energy(self) -> float:
        charge_diff = self.after_charge - self.before_charge
        return self.pc_correction_energy + charge_diff * self.ave_pot_diff

    @property
    def pc_correction_energy(self) -> float:
        charge_diff = self.after_charge - self.before_charge
        return ((self.after_ele_pc - self.before_ele_pc) -
                charge_diff * self.before_ion_pc)

