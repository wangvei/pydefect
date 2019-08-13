from monty.json import MontyEncoder, MSONable
from pydefect.corrections.corrections import Correction
from pydefect.corrections.efnv_corrections import Ewald, ExtendedFnvCorrection
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
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
                 ewald_electronic_part,
                 ewald_ionic_part,
                 before_electronic_pc_correction: float,
                 after_electronic_pc_correction: float,
                 before_ionic_pc_correction: float,
                 alignmentlike_correction: float):

        self.ewald_electronic_part = ewald_electronic_part
        self.ewald_ionic_part = ewald_ionic_part
        self.before_electronic_pc_correction = before_electronic_pc_correction
        self.after_electronic_pc_correction = after_electronic_pc_correction
        self.before_ionic_pc_correction = before_ionic_pc_correction
        self.alignmentlike_correction = alignmentlike_correction

    @classmethod
    def from_scratch(cls,
                     structure: Structure,
                     before_charge: int,
                     after_charge: int,
                     unitcell_dft: UnitcellCalcResults,
                     efnv_correction: ExtendedFnvCorrection):

        alignmentlike_correction = efnv_correction.alignment_correction_energy
        diele_ele = unitcell_dft.static_dielectric_tensor
        diele_ion = unitcell_dft.ionic_dielectric_tensor
        ewald_ele = Ewald.from_optimization(structure, diele_ele)
        ewald_ion = Ewald.from_optimization(structure, diele_ion)


