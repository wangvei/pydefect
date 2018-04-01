
import json
from monty.json import MontyEncoder
from monty.serialization import loadfn

from pydefect.core.supercell import Perfect

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "February 19, 2018"


class Unitcell(Perfect):
    """
    This class object holds some properties related to a unit cell calculation.
    Args:
        dft_results (UnitcellDftResults): SupercellDftResults class object
        defect_entry (DefectEntry): DefectEntry class object
            |- initial_structure
            |- removed_atoms
            |- inserted_atoms
            |- changes_of_num_elements
            |- in_name
            |- out_name
    """

    @property
    def static_dielectric_tensor(self):
        return self._dft_results._static_dielectric_tensor

    @property
    def ionic_dielectric_tensor(self):
        return self._dft_results._ionic_dielectric_tensor

    @property
    def total_dielectric_tensor(self):
        return self.static_dielectric_tensor + self.ionic_dielectric_tensor


    # def __init__(self, dft_results=None):
    #     self._dft_results = dft_results

    # def as_dict(self):
    #     """
    #     Dictionary representation of Defect class object.
    #     """
    #     d = {"dft_results":  self._dft_results.as_dict()}

        # return d

    # @classmethod
    # def from_dict(cls, d):
    #     """
    #     Constructs a class object from a dictionary.
    #     """
    #     dft_results = UnitcellDftResults.from_dict(d["dft_results"])

        # return cls(dft_results)

    # @property
    # def initial_structure(self):
    #     return self._defect_entry.initial_structure

    # @property
    # def removed_atoms(self):
    #     return self._defect_entry.removed_atoms

    # @property
    # def inserted_atoms(self):
    #     return self._defect_entry.inserted_atoms

    # @property
    # def changes_of_num_elements(self):
    #     return self._defect_entry.changes_of_num_elements

    # @property
    # def charge(self):
    #     return self._defect_entry.charge

    # @property
    # def in_name(self):
    #     return self._defect_entry.in_name

    # @property
    # def out_name(self):
    #     return self._defect_entry.out_name

    # @property
    # def relative_total_energy(self):
    #     if self._relative_total_energy:
    #         return self._relative_total_energy
    #     else:
    #         warnings.warn("relative_total_energy is not set.")
    #         return None

    # @property
    # def relative_potential(self):
    #     if self._relative_potential:
    #         return self._relative_potential
    #     else:
    #         warnings.warn("relative_potential is not set.")
    #         return None

#     #    @property
#     #    def relative_atomic_displacements(self):
#     #        if self._relative_atomic_displacements:
#     #            return self._relative_atomic_displacements
#     #        else:
#     #            warnings.warn("relative_atomic_displacements is not set.")
#     #            return None

    # def set_relative_values(self, perfect):
    #     """
    #     Set relative values with respect to those of perfect calc.
    #     Args:
    #         perfect (Perfect):
    #     """
    #     self._relative_total_energy = self.total_energy - perfect.total_energy
    #     self._relative_potential = list(map(operator.sub,
    #                                         self.electrostatic_potential,
    #                                         perfect.electrostatic_potential))
    #     # TODO: write below
    # #        self._relative_atomic_displacements = \
    # pass

