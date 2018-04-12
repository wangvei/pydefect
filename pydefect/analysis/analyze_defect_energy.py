# -*- coding: utf-8 -*-

# from pydefect.core.defect_set import DefectSet


class EnergyPlotter:


    def __init__(self, unitcell, chem_pot, perfect, defects):
        """

        is_lower_energy: If there are some defects with the same name and
        charge, use the lowest energy. False, raise error.

        """

        self._unitcell = unitcell
        self._chem_pot = chem_pot
        self._perfect = perfect
        self._defects = defects

    @classmethod
    def from_directories(cls, unitcell, chem_pot, perfect, defects):
        """
        """


