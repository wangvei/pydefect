# -*- coding: utf-8 -*-

# from pydefect.core.defect_set import DefectSet


class EnergyPlotter:


    def __init__(self, unitcell, chem_pot, perfect, defects):
        """

        is_lower_energy: If there are some defects with the same name and
        charge, use the lowest energy. False, raise error.

        1. Set relative energies of defects.
        2. Construct a set of data used for calculating the defect energies.
           {(name1, charge1): relative_energy, element_diff, correction_energy
            (name2, charge2): relative_energy, element_diff, correction_energy
            ....
        3. Calculate the charge transition level.
        4. Plot

        We can obtain the calculated spurious VBM and CBM from perfect

        """

        self._unitcell = unitcell
        self._chem_pot = chem_pot
        self._perfect = perfect
        self._defects = defects

    @classmethod
    def from_files(cls,
                   unitcell_json,
                   chem_pot_json):


