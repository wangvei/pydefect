# -*- coding: utf-8 -*-

from monty.json import MSONable
from pydefect.core.error_classes import NoConvergenceError

from typing import Union, Optional

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.analysis.defect_eigenvalues import DefectEigenvalue
from pydefect.analysis.defect_structure import DefectStructure
from pydefect.analysis.defect_energy import DefectEnergy
from pydefect.corrections.corrections import Correction

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class DefectProperties(MSONable):
    """
    """

    def __init__(self,
                 name: str,
                 is_converged: bool,
                 defect_eigenvalue: DefectEigenvalue,
                 defect_energy: DefectEnergy,
                 defect_structure: DefectStructure,
                 corrections: list,
                 is_shallow: Optional[bool] = None,):
        """
        Args:
            name (str):
                Name of a defect
            is_converged (bool):
                Whether converged or not.
            defect_eigenvalue (DefectEigenvalue):
            defect_energy (DefectEnergy):
            defect_structure (DefectStructure):
            is_shallow (bool):
                Whether is shallow or not.
        """
        self.name = name
        self.is_converged = is_converged
        self.defect_eigenvalue = defect_eigenvalue
        self.defect_energy = defect_energy
        self.defect_structure = defect_structure
        self._corrections = corrections
        self.is_shallow = is_shallow

    @classmethod
    def from_calc_results(cls,
                          defect_entry: Union[str, DefectEntry],
                          unitcell_results: Union[str, UnitcellCalcResults],
                          perfect_results: Union[str, SupercellCalcResults],
                          defect_results: Union[str, SupercellCalcResults],
                          find_deep_states: bool = True,
                          tolerance: float = DEFECT_SYMMETRY_TOLERANCE,
                          angle_tolerance: float = ANGLE_TOL):

        if isinstance(perfect_results, str):
            perfect_results = SupercellCalcResults.load_json(perfect_results)

        if not perfect_results.is_converged:
            raise NoConvergenceError("Perfect supercell is not converged.")

        if isinstance(defect_entry, str):
            defect_entry = DefectEntry.load_json(defect_entry)
        if isinstance(unitcell_results, str):
            unitcell_results = UnitcellCalcResults.load_json(unitcell_results)
        if isinstance(defect_results, str):
            defect_results = SupercellCalcResults.load_json(defect_results)

        name = defect_entry.name
        is_converged = defect_results.is_converged

        kpoints = defect_results.kpoints
        eigenvalues = defect_results.eigenvalues
        vbm, cbm = unitcell_results.band_edge
        band_gap = cbm - vbm
        _, supercell_cbm, supercell_vbm, _ = \
            perfect_results.eigenvalue_properties
        fermi_level = defect_results.fermi_level
        total_magnetization = defect_results.total_magnetization
        eigenvalue_correction = {}
        deep_states = []

        defect_eigenvalue = \
            DefectEigenvalue(name, kpoints, eigenvalues, vbm, cbm, band_gap,
                             supercell_vbm, supercell_cbm, fermi_level,
                             total_magnetization, eigenvalue_correction,
                             deep_states)

        relative_total_energy = \
            defect_results.total_energy - perfect_results.total_energy

        # get relative_potential
        mapping = defect_entry.atom_mapping_to_perfect
        relative_potential = []
        for d_atom, p_atom in enumerate(mapping):
            if p_atom is None:
                relative_potential.append(None)
            else:
                d_potential = defect_results.electrostatic_potential[d_atom]
                p_potential = perfect_results.electrostatic_potential[p_atom]
                relative_potential.append(d_potential - p_potential)

        changes_of_num_elements = defect_entry.changes_of_num_elements
        charge = defect_entry.charge
        energy_correction = {}

        defect_energy = \
            DefectEnergy(name, relative_total_energy, relative_potential,
                         changes_of_num_elements, charge, energy_correction)

        volume = defect_results.volume
        initial_structure = defect_entry.initial_structure
        perturbed_initial_structure = defect_entry.perturbed_initial_structure
        final_structure = defect_results.final_structure
        initial_site_symmetry = defect_entry.initial_site_symmetry
        initial_symmetry_multiplicity = \
            defect_entry.initial_symmetry_multiplicity
        perturbed_sites = defect_entry.perturbed_sites
        num_equiv_sites = defect_entry.num_equiv_sites

        sga = SpacegroupAnalyzer(final_structure, tolerance, angle_tolerance)
        final_site_symmetry = sga.get_point_group_symbol()
        final_symmetry_multiplicity = len(sga.get_point_group_operations())

        from pydefect.util.structure import get_displacements
        initial_distances, final_distances, displacements, angles = \
            get_displacements(final_structure, initial_structure, center, anchor_atom_index)

        defect_structure = \
            DefectStructure(name, volume, initial_structure,
                            perturbed_initial_structure, final_structure,
                            displacements, initial_site_symmetry,
                            final_site_symmetry, initial_symmetry_multiplicity,
                            final_symmetry_multiplicity, perturbed_sites,
                            num_equiv_sites)

#        is_shallow = defect_eigenvalue.is_shallow

        return cls(name, is_converged, defect_eigenvalue, defect_energy,
                   defect_structure, is_shallow)

    @property
    def corrections(self):
        return self._corrections

    @property.setter
    def corrections(self, corrections):
        if isinstance(corrections, list):
            self._corrections.extend(list(corrections))
        elif isinstance(corrections, Correction):
            self._corrections.append(corrections)


