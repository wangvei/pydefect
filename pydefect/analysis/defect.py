# -*- coding: utf-8 -*-
from enum import Enum, unique

import numpy as np
from monty.json import MontyEncoder, MSONable
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.corrections.corrections import Correction
from pydefect.database.num_symmetry_operation \
    import num_symmetry_operation as nsymop
from pydefect.util.logger import get_logger
from pydefect.vasp_util.util import calc_orbital_difference
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


@unique
class BandEdgeState(Enum):
    donor_phs = "Donor PHS"
    acceptor_phs = "Acceptor PHS"
    localized_state = "Localized state"
    no_in_gap = "No in-gap state"

    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value

    @classmethod
    def from_string(cls, s):
        for m in cls:
            if m.value == s or m.name == s:
                return m
        raise AttributeError(f"Band edge info {str(s)}  is not proper.\n",
                             f"Supported info:\n {cls.name_list()}")

    @classmethod
    def name_list(cls):
        return ', '.join([e.name for e in cls])

    @property
    def is_shallow(self):
        return self in [BandEdgeState.acceptor_phs, BandEdgeState.donor_phs]


def diagnose_band_edges(participation_ratio: dict,
                        orbital_character: dict,
                        perfect_orbital_character: dict,
                        band_edge_energies: dict,
                        supercell_vbm: float,
                        supercell_cbm: float,
                        different_criterion: float = 0.12,
                        localized_criterion: float = 0.4,
                        near_edge_energy_criterion: float = 0.5):
    """ Diagnose the band edge states in supercell with a defect.

    Args:
        participation_ratio (dict):
            keys: [spin][band_edge], where band_edge is "hob" or "lub",
                  meaning "highest-occupied band", "lowest-unoccupied band".
            values: participation ratio at the neighboring sites
        orbital_character (dict):
            keys: [spin][band_edge][position] (position = "top" or "bottom")
            values (dict): [element_name][orbital]
                           orbital = "s", "p", "d" or "f" (if exists)
        perfect_orbital_character (dict):
            Orbital character for perfect supercell for comparison.
        band_edge_energies:
            Averaged band energy over k-space as function of spin and band_edge
            keys: [spin][band_edge]
            values (float): eigenvalue in absolute scale.
        supercell_vbm (float):
            VBM in the perfect supercell.
        supercell_cbm (float):
            CBM in the perfect supercell.
        different_criterion:
            Criterion to judge if the eigenstate is different from a host state
        localized_criterion (float):
            Criterion to judge the orbital is localized if the participation
            ratio is larger than this value.
        near_edge_energy_criterion (float):
            Criterion to judge if orbital energy is close to band edge.
    """
    band_edges = {}

    for spin in orbital_character:
        # Spin.up state is substituted for Spin.down one if the latter is
        # absent in perfect
        try:
            perfect = perfect_orbital_character[spin]
        except KeyError:
            perfect = perfect_orbital_character[Spin.up]

        # highest-occupied band = hob
        vbm_difference = calc_orbital_difference(
            orbital_character[spin]["hob"]["top"], perfect["hob"]["top"])
        # The difference between the bottom of the highest-occupied
        # band and the conduction band minimum (= bottom of the
        # lowest-unoccupied band in the perfect supercell)
        donor_phs_difference = calc_orbital_difference(
            orbital_character[spin]["hob"]["bottom"], perfect["lub"]["bottom"])

        # lowest-unoccupied band = lub
        cbm_difference = calc_orbital_difference(
            orbital_character[spin]["lub"]["bottom"], perfect["lub"]["bottom"])
        acceptor_phs_difference = calc_orbital_difference(
            orbital_character[spin]["lub"]["top"], perfect["hob"]["top"])

        # When there is no any in-gap state,
        # 1. Vbm and cbm orbitals are similar to the perfect ones.
        # 2. Participation rations are small enough.
        # 3. Band-edge energies should be close to the host ones.
        # Note that the participation ratio depends on the supercell size.
        if (vbm_difference < different_criterion
                and cbm_difference < different_criterion
                and participation_ratio[spin]["hob"] < localized_criterion
                and participation_ratio[spin]["lub"] < localized_criterion
                and abs(band_edge_energies[spin]["lub"]["bottom"]
                        - supercell_cbm) < near_edge_energy_criterion
                and abs(band_edge_energies[spin]["hob"]["top"]
                        - supercell_vbm) < near_edge_energy_criterion):
            band_edges[spin] = BandEdgeState.no_in_gap

        # When the highest-occupied band (=hob) is a perturbed host state,
        # 1. Bottom of hob is similar to the cbm in perfect.
        # 2. Participation rations are small enough.
        # 3. Energy of bottom of hob is close to the supercell cbm.

        # To determine the difference of perturbed host state (phs),
        # the criterion is doubled from that of regular band edge, since
        # the
        elif (donor_phs_difference < different_criterion * 2
                and participation_ratio[spin]["hob"] < localized_criterion
                and abs(band_edge_energies[spin]["hob"]["bottom"]
                        - supercell_cbm) < near_edge_energy_criterion):
            band_edges[spin] = BandEdgeState.donor_phs

        elif (acceptor_phs_difference < different_criterion * 2
                and participation_ratio[spin]["lub"] < localized_criterion
                and abs(band_edge_energies[spin]["lub"]["top"] -
                        supercell_vbm) < near_edge_energy_criterion):
            band_edges[spin] = BandEdgeState.acceptor_phs

        else:
            band_edges[spin] = BandEdgeState.localized_state

    return band_edges


class Defect(MSONable):
    """ Holds full defect information related to analysis. """
    def __init__(self,
                 name: str,
                 charge: int,
                 annotation: str,
                 is_converged: bool,
                 relative_total_energy: float,
                 correction_energy: float,
                 initial_structure: Structure,
                 perturbed_initial_structure: Structure,
                 final_structure: Structure,
                 defect_center: list,
                 anchor_atom_index: int,
                 changes_of_num_elements: dict,
                 displacements: dict,
                 neighboring_sites: list,
                 defect_center_coords: list,
                 initial_symmetry: str,
                 final_symmetry: str,
                 final_volume: float,
                 num_equiv_sites: int,
                 multiplicity: int,
                 magnetization: float,
                 kpoint_coords: list,
                 eigenvalues: np.array,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 fermi_level: float,
                 orbital_character: dict,
                 band_edge_states: dict):
        """ See each object docstrings for details."""

        self.name = name
        self.charge = charge
        self.annotation = annotation
        self.is_converged = is_converged
        # energies
        self.relative_total_energy = relative_total_energy
        self.correction_energy = correction_energy
        # structures
        self.initial_structure = initial_structure
        self.perturbed_initial_structure = perturbed_initial_structure
        self.final_structure = final_structure
        self.final_volume = final_volume
        self.defect_center = defect_center
        self.anchor_atom_index = anchor_atom_index
        self.changes_of_num_elements = changes_of_num_elements
        self.displacements = displacements
        self.neighboring_sites = neighboring_sites
        self.defect_center_coords = defect_center_coords
        # symmetries
        self.initial_symmetry = initial_symmetry
        self.final_symmetry = final_symmetry
        self.num_equiv_sites = num_equiv_sites
        self.multiplicity = multiplicity
        # magnetization
        self.magnetization = magnetization
        # kpoint
        self.kpoint_coords = kpoint_coords
        # eigenvalues and eigenstates
        self.eigenvalues = eigenvalues
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.fermi_level = fermi_level
        self.orbital_character = orbital_character
        self.band_edge_states = band_edge_states

    @classmethod
    def from_objects(cls,
                     defect_entry: DefectEntry,
                     dft_results: SupercellCalcResults,
                     perfect_dft_results: SupercellCalcResults,
                     correction: Correction = None):
        """ Gather and generate full defect information related to analysis.

        While SupercellCalcResults class only collect the data from the
        first-principles calculation results, this method also extracts useful
        information from the data, e.g., band_edge_states.
        """
        relative_total_energy = \
            dft_results.total_energy - perfect_dft_results.total_energy

        initial_nsymop = nsymop(defect_entry.initial_site_symmetry)
        final_nsymop = nsymop(dft_results.site_symmetry)
        num_equiv_sites = defect_entry.num_equiv_sites
        multiplicity = num_equiv_sites * initial_nsymop / final_nsymop

        magnetization = dft_results.total_magnetization

        band_edge_states = \
            diagnose_band_edges(dft_results.participation_ratio,
                                dft_results.orbital_character,
                                perfect_dft_results.orbital_character,
                                dft_results.band_edge_energies,
                                dft_results.vbm, dft_results.cbm)

        if correction:
            correction_energy = correction.correction_energy
        else:
            correction_energy = None

        return cls(name=defect_entry.name,
                   charge=defect_entry.charge,
                   annotation=defect_entry.annotation,
                   is_converged=dft_results.is_converged,
                   relative_total_energy=relative_total_energy,
                   correction_energy=correction_energy,
                   initial_structure=defect_entry.initial_structure,
                   perturbed_initial_structure=
                   defect_entry.perturbed_initial_structure,
                   final_structure=dft_results.final_structure,
                   final_volume=dft_results.volume,
                   defect_center=dft_results.defect_center,
                   anchor_atom_index=defect_entry.anchor_atom_index,
                   changes_of_num_elements=
                   defect_entry.changes_of_num_elements,
                   displacements=dft_results.displacements,
                   neighboring_sites=dft_results.neighboring_sites_after_relax,
                   defect_center_coords=dft_results.defect_coords,
                   initial_symmetry=defect_entry.initial_site_symmetry,
                   final_symmetry=dft_results.site_symmetry,
                   num_equiv_sites=num_equiv_sites,
                   multiplicity=multiplicity,
                   magnetization=magnetization,
                   kpoint_coords=dft_results.kpoint_coords,
                   eigenvalues=dft_results.eigenvalues,
                   supercell_vbm=perfect_dft_results.vbm,
                   supercell_cbm=perfect_dft_results.cbm,
                   fermi_level=dft_results.fermi_level,
                   orbital_character=dft_results.orbital_character,
                   band_edge_states=band_edge_states)

#    def as_dict(self):


    @property
    def diagnose(self) -> str:
        band_edges = []
        for s, v in self.band_edge_states.items():
            band_edges.extend([s.name.upper(), "->", str(v).rjust(17)])

        outs = [f" conv. : {self.is_converged}",
                f"  mag. : {self.magnetization}",
                f"  sym. : {self.final_symmetry:>4}",
                f"  band edge : {'  '.join(band_edges)}"]

        return " ".join(outs)

    @property
    def is_shallow(self) -> bool:
        return any(i.is_shallow for i in self.band_edge_states.values())
