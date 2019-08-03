# -*- coding: utf-8 -*-
from enum import Enum, unique
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.corrections.corrections import Correction
from pydefect.util.structure_tools import get_displacements
from pydefect.vasp_util.util import calc_orbital_similarity
from pymatgen.electronic_structure.core import Spin
from pydefect.util.logger import get_logger

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
                        dissimilarity_criterion: float = 0.12,
                        localized_criterion: float = 0.4,
                        near_edge_energy_criterion: float = 0.5):
    """
    Args:
        participation_ratio (dict):
            keys: [spin][band_edge], where band_edge is "hob" or "lub",
                  meaning "highest-occupied band", "lowest-unoccupied band".
            values: participation ratio at the neighboring sites
        orbital_character (dict):
            keys: [s][band_edge][position], where position is "max" or "min"
            values (dict): [element_name][orbital]
                           orbital = "s", "p", "d" or "f" (if exists)
        perfect_orbital_character (dict):
            Orbital character for perfect supercell for comparison.
        band_edge_energies:
            keys: [spin][band_edge]
            values (float): eigenvalue in absolute scale.
        supercell_vbm (float):
            VBM in the perfect supercell.
        supercell_cbm (float):
            CBM in the perfect supercell.
        dissimilarity_criterion:
            Determines whether the eigenstate is a host state.
        localized_criterion (float):
            Judge the orbital is localized if the participation ratio is
            larger than this value.
        near_edge_energy_criterion (float):
            Judge whether the orbital energy is close to the band edge.
    """
    if all([orbital_character, perfect_orbital_character]) is False:
        logger.warning("Diagnosing shallow states is impossible.")
        return False

    band_edges = {}

    for spin in orbital_character:
        # Consider the situation where perfect has only Spin.up.
        try:
            perfect = perfect_orbital_character[spin]
        except KeyError:
            perfect = perfect_orbital_character[Spin.up]

        # hob
        vbm_dissimilarity = calc_orbital_similarity(
            orbital_character[spin]["hob"]["max"], perfect["hob"]["max"])
        donor_phs_dissimilarity = calc_orbital_similarity(
            orbital_character[spin]["hob"]["min"], perfect["lub"]["min"])
        # lub
        cbm_dissimilarity = calc_orbital_similarity(
            orbital_character[spin]["lub"]["min"], perfect["lub"]["min"])
        acceptor_phs_dissimilarity = calc_orbital_similarity(
            orbital_character[spin]["lub"]["max"], perfect["hob"]["max"])

        if vbm_dissimilarity < dissimilarity_criterion \
                and cbm_dissimilarity < dissimilarity_criterion \
                and participation_ratio[spin]["hob"] < localized_criterion \
                and participation_ratio[spin]["lub"] < localized_criterion \
                and supercell_cbm - band_edge_energies[spin]["lub"] \
                < near_edge_energy_criterion \
                and band_edge_energies[spin]["hob"] - supercell_vbm \
                < near_edge_energy_criterion:
            band_edges[spin] = BandEdgeState.no_in_gap

        # To determine the dissimilarity of perturbed host state (phs),
        # the criterion is doubled from that of regular band edge.
        elif donor_phs_dissimilarity < dissimilarity_criterion * 2 \
                and participation_ratio[spin]["hob"] < localized_criterion \
                and supercell_cbm - band_edge_energies[spin]["hob"] \
                < near_edge_energy_criterion:
            band_edges[spin] = BandEdgeState.donor_phs

        elif acceptor_phs_dissimilarity < dissimilarity_criterion * 2 \
                and participation_ratio[spin]["lub"] < localized_criterion \
                and band_edge_energies[spin]["lub"] - supercell_vbm \
                < near_edge_energy_criterion:
            band_edges[spin] = BandEdgeState.acceptor_phs

        else:
            band_edges[spin] = BandEdgeState.localized_state

    return band_edges


class Defect:
    def __init__(self,
                 defect_entry: DefectEntry,
                 dft_results: SupercellCalcResults,
                 perfect_dft_results: SupercellCalcResults,
                 correction: Correction = None):

        self.de = defect_entry
        self.dr = dft_results
        self.pdr = perfect_dft_results
        self.c = correction

        self.relative_total_energy = self.dr.total_energy - self.pdr.total_energy

        self.relative_potential = []
        for d_atom, p_atom in enumerate(self.de.atom_mapping_to_perfect):
            if p_atom is None:
                self.relative_potential.append(None)
            else:
                pot_defect = self.dr.electrostatic_potential[d_atom]
                pot_perfect = self.pdr.electrostatic_potential[p_atom]
                self.relative_potential.append(pot_defect - pot_perfect)

        self.displacements = get_displacements(self.dr.final_structure,
                                               self.de.initial_structure,
                                               self.dr.defect_center,
                                               self.de.anchor_atom_index)
        self.band_edge_states = \
            diagnose_band_edges(self.dr.participation_ratio,
                                self.dr.orbital_character,
                                self.pdr.orbital_character,
                                self.dr.band_edge_energies,
                                self.dr.vbm, self.dr.cbm)


    @property
    def is_shallow(self) -> bool:
        return any(i.is_shallow for i in self.band_edge_states.values())
