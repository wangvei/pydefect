# -*- coding: utf-8 -*-
import json
from enum import Enum, unique

import numpy as np
from monty.json import MSONable, MontyEncoder
from monty.serialization import loadfn
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.defect_name import DefectName
from pydefect.core.error_classes import StructureError
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.corrections.corrections import Correction
from pydefect.database.symmetry import num_symmetry_operation as nsymop
from pydefect.database.atom import rcore
from pydefect.util.logger import get_logger
from pydefect.util.tools import spin_key_to_str, str_key_to_spin
from pydefect.util.vasp_util import calc_orbital_difference
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


def too_close_atom_pairs(structure: Structure,
                         radius: float = 2,
                         too_close_criterion_factor: float = 0.7) -> bool:
    """Check whether too close atomic pairs exist in the structure or not.

    Raise a ValueError if the number of atoms in structure is 1.


    Args:
        structure (Structure):
            Input structure


    """

    distances = structure.get_all_neighbors(radius)
    for i, periodic_neighbors in enumerate(distances):
        elem1 = structure[i].species_string
        frac1 = structure[i].frac_coords
        for periodic_neighbor in periodic_neighbors:
            elem2 = periodic_neighbor.species_string
            frac2 = periodic_neighbor.frac_coords
            dist = periodic_neighbor.nn_distance
            rcore_sum = rcore[elem1] + rcore[elem2]
            if dist < rcore_sum * too_close_criterion_factor:
                logger.warning(
                    f"Element {elem1} (index: {i}) at {frac1} & Element {elem2}"
                    f" (index: {periodic_neighbor.index}) at {frac2} are too "
                    f"close with distance {round(dist, 4)}.")

                return True
    return False


def diagnose_band_edges(participation_ratio: dict,
                        orbital_character: dict,
                        perfect_orbital_character: dict,
                        band_edge_energies: dict,
                        perfect_supercell_vbm: float,
                        perfect_supercell_cbm: float,
                        similarity_criterion: float = 0.12,
                        localized_criterion: float = 0.5,
                        near_edge_energy_criterion: float = 0.3):
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
        perfect_supercell_vbm (float):
            VBM in the perfect supercell.
        perfect_supercell_cbm (float):
            CBM in the perfect supercell.
        similarity_criterion:
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
        # Large this value means the occupied band should be the VBM
        vbm_orbital_difference = calc_orbital_difference(
            orbital_character[spin]["hob"]["top"], perfect["hob"]["top"])

        # lowest-unoccupied band = lub
        # Large this value means the occupied band should be the CBM
        cbm_orbital_difference = calc_orbital_difference(
            orbital_character[spin]["lub"]["bottom"], perfect["lub"]["bottom"])

        hob_bottom_from_cbm = \
            band_edge_energies[spin]["hob"]["bottom"] - perfect_supercell_cbm
        lub_top_from_vbm = \
            band_edge_energies[spin]["lub"]["top"] - perfect_supercell_vbm

        # When there is no any in-gap state,
        # 1. Vbm and cbm orbitals are similar to the perfect ones.
        # 2. Participation rations are small enough.
        # 3. Band-edge energies should be close to the host ones.
        # Note that the participation ratio depends on the supercell size.
        if (vbm_orbital_difference < similarity_criterion
                and cbm_orbital_difference < similarity_criterion
                and participation_ratio[spin]["hob"] < localized_criterion
                and participation_ratio[spin]["lub"] < localized_criterion
                and abs(band_edge_energies[spin]["lub"]["bottom"]
                        - perfect_supercell_cbm) < near_edge_energy_criterion
                and abs(band_edge_energies[spin]["hob"]["top"]
                        - perfect_supercell_vbm) < near_edge_energy_criterion):
            band_edges[spin] = BandEdgeState.no_in_gap

        # When the highest-occupied band (=hob) is a perturbed host state,
        # 1. Participation rations are small enough.
        # 2. Energy of bottom of hob is close to or higher than supercell cbm.
        elif (participation_ratio[spin]["hob"] < localized_criterion
              and hob_bottom_from_cbm + near_edge_energy_criterion > 0) \
                or hob_bottom_from_cbm > 0:
            band_edges[spin] = BandEdgeState.donor_phs

        elif (participation_ratio[spin]["lub"] < localized_criterion
              and lub_top_from_vbm - near_edge_energy_criterion < 0)  \
                or lub_top_from_vbm < 0:
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
                 initial_multiplicity: int,
                 final_multiplicity: int,
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
        # multiplicity
        self.initial_multiplicity = initial_multiplicity
        self.final_multiplicity = final_multiplicity
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
                     correction: Correction = None,
                     check_structure: bool = True,
                     raise_error: bool = True):
        """ Gather and generate full defect information related to analysis.

        While SupercellCalcResults class only collect the data from the
        first-principles calculation results, this method also extracts useful
        information from the data, e.g., band_edge_states.
        """
        relative_total_energy = \
            dft_results.total_energy - perfect_dft_results.total_energy

        initial_nsymop = nsymop(defect_entry.initial_site_symmetry)
        final_nsymop = nsymop(dft_results.site_symmetry)
        initial_multiplicity = defect_entry.multiplicity
        final_multiplicity = \
            initial_multiplicity * initial_nsymop / final_nsymop

        if not final_multiplicity.is_integer() and raise_error:
            n = DefectName(
                defect_entry.name, defect_entry.charge, defect_entry.annotation)
            raise ValueError(
                f"Multiplicity of {n} is invalid. "
                f"initial symmetry: {defect_entry.initial_site_symmetry}, "
                f"final symmetry: {dft_results.site_symmetry}.")

        magnetization = round(dft_results.total_magnetization, 2)

        if check_structure:
            if too_close_atom_pairs(dft_results.final_structure):
                raise StructureError(
                    "Some atoms are too close. If you want to continue to "
                    "generate Defect object, set check_structure=False.")

        band_edge_states = \
            diagnose_band_edges(dft_results.participation_ratio,
                                dft_results.orbital_character,
                                perfect_dft_results.orbital_character,
                                dft_results.band_edge_energies,
                                perfect_dft_results.vbm,
                                perfect_dft_results.cbm)

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
                   initial_multiplicity=initial_multiplicity,
                   final_multiplicity=final_multiplicity,
                   magnetization=magnetization,
                   kpoint_coords=dft_results.kpoint_coords,
                   eigenvalues=dft_results.eigenvalues,
                   supercell_vbm=perfect_dft_results.vbm,
                   supercell_cbm=perfect_dft_results.cbm,
                   fermi_level=dft_results.fermi_level,
                   orbital_character=dft_results.orbital_character,
                   band_edge_states=band_edge_states)

    @classmethod
    def from_dict(cls, d):
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for s, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(s))] = np.array(v)

        initial_structure = d["initial_structure"]
        perturbed_initial_structure = d["perturbed_initial_structure"]
        final_structure = d["final_structure"]

        if isinstance(initial_structure, dict):
            initial_structure = Structure.from_dict(initial_structure)
        if isinstance(perturbed_initial_structure, dict):
            perturbed_initial_structure = \
                Structure.from_dict(perturbed_initial_structure)
        if isinstance(final_structure, dict):
            final_structure = Structure.from_dict(final_structure)

        orbital_character = str_key_to_spin(d["orbital_character"])
        band_edge_states = \
            str_key_to_spin(d["band_edge_states"], BandEdgeState.from_string)

        return cls(name=d["name"],
                   charge=d["charge"],
                   annotation=d["annotation"],
                   is_converged=d["is_converged"],
                   relative_total_energy=d["relative_total_energy"],
                   correction_energy=d["correction_energy"],
                   initial_structure=initial_structure,
                   perturbed_initial_structure=perturbed_initial_structure,
                   final_structure=final_structure,
                   final_volume=d["final_volume"],
                   defect_center=d["defect_center"],
                   anchor_atom_index=d["anchor_atom_index"],
                   changes_of_num_elements=d["changes_of_num_elements"],
                   displacements=d["displacements"],
                   neighboring_sites=d["neighboring_sites"],
                   defect_center_coords=d["defect_center_coords"],
                   initial_symmetry=d["initial_symmetry"],
                   final_symmetry=d["final_symmetry"],
                   initial_multiplicity=d["initial_multiplicity"],
                   final_multiplicity=d["final_multiplicity"],
                   magnetization=d["magnetization"],
                   kpoint_coords=d["kpoint_coords"],
                   eigenvalues=eigenvalues,
                   supercell_vbm=d["supercell_vbm"],
                   supercell_cbm=d["supercell_cbm"],
                   fermi_level=d["fermi_level"],
                   orbital_character=orbital_character,
                   band_edge_states=band_edge_states)

    def as_dict(self):
        eigenvalues = \
            {str(spin): v.tolist() for spin, v in self.eigenvalues.items()}
        orbital_character = spin_key_to_str(self.orbital_character)
        band_edge_states = \
            spin_key_to_str(self.band_edge_states, value_to_str=True)

        d = {"@module":                     self.__class__.__module__,
             "@class":                      self.__class__.__name__,
             "name":                        self.name,
             "charge":                      self.charge,
             "annotation":                  self.annotation,
             "is_converged":                self.is_converged,
             "relative_total_energy":       self.relative_total_energy,
             "correction_energy":           self.correction_energy,
             "initial_structure":           self.initial_structure,
             "perturbed_initial_structure": self.perturbed_initial_structure,
             "final_structure":             self.final_structure,
             "final_volume":                self.final_volume,
             "defect_center":               self.defect_center,
             "anchor_atom_index":           self.anchor_atom_index,
             "changes_of_num_elements":     self.changes_of_num_elements,
             "displacements":               self.displacements,
             "neighboring_sites":           self.neighboring_sites,
             "defect_center_coords":        self.defect_center_coords,
             "initial_symmetry":            self.initial_symmetry,
             "final_symmetry":              self.final_symmetry,
             "initial_multiplicity":        self.initial_multiplicity,
             "final_multiplicity":          self.final_multiplicity,
             "magnetization":               self.magnetization,
             "kpoint_coords":               self.kpoint_coords,
             "eigenvalues":                 eigenvalues,
             "supercell_vbm":               self.supercell_vbm,
             "supercell_cbm":               self.supercell_cbm,
             "fermi_level":                 self.fermi_level,
             "orbital_character":           orbital_character,
             "band_edge_states":            band_edge_states}

        return d

    def to_json_file(self, filename="defect.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename="defect.json"):
        return loadfn(filename)

    @property
    def diagnose(self) -> str:
        band_edges = []
        for s, v in self.band_edge_states.items():
            band_edges.extend([s.name.upper(), "->", str(v).rjust(17)])

        outs = [f" Convergence: {self.is_converged}",
                f" Magmom: {self.magnetization:>7}",
                f" Point group: {self.final_symmetry:>4}",
                f" Band edge: {'  '.join(band_edges)}"]

        return " ".join(outs)

    @property
    def is_shallow(self) -> bool:
        return any(i.is_shallow for i in self.band_edge_states.values())

    def set_band_edge_state(self, spin: Spin, state: str):
        self.band_edge_states[spin] = BandEdgeState.from_string(state)
