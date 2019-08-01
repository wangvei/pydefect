# -*- coding: utf-8 -*-

from collections import defaultdict
from enum import Enum, unique
import json
import numpy as np
import os
from typing import Union

from pathlib import Path

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from xml.etree.ElementTree import ParseError

from obadb.analyzer.band_gap import band_gap_properties

from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.surface import get_recp_symmetry_operation
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun, Procar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL, \
    CUTOFF_RADIUS
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.error_classes import NoConvergenceError, StructureError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import get_displacements
from pydefect.vasp_util.util import calc_participation_ratio, \
    calc_orbital_character, calc_orbital_similarity

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


class SupercellCalcResults(MSONable):
    """ Class with DFT results for supercell systems. """

    def __init__(self,
                 final_structure: Structure,
                 site_symmetry: str,
                 total_energy: float,
                 total_magnetization: float,
                 eigenvalues: np.array,
                 kpoint_coords: list,
                 kpoint_weights: list,
                 electrostatic_potential: list,
                 vbm: dict,
                 cbm: dict,
                 volume: float,
                 fermi_level: float,
                 is_converged: bool,
                 defect_center: Union[int, list],
                 neighboring_sites_after_relax: list,
                 band_edge_states: Union[dict, None] = None,
                 band_edge_energies: dict = None,
                 relative_total_energy: float = None,
                 relative_potential: list = None,
                 displacements: list = None,
                 symmetrized_structure: Structure = None,
                 symmops: list = None,
                 participation_ratio: dict = None,
                 orbital_character: dict = None):
        """
        None is set for some properties of perfect supercell.

        Args:
            final_structure (Structure):
                pmg Structure class object. Usually relaxed structures
            site_symmetry (str):
                Site symmetry after structure optimization.
            total_energy (float):
                Final total energy in eV.
            total_magnetization (float):
                Total total_magnetization in mu_B
            eigenvalues (N_spin x N_kpoint x N_band np.array):
                Numpy array for the electron eigenvalues in absolute scale.
            kpoint_coords (list):
                List of k-point coordinates
            kpoint_weights (list):
                List of k-point weights.
            electrostatic_potential (list):
                Atomic site electrostatic potential.
            volume (float):
                Volume of the supercell.
            fermi_level (float):
               Fermi level in the absolute scale.
            is_converged (bool):
                Whether the calculation is converged or not.
            defect_center (int/list):
                Show the defect center. When the defect center is an
                atomic position, atom index number is set. Contrary, when
                the defect center is defined by the initial position, such as
                vacancies and complex defects, the fractional coordinates are
                set.
            neighboring_sites_after_relax (list):
                Atomic indices of the neighboring sites within cutoff radius
                after structure optimization.
            band_edge_states (dict):
                Band edge states at each spin channel.
                None: no in gap state.
                "Donor PHS": Donor-type perturbed host state (PHS).
                "Acceptor PHS": Acceptor-type PHS.
                "Localized state": With in-gap localized state.
                    ex. {Spin.up: None, Spin.down:"Localized state}
            band_edge_energies (dict):
            symmetrized_structure (Structure):
                Symmetrized structure with a defect.
            symmops (list):
                Point-group symmetry operation.
            participation_ratio (dict):
                ex. {Spin.up: {"hob": 0.72, "lub": 0.34,
                     Spin.down: {"hob": 0.32, "lub": 0.14}}
            orbital_character (dict):
                ex. {Spin.up: {"hob": {"Mg": {"s": 0.1, ...}, "O": {...},
                               "lub": {...},
                     Spin.down: {...}}
        """
        self.final_structure = final_structure
        self.site_symmetry = site_symmetry
        self.total_energy = total_energy
        self.total_magnetization = total_magnetization
        self.eigenvalues = eigenvalues
        self.kpoint_coords = kpoint_coords
        self.kpoint_weights = kpoint_weights
        self.electrostatic_potential = electrostatic_potential
        self.vbm = vbm
        self.cbm = cbm
        self.volume = volume
        self.fermi_level = fermi_level
        self.is_converged = is_converged
        self.defect_center = defect_center
        self.neighboring_sites_after_relax = neighboring_sites_after_relax
        self.band_edge_states = band_edge_states
        self.band_edge_energies = band_edge_energies
        self.relative_total_energy = relative_total_energy
        self.relative_potential = relative_potential
        self.displacements = displacements
        self.symmetrized_structure = symmetrized_structure
        self.symmops = symmops
        self.participation_ratio = participation_ratio
        self.orbital_character = orbital_character

    def __repr__(self):
        outs = ["total energy (eV): " + str(self.total_energy),
                "total total_magnetization (mu_B): "
                + str(self.total_magnetization),
                "electrostatic potential: "
                + str(self.electrostatic_potential),
                "vbm: " + str(self.vbm),
                "cbm: " + str(self.cbm),
                "final structure: \n" + str(self.final_structure),
                "site_symmetry: " + str(self.site_symmetry),
                "volume: \n" + str(self.volume),
                "Fermi level (eV): \n" + str(self.fermi_level),
                "Orbital character (eV): \n" + str(self.orbital_character)]

        return "\n".join(outs)

    @property
    def diagnose(self) -> str:

        band_edges = []
        for s, v in self.band_edge_states.items():
            band_edges.extend([s.name.upper(), "->", str(v).rjust(17)])

        band_edges = "  ".join(band_edges)

        outs = [" conv. : " + str(self.is_converged)[0],
                "  mag. : {}".format(str(round(self.total_magnetization, 1))),
                "  sym. : {:>4}".format(self.site_symmetry),
                "  band edge : " + band_edges]

        return " ".join(outs)

    @classmethod
    def from_vasp_files(cls,
                        directory_path: str,
                        defect_entry: DefectEntry = None,
                        vasprun: str = None,
                        contcar: str = None,
                        outcar: str = None,
                        procar: Union[str, bool] = False,
                        referenced_dft_results=None,
                        symmetrize: bool = True,
                        cutoff: float = CUTOFF_RADIUS,
                        symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                        angle_tolerance: float = ANGLE_TOL):
        """ Constructs class object from vasp output files.

        Args:
            directory_path (str):
                path to the directory storing calc results.
            defect_entry (DefectEntry):
            vasprun (str):
                Name of the vasprun.xml file.
            contcar (str):
                Name of the converged CONTCAR file.
            outcar (str):
                Name of the OUTCAR file.
            procar (str):
                Name of the PROCAR file.
                If True, parse the PROCAR file but file name is determined from
                timestamp.
            referenced_dft_results (SupercellCalcResults):
                SupercellDftResults object for referenced supercell dft results.
            symmetrize (bool):
                Whether to obtain the symmetrized information.
            cutoff (float):
                Cutoff Radius determining the neighboring sites.
            symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
        """
        p = Path(directory_path)

        # get the names of the latest files in the directory_path
        if vasprun is None:
            vasprun = str(max(p.glob("*vasprun*"), key=os.path.getctime))
        if contcar is None:
            contcar = \
                str(max(list(p.glob("*CONTCAR*")) +
                        list(p.glob("*POSCAR*")), key=os.path.getctime))
        if outcar is None:
            outcar = str(max(p.glob("*OUTCAR*"), key=os.path.getctime))

        def parse_file(classmethod_name, parsed_filename):
            try:
                logger.info("Parsing {}...".format(parsed_filename))
                return classmethod_name(parsed_filename)
            except ParseError:
                logger.warning("Parsing {} failed.".format(parsed_filename))
                raise ParseError
            except FileNotFoundError:
                logger.warning("File {} doesn't exist.".format(parsed_filename))
                raise FileNotFoundError

        vasprun = parse_file(Vasprun, vasprun)
        eigenvalues = vasprun.eigenvalues
        fermi_level = vasprun.efermi
        gap_properties = band_gap_properties(vasprun)

        if gap_properties:
            vbm = gap_properties[1]["energy"]
            cbm = gap_properties[2]["energy"]
        else:
            vbm = fermi_level
            cbm = fermi_level

        # Check if the electronic and ionic steps are converged.
        if vasprun.converged_electronic is False:
            raise NoConvergenceError("Electronic step is not converged.")

        if vasprun.converged_ionic is False:
            logger.warning("Ionic step is not converged.")

        kpoint_coords = vasprun.actual_kpoints
        kpoint_weights = vasprun.actual_kpoints_weights

        contcar = parse_file(Poscar.from_file, contcar)
        final_structure = contcar.structure
        volume = contcar.structure.volume

        if symmetrize:
            sga = SpacegroupAnalyzer(final_structure, symprec, angle_tolerance)
            site_symmetry = sga.get_point_group_symbol()
            symmetrized_structure = sga.get_refined_structure()
            symmops = get_recp_symmetry_operation(symmetrized_structure)
        else:
            site_symmetry = None
            symmetrized_structure = None
            symmops = None

        outcar = parse_file(Outcar, outcar)
        total_energy = outcar.final_energy
        magnetization = 0.0 if outcar.total_mag is None else outcar.total_mag
        electrostatic_potential = outcar.electrostatic_potential

        neighboring_sites_after_relax = []
        if defect_entry:
            if defect_entry.defect_type.is_defect_center_atom:
                defect_center = list(defect_entry.inserted_atoms.keys())[0]
                defect_coords = final_structure[defect_center].frac_coords
            else:
                defect_center = defect_coords = \
                    defect_entry.defect_center_coords
            for i, site in enumerate(final_structure):
                distance = \
                    site.distance_and_image_from_frac_coords(defect_coords)[0]
                # Defect itself is not included to the neighboring sites as
                # neighboring_sites property in DefectEntry.
                if 1e-5 < distance < cutoff:
                    neighboring_sites_after_relax.append(i)

            if not neighboring_sites_after_relax:
                raise StructureError(f"No neighboring site detected. Increase "
                                     f"cutoff radius from {cutoff}.")

        else:
            defect_center = None

        # Perfect supercell does not need the below.
        relative_total_energy = None
        relative_potential = None
        displacements = None
        band_edge_states = None

        if procar:
            if procar is True:
                procar = str(max(p.glob("*PROCAR*"), key=os.path.getctime))
            procar = parse_file(Procar, procar)
            # The k-point indices at the band edges in defect calculations.
            # hob (lub) means highest (lowest) (un)occupied state
            # The tiny number needs to be added to avoid magnetization = 0.0
            if magnetization == 0.0:
                magnetization = 0.001
            hob_index = \
                {Spin.up: round((outcar.nelect + magnetization) / 2) - 1,
                 Spin.down: round((outcar.nelect - magnetization) / 2) - 1}

            band_edge_energies = defaultdict(dict)
            participation_ratio = defaultdict(dict)
            orbital_character = defaultdict(dict)

            for s in eigenvalues.keys():
                orbital_character[s] = {}
                for i, band_edge in enumerate(["hob", "lub"]):

                    orbital_character[s][band_edge] = {}

                    band_index = hob_index[s] + i
                    band_edge_energies[s][band_edge] = \
                        eigenvalues[s][0, band_index, 0]

                    max_eigenvalue = np.amax(eigenvalues[s][:, band_index, 0])
                    max_k = int(np.where(eigenvalues[s][:, band_index, 0]
                                         == max_eigenvalue)[0][0])

                    min_eigenvalue = np.amax(eigenvalues[s][:, band_index, 0])
                    min_k = int(np.where(eigenvalues[s][:, band_index, 0]
                                         == min_eigenvalue)[0][0])

                    for position, k in zip(["max", "min"], [max_k, min_k]):

                        participation_ratio[s][band_edge] = \
                            calc_participation_ratio(
                                procar=procar,
                                spin=s,
                                band_index=band_index,
                                kpoint_index=k,
                                atom_indices=neighboring_sites_after_relax)

                        orbital_character[s][band_edge][position] = \
                            calc_orbital_character(
                                procar=procar,
                                structure=final_structure,
                                spin=s,
                                band_index=band_index,
                                kpoint_index=k)

            if participation_ratio:
                participation_ratio = dict(participation_ratio)
            else:
                participation_ratio = None
            orbital_character = dict(orbital_character)
            band_edge_energies = dict(band_edge_energies)

        else:
            participation_ratio = None
            orbital_character = None
            band_edge_energies = None

        if referenced_dft_results:
            relative_total_energy = \
                total_energy - referenced_dft_results.total_energy

            if defect_entry is None:
                raise ValueError("DefectEntry is necessary for analyzing "
                                 "relative values.")

            mapping = defect_entry.atom_mapping_to_perfect
            relative_potential = []

            for d_atom, p_atom in enumerate(mapping):

                if p_atom is None:
                    relative_potential.append(None)
                else:
                    potential_defect = electrostatic_potential[d_atom]
                    potential_perfect = \
                        referenced_dft_results.electrostatic_potential[p_atom]
                    relative_potential.append(
                        potential_defect - potential_perfect)

            initial_structure = defect_entry.initial_structure

            displacements = \
                get_displacements(final_structure,
                                  initial_structure,
                                  defect_center,
                                  defect_entry.anchor_atom_index)

            if participation_ratio and orbital_character:
                perfect_orbital_character = \
                    referenced_dft_results.orbital_character

                band_edge_states = \
                    cls.diagnose_band_edges(participation_ratio,
                                            orbital_character,
                                            perfect_orbital_character,
                                            band_edge_energies,
                                            vbm, cbm)

        return cls(final_structure=final_structure,
                   site_symmetry=site_symmetry,
                   total_energy=total_energy,
                   total_magnetization=magnetization,
                   eigenvalues=eigenvalues,
                   kpoint_coords=kpoint_coords,
                   kpoint_weights=kpoint_weights,
                   electrostatic_potential=electrostatic_potential,
                   vbm=vbm,
                   cbm=cbm,
                   volume=volume,
                   fermi_level=fermi_level,
                   is_converged=vasprun.converged_ionic,
                   defect_center=defect_center,
                   neighboring_sites_after_relax=neighboring_sites_after_relax,
                   band_edge_states=band_edge_states,
                   band_edge_energies=band_edge_energies,
                   relative_total_energy=relative_total_energy,
                   relative_potential=relative_potential,
                   displacements=displacements,
                   symmetrized_structure=symmetrized_structure,
                   symmops=symmops,
                   participation_ratio=participation_ratio,
                   orbital_character=orbital_character)

    @staticmethod
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

    def set_band_edges(self, spin: Spin, state: Union[str, None]):
        self.band_edge_states[spin] = BandEdgeState.from_string(state)
        return True

    @classmethod
    def from_dict(cls, d):
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for s, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(s))] = np.array(v)

        final_structure = d["final_structure"]
        if isinstance(final_structure, dict):
            final_structure = Structure.from_dict(final_structure)

        symmops = []
        for symmop in d["symmops"]:
            if isinstance(symmop, dict):
                symmops.append(SymmOp.from_dict(symmop))
            else:
                symmops.append(symmop)

        def str_key_to_spin(arg: dict, value_to_band_edges=False):
            if arg is not None:
                x = {}
                for spin, value in arg.items():
                    if value_to_band_edges:
                        x[Spin(int(spin))] = BandEdgeState.from_string(value)
                    else:
                        x[Spin(int(spin))] = value
                return x
            else:
                return

        band_edge_states = \
            str_key_to_spin(d["band_edge_states"], value_to_band_edges=True)
        band_edge_energies = str_key_to_spin(d["band_edge_energies"])
        participation_ratio = str_key_to_spin(d["participation_ratio"])
        orbital_character = str_key_to_spin(d["orbital_character"])

        return cls(final_structure=final_structure,
                   site_symmetry=d["site_symmetry"],
                   total_energy=d["total_energy"],
                   total_magnetization=d["total_magnetization"],
                   eigenvalues=eigenvalues,
                   kpoint_coords=d["kpoint_coords"],
                   kpoint_weights=d["kpoint_weights"],
                   electrostatic_potential=d["electrostatic_potential"],
                   vbm=d["vbm"],
                   cbm=d["cbm"],
                   volume=d["volume"],
                   fermi_level=d["fermi_level"],
                   is_converged=d["is_converged"],
                   defect_center=d["defect_center"],
                   neighboring_sites_after_relax
                   =d["neighboring_sites_after_relax"],
                   band_edge_states=band_edge_states,
                   band_edge_energies=band_edge_energies,
                   relative_total_energy=d["relative_total_energy"],
                   relative_potential=d["relative_potential"],
                   displacements=d["displacements"],
                   symmetrized_structure=d["symmetrized_structure"],
                   symmops=symmops,
                   participation_ratio=participation_ratio,
                   orbital_character=orbital_character)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    def as_dict(self):
        # Spin object must be converted to string for to_json_file.
        eigenvalues = \
            {str(spin): v.tolist() for spin, v in self.eigenvalues.items()}

        def spin_key_to_str(arg: dict, value_to_str=False):
            if arg is not None:
                if value_to_str:
                    return {str(spin): str(v) for spin, v in arg.items()}
                else:
                    return {str(spin): v for spin, v in arg.items()}
            else:
                return

        band_edge_states = \
            spin_key_to_str(self.band_edge_states, value_to_str=True)
        band_edge_energies = spin_key_to_str(self.band_edge_energies)
        participation_ratio = spin_key_to_str(self.participation_ratio)
        orbital_character = spin_key_to_str(self.orbital_character)

        d = {"@module":                 self.__class__.__module__,
             "@class":                  self.__class__.__name__,
             "final_structure":         self.final_structure,
             "site_symmetry":           self.site_symmetry,
             "total_energy":            self.total_energy,
             "total_magnetization":     self.total_magnetization,
             "eigenvalues":             eigenvalues,
             "kpoint_coords":           self.kpoint_coords,
             "kpoint_weights":          self.kpoint_weights,
             "electrostatic_potential": self.electrostatic_potential,
             "vbm":                     self.vbm,
             "cbm":                     self.cbm,
             "volume":                  self.volume,
             "fermi_level":             self.fermi_level,
             "is_converged":            self.is_converged,
             "defect_center":           self.defect_center,
             "neighboring_sites_after_relax":
                 self.neighboring_sites_after_relax,
             "band_edge_states":        band_edge_states,
             "band_edge_energies":      band_edge_energies,
             "relative_total_energy":   self.relative_total_energy,
             "relative_potential":      self.relative_potential,
             "displacements":           self.displacements,
             "symmetrized_structure":   self.symmetrized_structure,
             "symmops":                 self.symmops,
             "participation_ratio":     participation_ratio,
             "orbital_character":       orbital_character}

        return d

    def to_json_file(self, filename):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def relative_total_energy(self, referenced_dft_results):
        """ Return a relative total energy w.r.t. referenced supercell.

        Args:
            referenced_dft_results (SupercellCalcResults):
                SupercellDftResults object for referenced supercell dft results.
        """
        return self.total_energy - referenced_dft_results.total_energy

    def relative_potential(self, referenced_dft_results, defect_entry):
        """
        Return a list of relative site potential w.r.t. the perfect supercell.
        None is inserted for interstitial sites.

        Args:
            referenced_dft_results (SupercellCalcResults):
                SupercellDftResults object for referenced supercell dft results.
                Usually it is for perfect supercell.
            defect_entry (DefectEntry):
                DefectEntry class object corresponding the SupercellDftResuts.
        """
        mapping = defect_entry.atom_mapping_to_perfect
        relative_potential = []

        for d_atom, p_atom in enumerate(mapping):

            if p_atom is None:
                relative_potential.append(None)
            else:
                potential_defect = self.electrostatic_potential[d_atom]
                potential_perfect = \
                    referenced_dft_results.electrostatic_potential[p_atom]
                relative_potential.append(potential_defect - potential_perfect)

        return relative_potential

#    def inserted_atom_displacements(self, defect_entry):
#        """
#        Returns coordinates of defect center by calculating the averaged
#        coordinates. If len(defect_coords) == 1, returns defect_coords[0].
#        Args:
#            defect_entry (DefectEntry):
#                related DefectEntry class object
#        """
#        displacements = []
#
#        for k in defect_entry.inserted_atoms.keys:
#            before_relaxation = defect_entry.initial_structure.frac_coords[k]
#            after_relaxation = self.final_structure.frac_coords[k]
#            displacements.append(
#                min_distance_and_its_v2coord(before_relaxation,
#                                             after_relaxation,
#                                             self.final_structure.axis))
#        return displacements

