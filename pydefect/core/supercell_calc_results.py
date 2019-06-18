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

from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.surface import get_recp_symmetry_operation
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun, Procar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.error_classes import NoConvergenceError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import get_displacements
from pydefect.vasp_util.util import calc_participation_ratio, \
    calc_orbital_character, calc_orbital_similarity

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


@unique
class BandEdges(Enum):
    donor_phs = "Donor PHS"
    acceptor_phs = "Acceptor PHS"
    localized_state = "Localized state"
    no_in_gap = "No in-gap state"

    def __str__(self):
        return self.value

    @classmethod
    def from_string(cls, s):
        for m in cls:
            if m.value == s or m.name == s:
                return m
        raise AttributeError("Band edge info: " + str(s) + " is not proper.\n" +
                             "Supported info:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.name for e in cls])

    @property
    def is_shallow(self):
        if self in [BandEdges.acceptor_phs, BandEdges.donor_phs]:
            return True
        else:
            return False


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
                 eigenvalue_properties,
                 volume: float,
                 fermi_level: float,
                 is_converged: bool,
                 band_edges: dict = None,
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
            band_edges (dict):
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
        self.eigenvalue_properties = eigenvalue_properties
        self.volume = volume
        self.fermi_level = fermi_level
        self.is_converged = is_converged
        self.band_edges = band_edges
        self.band_edge_energies = band_edge_energies
        self.relative_total_energy = relative_total_energy
        self.relative_potential = relative_potential
        self.displacements = displacements
        self.symmetrized_structure = symmetrized_structure
        self.symmops = symmops
        self.participation_ratio = participation_ratio
        self.orbital_character = orbital_character

    def __str__(self):
        outs = ["total energy (eV): " + str(self.total_energy),
                "total total_magnetization (mu_B): " +
                str(self.total_magnetization),
                "electrostatic potential: " + str(self.electrostatic_potential),
                "eigenvalues_properties: " + str(self.eigenvalue_properties),
                "final structure: \n" + str(self.final_structure),
                "site_symmetry: " + str(self.site_symmetry),
                "volume: \n" + str(self.volume),
                "Fermi level (eV): \n" + str(self.fermi_level),
                "Orbital character (eV): \n" + str(self.orbital_character)]

        return "\n".join(outs)

    @property
    def diagnose(self):
        band_edges = []
        for s, v in self.band_edges.items():
            band_edges.extend([s.name.upper(), ":", str(v).rjust(17)])

        band_edges = "  ".join(band_edges)

        outs = ["convergence : " + str(self.is_converged)[0],
                "  band edge : " + band_edges]

        print("  ".join(outs))

    @classmethod
    def from_vasp_files(cls,
                        directory_path: str,
                        vasprun: str = None,
                        contcar: str = None,
                        outcar: str = None,
                        procar: Union[str, bool] = False,
                        referenced_dft_results=None,
                        defect_entry: DefectEntry = None,
                        symmetrize: bool = True,
                        symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                        angle_tolerance: float = ANGLE_TOL):
        """ Constructs class object from vasp output files.

        Args:
            directory_path (str):
                path to the directory storing calc results.
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
            defect_entry (DefectEntry):
            symmetrize (bool):
                Whether to obtain the symmetrized information.
            symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
        """
        p = Path(directory_path)

        # get the names of the latest files in the directory_path
        if vasprun is None:
            vasprun = str(max(p.glob("**/*vasprun*"), key=os.path.getctime))
        if contcar is None:
            contcar = \
                str(max(list(p.glob("**/*CONTCAR*")) +
                        list(p.glob("**/*POSCAR*")), key=os.path.getctime))
        if outcar is None:
            outcar = str(max(p.glob("**/*OUTCAR*"), key=os.path.getctime))

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
        # (band gap, cbm, vbm, is_band_gap_direct)
        eigenvalue_properties = vasprun.eigenvalue_band_properties
        fermi_level = vasprun.efermi

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
            symmetrized_structure = sga.get_symmetrized_structure()
            symmops = get_recp_symmetry_operation(symmetrized_structure)
        else:
            site_symmetry = None
            symmetrized_structure = None
            symmops = None

        outcar = parse_file(Outcar, outcar)
        total_energy = outcar.final_energy
        magnetization = 0.0 if outcar.total_mag is None else outcar.total_mag
        electrostatic_potential = outcar.electrostatic_potential

        if procar:
            if procar is True:
                procar = str(max(p.glob("**/*PROCAR*"), key=os.path.getctime))
            procar = parse_file(Procar, procar)

            neighboring_sites = \
                defect_entry.neighboring_sites if defect_entry else None

            # The k-point indices at the band edges in defect calculations.
            # hob (lub) means highest (lowest) (un)occupied state
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

                        if neighboring_sites:
                            participation_ratio[s][band_edge] = \
                                calc_participation_ratio(
                                    procar=procar,
                                    spin=s,
                                    band_index=band_index,
                                    kpoint_index=k,
                                    atom_indices=neighboring_sites)

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

        # Perfect supercell does not need the below.
        relative_total_energy = None
        relative_potential = None
        displacements = None
        band_edges = None

        if referenced_dft_results:
            relative_total_energy = \
                total_energy - referenced_dft_results.total_energy

            if defect_entry is None:
                raise ValueError("DefectEntry is necessary for analyzing"
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
                get_displacements(final_structure, initial_structure,
                                  defect_entry.defect_center,
                                  defect_entry.anchor_atom_index)

            if participation_ratio and orbital_character:
                perfect_orbital_character = \
                    referenced_dft_results.orbital_character

                supercell_vbm = referenced_dft_results.eigenvalue_properties[2]
                supercell_cbm = referenced_dft_results.eigenvalue_properties[1]

                band_edges = \
                    cls.diagnose_band_edges(participation_ratio,
                                            orbital_character,
                                            perfect_orbital_character,
                                            band_edge_energies,
                                            supercell_vbm,
                                            supercell_cbm)

        return cls(final_structure=final_structure,
                   site_symmetry=site_symmetry,
                   total_energy=total_energy,
                   total_magnetization=magnetization,
                   eigenvalues=eigenvalues,
                   kpoint_coords=kpoint_coords,
                   kpoint_weights=kpoint_weights,
                   electrostatic_potential=electrostatic_potential,
                   eigenvalue_properties=eigenvalue_properties,
                   volume=volume,
                   fermi_level=fermi_level,
                   is_converged=vasprun.converged_ionic,
                   band_edges=band_edges,
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
            orbital_character (dict):
            perfect_orbital_character (dict):
                Orbital character for perfect supercell for comparison.
            dissimilarity_criterion:
                Determines whether the eigenstate is a host state.
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

            # print("vbm_dissimilarity, donor_phs_dissimilarity, "
            #       "cbm_dissimilarity, acceptor_phs_dissimilarity")
            # print(vbm_dissimilarity, donor_phs_dissimilarity, cbm_dissimilarity,
            #       acceptor_phs_dissimilarity)
            # print("vbm_participation_ratio", "cbm_participation_ratio")
            # print(participation_ratio[spin]["hob"],
            #       participation_ratio[spin]["lub"])

            if vbm_dissimilarity < dissimilarity_criterion \
                    and cbm_dissimilarity < dissimilarity_criterion \
                    and participation_ratio[spin]["hob"] < localized_criterion \
                    and participation_ratio[spin]["lub"] < localized_criterion \
                    and supercell_cbm - band_edge_energies[spin]["lub"] \
                    < near_edge_energy_criterion \
                    and band_edge_energies[spin]["hob"] - supercell_vbm \
                    < near_edge_energy_criterion:
                band_edges[spin] = BandEdges.no_in_gap

            elif donor_phs_dissimilarity < dissimilarity_criterion * 2 \
                    and participation_ratio[spin]["hob"] < localized_criterion \
                    and supercell_cbm - band_edge_energies[spin]["hob"] \
                    < near_edge_energy_criterion:
                band_edges[spin] = BandEdges.donor_phs

            elif acceptor_phs_dissimilarity < dissimilarity_criterion * 2 \
                    and participation_ratio[spin]["lub"] < localized_criterion \
                    and band_edge_energies[spin]["lub"] - supercell_vbm \
                    < near_edge_energy_criterion:
                band_edges[spin] = BandEdges.acceptor_phs

            else:
                band_edges[spin] = BandEdges.localized_state

        return band_edges

    def set_band_edges(self, spin: Spin, state: Union[str, None]):
        """ Set band edges manually."""
        self.band_edges[spin] = BandEdges.from_string(state)
        return True

    @classmethod
    def from_dict(cls, d):
        """ Construct a class object from a dictionary. """

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
                d = {}
                for spin, value in arg.items():
                    if value_to_band_edges:
                        d[Spin(int(spin))] = BandEdges.from_string(value)
                    else:
                        d[Spin(int(spin))] = value
                return d
            else:
                return

        band_edges = str_key_to_spin(d["band_edges"], value_to_band_edges=True)
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
                   eigenvalue_properties=d["eigenvalue_properties"],
                   volume=d["volume"],
                   fermi_level=d["fermi_level"],
                   is_converged=d["is_converged"],
                   band_edges=band_edges,
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

        band_edges = spin_key_to_str(self.band_edges, value_to_str=True)
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
             "eigenvalue_properties":   self.eigenvalue_properties,
             "volume":                  self.volume,
             "fermi_level":             self.fermi_level,
             "is_converged":            self.is_converged,
             "band_edges":              band_edges,
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

