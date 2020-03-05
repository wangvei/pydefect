# -*- coding: utf-8 -*-

import json
from pathlib import Path
from pprint import pformat
from typing import Union, Optional

import numpy as np

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from pydefect.core.config import (
    DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL, CUTOFF_FACTOR)
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.error_classes import NoConvergenceError, StructureError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import (
    get_displacements, min_distance_from_coords)
from pydefect.util.tools import (
    spin_key_to_str, str_key_to_spin, defaultdict_to_dict,
    mod_defaultdict)
from pydefect.util.vasp_util import (
    calc_participation_ratio, calc_orbital_character)

from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun, Procar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vise.analyzer.band_gap import band_gap_properties
from vise.util.tools import parse_file

logger = get_logger(__name__)


class ProcarDefectProperty(MSONable):
    """ Class with DFT results for supercell systems. """

    def __init__(self,
                 band_edge_energies: dict,
                 orbital_character: dict,
                 orbital_character_indices: dict,
                 participation_ratio: dict):
        """
        Args:
            band_edge_energies (dict):
                Averaged band energy over k-space as functions of spin and
                band_edge.
            orbital_character (dict):
                Orbital character at the eigenstate of each spin, band_edge
                (="hob" or "lub"), and energy_position = (="top" or "bottom")
                ex. {Spin.up: {"hob": {"top": {"Mg": {"s": 0.1, ...},
                                               "O": {...},..
                               "lub": {...},
                     Spin.down: {...}}
            orbital_character_indices (dict):
                Band indices used for orbital_character.
            participation_ratio (dict):
               Participation ratio averaged over k-space as function of Spin and
               band_edge determined from the neighboring sites.
        """
        self.band_edge_energies = band_edge_energies
        self.orbital_character = orbital_character
        self.orbital_character_indices = orbital_character_indices
        self.participation_ratio = participation_ratio

    @classmethod
    def analyze_procar(cls,
                       hob_index: dict,
                       procar: Procar,
                       eigenvalues: dict,
                       structure: Structure,
                       neighboring_sites: list
                       ) -> "ProcarDefectProperty":
        """ Analyze Procar to investigate defect properties

        Args:
            hob_index (dict):
               Highest occupied band (HOB) index for each spin channel
               {Spin.up: 100, Spin.down: 100}.
               Note that the index of the lowest unoccupied band (LUB) is
               incremented from HOB index by one.
            procar (Procar):
                Pymatgen Procar class object.
            eigenvalues:
               eigenvalues[Spin][k-index][band-index] = [energy, occupation]
            structure:
                Structure used for extracting symbol_set
            neighboring_sites:
                Atomic site indices neighboring a defect determining the
                participation ratio. When None, participation_ratio is None.
        """

        band_edge_energies = mod_defaultdict(depth=3)
        orbital_character = mod_defaultdict(depth=3)
        orbital_character_indices = mod_defaultdict(depth=3)
        participation_ratio = mod_defaultdict(depth=2)

        for spin in eigenvalues.keys():
            # index i is used to increment band index from hob to lub
            for i, band_edge in enumerate(["hob", "lub"]):

                # The band index of LUB is incremented from HOB by 1.
                band_index = hob_index[spin] + i

                # participation_ratio is not calculated for perfect supercell.
                if neighboring_sites:
                    participation_ratio[spin][band_edge] = \
                        calc_participation_ratio(
                            procar, spin, band_index, neighboring_sites)

                top_eigenvalue = np.amax(eigenvalues[spin][:, band_index, 0])
                band_edge_energies[spin][band_edge]["top"] = top_eigenvalue
                top_k_index = int(
                    np.where(eigenvalues[spin][:, band_index, 0]
                             == top_eigenvalue)[0][0])

                bottom_eigenvalue = np.amin(eigenvalues[spin][:, band_index, 0])
                band_edge_energies[spin][band_edge]["bottom"] \
                    = bottom_eigenvalue
                bottom_k_index = int(
                    np.where(eigenvalues[spin][:, band_index, 0]
                             == bottom_eigenvalue)[0][0])

                for energy_position, k_index in \
                        zip(["top", "bottom"], [top_k_index, bottom_k_index]):

                    orbital_character_indices[spin][band_edge][energy_position] = \
                        {"band_index": band_index, "k_index": k_index}

                    orbital_character[spin][band_edge][energy_position] = \
                        calc_orbital_character(procar=procar,
                                               structure=structure,
                                               spin=spin,
                                               band_index=band_index,
                                               kpoint_index=k_index)

        # participation_ratio must be None for perfect supercell.
        if participation_ratio:
            participation_ratio = defaultdict_to_dict(participation_ratio)
        else:
            participation_ratio = None

        return cls(defaultdict_to_dict(band_edge_energies),
                   defaultdict_to_dict(orbital_character),
                   defaultdict_to_dict(orbital_character_indices),
                   participation_ratio)


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
                 defect_coords: list,
                 displacements: dict,
                 neighboring_sites_after_relax: list,
                 band_edge_energies: dict = None,
                 orbital_character: dict = None,
                 orbital_character_indices: dict = None,
                 participation_ratio: dict = None):
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
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
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
            defect_coords (list):
                Defect center coordinates.
            displacements (dict):
                See get_displacements function docstring.
            neighboring_sites_after_relax (list):
                Atomic indices of the neighboring sites within cutoff radius
                in the final structure.
            band_edge_energies (dict):
            orbital_character (dict):
            orbital_character_indices (dict):
            participation_ratio (dict):
                See analyze_procar function for details
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
        self.defect_coords = defect_coords
        self.displacements = displacements
        self.neighboring_sites_after_relax = neighboring_sites_after_relax
        self.band_edge_energies = band_edge_energies
        self.orbital_character = orbital_character
        self.orbital_character_indices = orbital_character_indices
        self.participation_ratio = participation_ratio

#        self.manually_set_defect_center = None

    def __repr__(self):
        outs = [f"total energy (eV): {self.total_energy}",
                f"total total_magnetization (mu_B): {self.total_magnetization}",
                f"VBM: {self.vbm}",
                f"CBM: {self.cbm}",
                f"site_symmetry: {self.site_symmetry}",
                f"volume: {self.volume}",
                f"Fermi level (eV): {self.fermi_level}",
                f"Orbital character indices:",
                pformat(self.orbital_character_indices),
                f"Orbital character :",
                pformat(self.orbital_character)]

        return "\n".join(outs)

    @classmethod
    def from_vasp_files(cls,
                        directory_path: str,
                        defect_entry: DefectEntry = None,
                        vasprun: str = "vasprun.xml",
                        contcar: str = "CONTCAR",
                        outcar: str = "OUTCAR",
                        procar: Optional[str] = "PROCAR",
                        cutoff: Optional[float] = None,
                        defect_symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                        angle_tolerance: float = ANGLE_TOL):
        """ Constructs class object from vasp output files.

        Note1: Do not analyze data here, but in Defect.from_object classmethod.
        Note2: Do not change "defect_symprec" to "symprec" as the attribute
               defaults are used in main.py

        Args:
            directory_path (str):
                path to the directory storing calc results.
            vasprun (str):
                Name of the vasprun.xml file.
            contcar (str):
                Name of the converged CONTCAR file.
            outcar (str):
                Name of the OUTCAR file.
            defect_entry (DefectEntry):
            procar (str):
                Name of the PROCAR file.
                If True, parse the PROCAR file but file name is determined from
                timestamp.
            cutoff (float):
                Cutoff Radius determining the neighboring sites.
            defect_symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
        """
        # vasprun related
        vasprun = parse_file(Vasprun, Path(directory_path) / vasprun)
        if vasprun.converged_electronic is False:
            raise NoConvergenceError("Electronic step is not converged.")
        if vasprun.converged_ionic is False:
            logger.warning("Ionic step is not converged.")

        # outcar related
        outcar = parse_file(Outcar, Path(directory_path) / outcar)
        magnetization = outcar.total_mag or 0.0

        try:
            _, vbm, cbm = \
                (i["energy"] for i in band_gap_properties(vasprun, outcar))
        except TypeError:
            vbm = cbm = vasprun.efermi

        # contcar related
        contcar = parse_file(Poscar.from_file, Path(directory_path) / contcar)
        final_structure = contcar.structure
        volume = contcar.structure.volume
        sga = SpacegroupAnalyzer(
            final_structure, defect_symprec, angle_tolerance)
        site_symmetry = sga.get_point_group_symbol()

        # If defect_entry is None, system is regarded as perfect supercell.
        if not defect_entry:
            center = neighboring_sites = defect_coords = displacements = None
        else:
            if defect_entry.defect_type.is_defect_center_atom:
                center = defect_entry.inserted_atoms[0]["index"]
                defect_coords = list(final_structure[center].frac_coords)
            else:
                center = defect_coords = defect_entry.defect_center_coords

            if cutoff is None:
                min_d = min_distance_from_coords(final_structure, defect_coords)
                cutoff = round(min_d * CUTOFF_FACTOR, 2)

            neighboring_sites = []
            for i, site in enumerate(final_structure):
                # Calculate the distance between defect and site
                d, _ = site.distance_and_image_from_frac_coords(defect_coords)
                # Defect itself is excluded from neighboring sites in Pydefect.
                if 1e-5 < d < cutoff:
                    neighboring_sites.append(i)
            if not neighboring_sites:
                raise StructureError(f"No neighboring site detected. Increase "
                                     f"cutoff radius from {cutoff}.")

            displacements = get_displacements(final_structure,
                                              defect_entry.initial_structure,
                                              center,
                                              defect_entry.anchor_atom_index)
        # procar related
        kwargs = {}
        if procar:
            procar = parse_file(Procar, Path(directory_path) / procar)
            # The k-point indices at the band edges in defect calculations.
            # hob (lub) = highest (lowest) (un)occupied state
            # The small number (0.1) must be added to avoid magnetization = 0.0
            up_index = round((outcar.nelect + (magnetization + 0.1)) / 2) - 1
            down_index = round((outcar.nelect - (magnetization + 0.1)) / 2) - 1
            hob_index = {Spin.up: up_index, Spin.down: down_index}

            prop = ProcarDefectProperty. \
                analyze_procar(hob_index, procar, vasprun.eigenvalues,
                               final_structure, neighboring_sites)
            for k, v in prop.as_dict().items():
                if k[0] != "@":
                    kwargs[k] = v

        return cls(final_structure=final_structure,
                   site_symmetry=site_symmetry,
                   total_energy=outcar.final_energy,
                   total_magnetization=magnetization,
                   eigenvalues=vasprun.eigenvalues,
                   kpoint_coords=vasprun.actual_kpoints,
                   kpoint_weights=vasprun.actual_kpoints_weights,
                   electrostatic_potential=outcar.electrostatic_potential,
                   vbm=vbm,
                   cbm=cbm,
                   volume=volume,
                   fermi_level=vasprun.efermi,
                   is_converged=vasprun.converged_ionic,
                   defect_center=center,
                   defect_coords=defect_coords,
                   displacements=displacements,
                   neighboring_sites_after_relax=neighboring_sites,
                   **kwargs)

    @classmethod
    def from_dict(cls, d):
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for s, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(s))] = np.array(v)

        final_structure = d["final_structure"]
        if isinstance(final_structure, dict):
            final_structure = Structure.from_dict(final_structure)

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
                   defect_coords=d["defect_coords"],
                   displacements=d["displacements"],
                   neighboring_sites_after_relax
                   =d["neighboring_sites_after_relax"],
                   band_edge_energies=band_edge_energies,
                   participation_ratio=participation_ratio,
                   orbital_character=orbital_character)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    def as_dict(self):
        # Spin object must be converted to string for to_json_file as Enum
        # is not compatible with MSONable.
        eigenvalues = \
            {str(spin): v.tolist() for spin, v in self.eigenvalues.items()}

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
             "defect_coords":           self.defect_coords,
             "displacements":           self.displacements,
             "neighboring_sites_after_relax":
                 self.neighboring_sites_after_relax,
             "band_edge_energies":      band_edge_energies,
             "orbital_character":       orbital_character,
             "participation_ratio":     participation_ratio}

        return d

    def to_json_file(self, filename):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

