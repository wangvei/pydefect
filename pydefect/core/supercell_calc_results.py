# -*- coding: utf-8 -*-

import json
import numpy as np
import os

from collections import defaultdict
from pathlib import Path

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from xml.etree.cElementTree import ParseError

from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.surface import get_recp_symmetry_operation
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun, Procar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.defect import DefectEntry
from pydefect.core.error_classes import NoConvergenceError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import get_displacements
from pydefect.vasp_util.util import calc_participation_ratio, \
    calc_orbital_character

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

MAGNETIZATION_THRESHOLD = 0.1


class SupercellCalcResults(MSONable):
    """ Container class with DFT results for supercell systems. """

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
                 shallow: bool = None,
                 relative_total_energy: float = None,
                 relative_potential: list = None,
                 displacements: list = None,
                 symmetrized_structure: Structure = None,
                 symmops: list = None,
                 participation_ratio: dict = None,
                 orbital_character: dict = None):
        """
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
            symmops (list):
                Point-group symmetry operation.
            shallow (list):
                If we don't know Whether the defect is shallow, None.
                If the defect is detected not to be shallow, [].
                If the defect is detected to be shallow, store the detected ways
                like, ["not_integer_mag", "fermi_is_higher_than_cbm"]
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
        self.shallow = shallow
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
                "Fermi level (eV): \n" + str(self.fermi_level)]

        return "\n".join(outs)

    @classmethod
    def from_vasp_files(cls,
                        directory_path: str,
                        vasprun: str = None,
                        contcar: str = None,
                        outcar: str = None,
                        procar: str = None,
                        parse_procar: bool = True,
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
            parse_procar (bool):
                Whether to parse the PROCAR file.
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
            symmetrized_structure = None
            symmops = None
            site_symmetry = None

        outcar = parse_file(Outcar, outcar)
        total_energy = outcar.final_energy
        magnetization = 0.0 if outcar.total_mag is None else outcar.total_mag
        electrostatic_potential = outcar.electrostatic_potential

        if parse_procar and defect_entry is not None:
            if procar is None:
                procar = str(max(p.glob("**/*PROCAR*"), key=os.path.getctime))
            procar = Procar(procar)

            neighboring_atoms = defect_entry.perturbed_sites

            vbm_index = \
                {Spin.up: round((outcar.nelect + magnetization) / 2) - 1,
                 Spin.down: round((outcar.nelect - magnetization) / 2) - 1}

            participation_ratio = defaultdict(dict)
            orbital_character = defaultdict(dict)

            for s in Spin.up, Spin.down:
                if s not in eigenvalues.keys():
                    continue
                # The k-point indices at the band edges in defect calculations.
                for i, band_edge in enumerate(["vbm", "cbm"]):
                    band_index = vbm_index[s] + i

                    max_eigenvalue = np.amax(eigenvalues[s][:, band_index, 0])
                    kpoint_index = \
                        int(np.where(eigenvalues[s][:, band_index, 0]
                                     == max_eigenvalue)[0][0])

                    participation_ratio[s][band_edge] = \
                        calc_participation_ratio(
                            procar=procar,
                            spin=s,
                            band_index=band_index,
                            kpoint_index=kpoint_index,
                            atom_indices=neighboring_atoms)

                    orbital_character[s][band_edge] = \
                        calc_orbital_character(
                            procar=procar,
                            structure=final_structure,
                            spin=s,
                            band_index=band_index,
                            kpoint_index=kpoint_index)

            participation_ratio = dict(participation_ratio)
            orbital_character = dict(orbital_character)

        else:
            participation_ratio = None
            orbital_character = None

        if referenced_dft_results:
            relative_total_energy = \
                total_energy - referenced_dft_results.total_energy

            if defect_entry is None:
                raise FileNotFoundError("DefectEntry is necessary.")

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

        else:
            relative_total_energy = None
            relative_potential = None
            displacements = None

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
                   relative_total_energy=relative_total_energy,
                   relative_potential=relative_potential,
                   displacements=displacements,
                   symmetrized_structure=symmetrized_structure,
                   symmops=symmops,
                   participation_ratio=participation_ratio,
                   orbital_character=orbital_character)

    @classmethod
    def from_dict(cls, d):
        """ Construct a class object from a dictionary. """

        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        final_structure = d["final_structure"]
        if isinstance(final_structure, dict):
            final_structure = Structure.from_dict(final_structure)

        symmops = []
        for symmop in d["symmops"]:
            if isinstance(symmop, dict):
                symmops.append(SymmOp.from_dict(symmop))
            else:
                symmops.append(symmop)

        if d["participation_ratio"] is not None:
            participation_ratio = {}
            for spin, v in d["participation_ratio"].items():
                participation_ratio[Spin(int(spin))] = v
        else:
            participation_ratio = None

        if d["orbital_character"] is not None:
            orbital_character = {}
            for spin, v in d["participation_ratio"].items():
                orbital_character[Spin(int(spin))] = v
        else:
            orbital_character = None

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
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self.eigenvalues.items()}
        if self.participation_ratio is not None:
            participation_ratio = {str(spin): v
                                   for spin, v in self.participation_ratio.items()}
        else:
            participation_ratio = None

        if self.orbital_character is not None:
            orbital_character = {str(spin): v
                                 for spin, v in self.orbital_character.items()}
        else:
            orbital_character = None

        d = {"@module":                 self.__class__.__module__,
             "@class":                  self.__class__.__name__,
             "final_structure":         self.final_structure,
             "site_symmetry":           self.site_symmetry,
             "total_energy":            self.total_energy,
             "total_magnetization":     self.total_magnetization,
             "eigenvalues":             eigenvalues,
             "kpoint_coords":          self.kpoint_coords,
             "kpoint_weights":         self.kpoint_weights,
             "electrostatic_potential": self.electrostatic_potential,
             "eigenvalue_properties":   self.eigenvalue_properties,
             "volume":                  self.volume,
             "fermi_level":             self.fermi_level,
             "is_converged":            self.is_converged,
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

