# -*- coding: utf-8 -*-

import json
import numpy as np
import os

from pathlib import Path

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from xml.etree.cElementTree import ParseError

from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.defect import DefectEntry
from pydefect.core.error_classes import NoConvergenceError
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import get_displacements


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
                 kpoints: Kpoints,
                 electrostatic_potential: list,
                 eigenvalue_properties,
                 volume: float,
                 fermi_level: float,
                 is_converged: bool,
                 shallow: bool = None,
                 relative_total_energy: float = None,
                 relative_potential: list = None,
                 displacements: list = None,
                 symmetrized_structure: Structure = None):
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
            kpoints (Kpoints): 
                parsed IBZKPT file.
            electrostatic_potential (list):
                Atomic site electrostatic potential.
            volume (float):
                Volume of the supercell.
            fermi_level (float):
               Fermi level in the absolute scale.
            is_converged (bool):
                Whether the calculation is converged or not.
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
        self.kpoints = kpoints
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
                "shallow: \n" + str(self.shallow)]

        if self.kpoints:
            outs.append("IBZKPT is parsed \n")
        return "\n".join(outs)

    @classmethod
    def from_vasp_files(cls,
                        directory_path: str,
                        vasprun: str = None,
                        ibzkpt: str = None,
                        contcar: str = None,
                        outcar: str = None,
                        check_shallow: bool = True,
                        referenced_dft_results=None,
                        defect_entry: DefectEntry = None,
                        symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                        angle_tolerance: float = ANGLE_TOL):
        """ Constructs class object from vasp output files.

        Args:
            directory_path (str):
                path to the directory storing calc results.
            vasprun (str):
                Name of the vasprun.xml file.
            ibzkpt (str):
                Name of the IBZKPT file.
            contcar (str):
                Name of the converged CONTCAR file.
            outcar (str):
                Name of the OUTCAR file.
            check_shallow (bool):
                check if the defect is shallow or not.
            referenced_dft_results (SupercellCalcResults):
                SupercellDftResults object for referenced supercell dft results.
            defect_entry (DefectEntry):

            symprec (float):
                Distance precision used for symmetry analysis.
            angle_tolerance (float):
                Angle precision used for symmetry analysis.
        """
        p = Path(directory_path)

        # get the names of the latest files in the directory_path
        if vasprun is None:
            vasprun = str(max(p.glob("**/*vasprun*"), key=os.path.getctime))
        if ibzkpt is None:
            ibzkpt = str(max(p.glob("**/*IBZKPT*"), key=os.path.getctime))
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

        is_converged = True
        vasprun = parse_file(Vasprun, vasprun)
        eigenvalues = vasprun.eigenvalues
        # (band gap, cbm, vbm, is_band_gap_direct).
        eigenvalue_properties = vasprun.eigenvalue_band_properties
        fermi_level = vasprun.efermi

        # Check if the electronic and ionic steps are converged.
        if vasprun.converged_electronic is False:
            raise NoConvergenceError("Electronic step is not converged.")

        if vasprun.converged_ionic is False:
            logger.warning("Ionic step is not converged.")
            is_converged = False

        kpoints = parse_file(Kpoints.from_file, ibzkpt)

        contcar = parse_file(Poscar.from_file, contcar)
        final_structure = contcar.structure
        volume = contcar.structure.volume

        sga = SpacegroupAnalyzer(final_structure, symprec, angle_tolerance)
        site_symmetry = sga.get_point_group_symbol()

        outcar = parse_file(Outcar, outcar)
        total_energy = outcar.final_energy
        magnetization = outcar.total_mag

        if magnetization is None:
            magnetization = 0.0

        electrostatic_potential = outcar.electrostatic_potential

        shallow = None
        if referenced_dft_results:
            if check_shallow:
                shallow = []
                if abs(magnetization - round(magnetization))\
                        > MAGNETIZATION_THRESHOLD:
                    shallow.append("not_integer_mag")

                defect_supercell_vbm = eigenvalue_properties[2]
                defect_supercell_cbm = eigenvalue_properties[1]

                perfect_supercell_vbm = \
                    referenced_dft_results.eigenvalue_properties[2]
                perfect_supercell_cbm = \
                    referenced_dft_results.eigenvalue_properties[1]
                if defect_supercell_cbm > perfect_supercell_cbm - 0.1:
                    shallow.append("fermi_is_higher_than_cbm")
                if defect_supercell_vbm < perfect_supercell_vbm + 0.1:
                    shallow.append("fermi_is_lower_than_vbm")

                # In case of LDA or GGA, the defect orbital might be partially
                # occupied.
                # nelect = outcar.nelect
                # if round(nelect) % 2 != round(magnetization) % 2:
                #     shallow.append("not_localized_orbital")

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
                   kpoints=kpoints,
                   electrostatic_potential=electrostatic_potential,
                   eigenvalue_properties=eigenvalue_properties,
                   volume=volume,
                   fermi_level=fermi_level,
                   is_converged=is_converged,
                   shallow=shallow,
                   relative_total_energy=relative_total_energy,
                   relative_potential=relative_potential,
                   displacements=displacements)

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

        return cls(final_structure=final_structure,
                   site_symmetry=d["site_symmetry"],
                   total_energy=d["total_energy"],
                   total_magnetization=d["total_magnetization"],
                   eigenvalues=eigenvalues,
                   kpoints=d["kpoints"],
                   electrostatic_potential=d["electrostatic_potential"],
                   eigenvalue_properties=d["eigenvalue_properties"],
                   volume=d["volume"],
                   fermi_level=d["fermi_level"],
                   is_converged=d["is_converged"],
                   shallow=d["shallow"],
                   relative_total_energy=d["relative_total_energy"],
                   relative_potential=d["relative_potential"],
                   displacements=d["displacements"],
                   symmetrized_structure=d["symmetrized_structure"])

    @classmethod
    def load_json(cls, filename):
        """ Construct a class object from a json file.

        defaultdict is imperative to keep the backward compatibility.
        For instance, when one adds new attributes, they do not exist in old
        json files. Then, the corresponding values are set to None.
        """
#        dd = defaultdict(lambda: None, loadfn(filename))
#        return cls.from_dict(dd)
        return loadfn(filename)

    def as_dict(self):
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self.eigenvalues.items()}
        d = {"@module":                 self.__class__.__module__,
             "@class":                  self.__class__.__name__,
             "final_structure":         self.final_structure,
             "site_symmetry":           self.site_symmetry,
             "total_energy":            self.total_energy,
             "total_magnetization":     self.total_magnetization,
             "eigenvalues":             eigenvalues,
             "kpoints":                 self.kpoints,
             "electrostatic_potential": self.electrostatic_potential,
             "eigenvalue_properties":   self.eigenvalue_properties,
             "volume":                  self.volume,
             "fermi_level":             self.fermi_level,
             "is_converged":            self.is_converged,
             "shallow":                 self.shallow,
             "relative_total_energy":   self.relative_total_energy,
             "relative_potential":      self.relative_potential,
             "displacements":           self.displacements,
             "symmetrized_structure":   self.symmetrized_structure}

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
        #        print(defect_entry)
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

    @property
    def is_shallow(self):
        if self.shallow:
            return True
        else:
            return self.shallow

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
