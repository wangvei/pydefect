# -*- coding: utf-8 -*-

from collections import defaultdict
import json
import numpy as np
import os

from pathlib import Path

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from xml.etree.cElementTree import ParseError

from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.electronic_structure.core import Spin
from pymatgen.core import Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pydefect.core.defect_entry import DefectEntry
from pydefect.util.logger import get_logger
from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class SupercellCalcResults(MSONable):
    """ Holds DFT results for supercell systems both w/ and w/o a defect.
    """

    def __init__(self,
                 final_structure: Structure,
                 total_energy: float,
                 total_magnetization: float,
                 eigenvalues: np.array,
                 kpoints: Kpoints,
                 electrostatic_potential: list,
                 eigenvalue_properties,
                 volume: float,
                 fermi_level: float,
                 is_converged: bool,
                 does_symmetrize: bool = True):
        """
        Args:
            final_structure (Structure):
                pmg Structure class object. Usually relaxed structures
            total_energy (float):
                Final total energy in eV.
            total_magnetization (float):
                Total total_magnetization in \mu_B
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
        """
        self.final_structure = final_structure
        self.total_energy = total_energy
        self.total_magnetization = total_magnetization
        self.eigenvalues = eigenvalues
        self.kpoints = kpoints
        self.electrostatic_potential = electrostatic_potential
        self.eigenvalue_properties = eigenvalue_properties
        self.volume = volume
        self.fermi_level = fermi_level
        self.is_converged = is_converged

        self.symmetrized_structure = None
        if does_symmetrize:
            sga = SpacegroupAnalyzer(self.final_structure,
                                     DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL)

            self.symmetrized_structure = sga.get_symmetrized_structure()

    def __str__(self):
        outs = ["total energy (eV): " + str(self.total_energy),
                "total total_magnetization (mu_B): " +
                str(self.total_magnetization),
                "electrostatic potential: " + str(self.electrostatic_potential),
                "eigenvalues_properties: " + str(self.eigenvalue_properties),
                "final structure: \n" + str(self.final_structure),
                "volume: \n" + str(self.volume),
                "Fermi level (eV): \n" + str(self.fermi_level)]

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
                        does_symmetrize: bool = True):
        """
        Construct class object from vasp output files.
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
            does_symmetrize (bool):
                Whether to symmetrize the structure.
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
        eigenvalue_properties = vasprun.eigenvalue_band_properties
        fermi_level = vasprun.efermi
        # Check if the electronic and ionic steps are converged.
        if vasprun.converged_electronic is False:
            logger.warning("Electronic step is not converged.")
            is_converged = False
        if vasprun.converged_ionic is False:
            logger.warning("Electronic step is not converged.")
            is_converged = False

        kpoints = parse_file(Kpoints.from_file, ibzkpt)

        contcar = parse_file(Poscar.from_file, contcar)
        final_structure = contcar.structure
        volume = contcar.structure.volume

        outcar = parse_file(Outcar, outcar)
        total_energy = outcar.final_energy
        magnetization = outcar.total_mag
        if not magnetization:
            magnetization = 0.0
        electrostatic_potential = outcar.electrostatic_potential

        return cls(final_structure, total_energy, magnetization, eigenvalues,
                   kpoints, electrostatic_potential, eigenvalue_properties,
                   volume, fermi_level, is_converged, does_symmetrize)

    @classmethod
    def from_dict(cls, d):
        """
        Construct a class object from a dictionary.
        """
        eigenvalues = {}
        # Programmatic access to enumeration members in Enum class.
        for spin, v in d["eigenvalues"].items():
            eigenvalues[Spin(int(spin))] = np.array(v)

        final_structure = d["final_structure"]
        if isinstance(final_structure, dict):
            final_structure = Structure.from_dict(final_structure)

        return cls(final_structure=final_structure,
                   total_energy=d["total_energy"],
                   total_magnetization=d["total_magnetization"],
                   eigenvalues=eigenvalues,
                   kpoints=d["kpoints"],
                   electrostatic_potential=d["electrostatic_potential"],
                   eigenvalue_properties=d["eigenvalue_properties"],
                   volume=d["volume"],
                   fermi_level=d["fermi_level"],
                   is_converged=d["is_converged"])

    @classmethod
    def load_json(cls, filename):
        """
        Construct a class object from a json file.
        defaultdict is imperative to keep the backward compatibility.
        For instance, when one adds new attributes, they do not exist in old
        json files. Then, the corresponding values are set to None.
        """
#        dd = defaultdict(lambda: None, loadfn(filename))
#        return cls.from_dict(dd)
        print(filename)
        return loadfn(filename)

    def as_dict(self):
        """
        Dict representation of DefectInitialSetting class object.
        Json-serializable dict representation.
        """
        # Spin object must be converted to string for to_json_file.
        eigenvalues = {str(spin): v.tolist()
                       for spin, v in self.eigenvalues.items()}
        d = {"@module":                 self.__class__.__module__,
             "@class":                  self.__class__.__name__,
             "final_structure":         self.final_structure,
             "total_energy":            self.total_energy,
             "total_magnetization":     self.total_magnetization,
             "eigenvalues":             eigenvalues,
             "kpoints":                 self.kpoints,
             "electrostatic_potential": self.electrostatic_potential,
             "eigenvalue_properties":   self.eigenvalue_properties,
             "volume":                  self.volume,
             "fermi_level":             self.fermi_level,
             "is_converged":            self.is_converged}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file, named dft_results.json.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def relative_total_energy(self, referenced_dft_results):
        """
        Return a relative total energy w.r.t. referenced supercell.

        Args:
            referenced_dft_results (SupercellCalcResults):
                SupercellDftResults object for referenced supercell dft results.
        """
        return self.total_energy - referenced_dft_results.total_energy

    def relative_potential(self, referenced_dft_results, defect_entry):
        """
        Return a list of relative site potential w.r.t. the perfect supercell.
        Note that None is inserted for interstitial sites.

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
