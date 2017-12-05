#!/usr/bin/env python

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, get_el_sp
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun, Outcar

__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

class DefectProperty:
    """
    Object for properties of defect.
    This information is used for correction and analyzation.
    Args:
        structure (Structure, of pymatgen package): Structure object, indicating relaxed structure including a point defect (when using from_file, read file named CONTCAR, not POSCAR-finish (for Oba-lab) ).
        energy (float): energy of VASP result ()
        atomic_site_pot: potential of atomic site
    """

    def __init__(self, structure, energy, atomic_site_pot):
        self.__structure = structure
        self.__energy = energy
        self.__atomic_site_pot = atomic_site_pot
    
    @property
    def structure(self):
        return self.__structure

    @property
    def energy(self):
        return self.__energy

    @property
    def atomic_site_pot(self):
        return self.__atomic_site_pot

    @staticmethod
    def from_directory(dirname):
        """
        Reads a DefectProperty object from directory.

        Args:
            dirname (str): Directory name of VASP result.

        Returns:
            DefectProperty object
        """

        poscar_path = dirname + "/CONTCAR" 
        outcar_path = dirname + "/OUTCAR"
        vasprun_path = dirname + "/vasprun.xml"
        poscar = Poscar.from_file(poscar_path)
        outcar = Outcar(outcar_path)
        vasprun = Vasprun(vasprun_path)
        structure = poscar.structure
        energy = vasprun.final_energy
        atomic_site_pot = outcar.electrostatic_potential
        return DefectProperty(structure, energy, atomic_site_pot)

    def as_dict(self):
        d = {"energy" : self.energy,
             "structure" : self.structure,
             "atomic_site_pot" : self.atomic_site_pot}
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(__energy=d["energy"],
                   __structure=d["structure"],
                   __atomic_site_pot=d["atomic_site_pot"])

    def write_file(self):
        pass

    @classmethod
    def from_file(cls, filename):
        pass
