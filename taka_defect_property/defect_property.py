#!/usr/bin/env python

import sys
import json
import copy

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

class DefectProperty():
    """
    Object for properties of defect.
    This information is used for correction and analyzation.
    Args:
        structure (Structure, of pymatgen package): Structure object, indicating relaxed structure including a point defect (when using from_file, read file named CONTCAR, not POSCAR-finish (for Oba-lab) ).
        energy (float): energy of VASP result ()
        atomic_site_pot: potential of atomic site
    """

    def __init__(self, structure, energy, atomic_site_pot, charge):
        self._structure = structure
        self._energy = energy
        self._atomic_site_pot = atomic_site_pot
        self._charge = charge
    
    @property
    def structure(self):
        return self._structure

    @property
    def energy(self):
        return self._energy

    @property
    def atomic_site_pot(self):
        return copy.copy(self._atomic_site_pot)

    @property
    def charge(self):
        return self._charge

    @staticmethod
    def from_directory(dirname, charge):
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
        charge = -2 #TODO temporary
        return DefectProperty(structure, energy, atomic_site_pot, charge)

    def as_dict(self):
        d = {"energy" : self.energy,
             "structure" : self.structure,
             "atomic_site_pot" : self.atomic_site_pot,
             "charge" : charge}
        return d

    @classmethod
    def from_dict(cls, d):
        st_item = d["structure"]
        if isinstance(st_item, Structure):
            structure = st_item
        elif isinstance(st_item, str):
            structure = Structure.from_str(st_item, fmt="json")
        elif isinstance(st_item, dict):
            structure = Structure.from_dict(st_item)
        else:
            raise TypeError("Failed to convert an element of dictionary named 'structure', or no item named 'structure' input. ")
        return cls(energy=d["energy"],
                   structure=structure,
                   atomic_site_pot=d["atomic_site_pot"],
                   charge=d["charge"])

    @classmethod
    def from_str(cls, s): 
        print(s)
        d = json.loads(s)
        #print(d["structure"])
        #d["structure"] = Structure.from_dict(d["structure"])
        return(cls.from_dict(d))

    def to(self, filename=None):
        d = self.as_dict()
        d["structure"] = d["structure"].to(fmt="json")
        s = json.dumps(d, indent=4)
        if filename:
            with open(filename, 'w') as f:
                f.write(s)
        else:
            return s

    @classmethod
    def from_file(cls, filename):
        with open(filename) as f:
            d = json.load(f)
        return(cls.from_dict(d))
