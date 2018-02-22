import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 19, 2018"


class Unitcell:
    """
    DFT result of unitcell without any defect.
    """

    def __init__(self, directory_path,
                 poscar_name = "/POSCAR",
                 outcar_name = "/OUTCAR",
                 vasprun_name = "/vasprun.xml"):
        """
        Args:
            directory_path (str): path of directory. 
            poscar_name (str): Name of POSCAR file. Defaults to POSCAR.
            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
            vasprun_name (str): Name of vasprun.xml file.
                                Defaults to vasprun.xml.
        """
        path_poscar = directory_path + poscar_name
        path_outcar = directory_path + outcar_name
        path_vasprun = directory_path + vasprun_name
        poscar = Poscar.from_file(path_poscar)
        outcar = Outcar(path_outcar)
        vasprun = Vasprun(path_vasprun)
        self._structure = poscar.structure
        self._energy = outcar.final_energy
        self._dielectric_tensor \
            = np.array(outcar.dielectric_tensor) + \
            np.array(outcar.dielectric_ionic_tensor)
        self._eigen_value = vasprun.eigenvalues
        self._ewald_param = None

    @property
    def structure(self):
        return self._structure

    @property
    def energy(self):
        return self._energy

    @property
    def dielectric_tensor(self):
        return self._dielectric_tensor

    @property
    def eigen_value(self):
        return self._eigen_value

