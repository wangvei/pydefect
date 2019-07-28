# -*- coding: utf-8 -*-

import json
import numpy as np
import os

from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from obadb.analyzer.band_gap import band_gap_properties

from pymatgen.io.vasp.outputs import Outcar, Vasprun, Poscar
from pymatgen.electronic_structure.core import Spin

from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


def check_attribute(name, attr):
    if attr is None:
        logger.warning("{}: is None.".format(name))
        return
    return attr


def make_symmetric_matrix(d):
    """
    d (list or float):
        len(d) == 1: Suppose cubic system
        len(d) == 3: Suppose tetragonal or orthorhombic system
        len(d) == 6: Suppose the other system
    """
    if isinstance(d, float):
        tensor = np.array([[d, 0, 0],
                           [0, d, 0],
                           [0, 0, d]])
    elif len(d) == 1:
        tensor = np.array([[d[0], 0,  0],
                           [0,  d[0], 0],
                           [0,  0,  d[0]]])
    elif len(d) == 3:
        tensor = np.array([[d[0], 0, 0],
                           [0, d[1], 0],
                           [0, 0, d[2]]])
    elif len(d) == 6:
        from pymatgen.util.num import make_symmetric_matrix_from_upper_tri
        """ 
        Given a symmetric matrix in upper triangular matrix form as flat array 
        indexes as:
        [A_xx, A_yy, A_zz, A_xy, A_xz, A_yz]
        This will generate the full matrix:
        [[A_xx, A_xy, A_xz], [A_xy, A_yy, A_yz], [A_xz, A_yz, A_zz]
        """
        tensor = make_symmetric_matrix_from_upper_tri(d)
    else:
        raise ValueError("{} is not valid to make symmetric matrix".format(d))

    return tensor


class UnitcellCalcResults(MSONable):
    """ DFT results for a unitcell. """

    def __init__(self,
                 band_edge: list = None,
                 static_dielectric_tensor: np.array = None,
                 ionic_dielectric_tensor: np.array = None,
                 total_dos: list = None,
                 volume: float = None,
                 is_direct: bool = None):
        """
        Args:
            band_edge (list): [VBM, CBM].
            static_dielectric_tensor (3x3 numpy array):
            ionic_dielectric_tensor (3x3 numpy array):
            total_dos (list):
                [[energy1, energy2, ...], [dos1, dos2, ...]]
            volume (float): Volume in A-3.
        """
        self._band_edge = None if band_edge is None else list(band_edge)
        self._static_dielectric_tensor = static_dielectric_tensor
        self._ionic_dielectric_tensor = ionic_dielectric_tensor
        self._total_dos = total_dos
        self._volume = volume
        self.is_direct = is_direct

    def __repr__(self):

        def xstr(s):
            return 'None' if s is None else str(s)

        band_edge = [None, None] if self._band_edge is None else self._band_edge

        outs = ["vbm (eV): " + xstr(band_edge[0]),
                "cbm (eV): " + xstr(band_edge[1]),
                "static dielectric tensor:",
                xstr(self._static_dielectric_tensor),
                "ionic dielectric tensor:",
                xstr(self._ionic_dielectric_tensor),
                "total dielectric tensor:",
                xstr(self.total_dielectric_tensor),
                "volume (A^3): " + xstr(self._volume)]

        is_total_dos = "None" if self._total_dos is None else "Exists"
        outs.append("Total DOS: {}".format(is_total_dos))

        return "\n".join(outs)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    @property
    def band_edge(self):
        return check_attribute("band edge", self._band_edge)

    @property
    def static_dielectric_tensor(self):
        return check_attribute("static dielectric tensor",
                               self._static_dielectric_tensor)

    @property
    def ionic_dielectric_tensor(self):
        return check_attribute("ionic dielectric tensor",
                               self._ionic_dielectric_tensor)

    @property
    def total_dielectric_tensor(self):
        try:
            return (np.array(self.static_dielectric_tensor) +
                    np.array(self.ionic_dielectric_tensor)).tolist()
        except TypeError:
            return

    @property
    def total_dos(self):
        return check_attribute("total dos", self._total_dos)

    @property
    def volume(self):
        return check_attribute("volume", self._volume)

    @property
    def is_set_all(self):
        for i in self._band_edge, self._static_dielectric_tensor, \
                 self._ionic_dielectric_tensor, self._total_dos, self._volume:
            if i is None:
                return False

        return True

    @band_edge.setter
    def band_edge(self, band_edge):
        self._band_edge = band_edge

    @static_dielectric_tensor.setter
    def static_dielectric_tensor(self, d):
        self._static_dielectric_tensor = make_symmetric_matrix(d)

    @ionic_dielectric_tensor.setter
    def ionic_dielectric_tensor(self, d):
        self._ionic_dielectric_tensor = make_symmetric_matrix(d)

    @total_dos.setter
    def total_dos(self, total_dos):
        self._total_dos = total_dos

    @volume.setter
    def volume(self, volume):
        self._volume = volume

    # setter from vasp results
    def set_band_edge_from_vasp(self, directory_path,
                                vasprun_name="vasprun.xml"):
        vasprun = Vasprun(os.path.join(directory_path, vasprun_name))

        # 2019/7/13 NEVER USE Vasprun.eigenvalue_band_properties
        # THERE IS A BUG TO ESTIMATE VBM AND CBM
        _, vbm_info, cbm_info = band_gap_properties(vasprun)
        self.is_direct = vbm_info["kpoints"] == cbm_info["kpoints"]
        self._band_edge = [vbm_info["energy"], cbm_info["energy"]]

    def set_static_dielectric_tensor_from_vasp(self, directory_path,
                                               outcar_name="OUTCAR"):
        outcar = Outcar(os.path.join(directory_path, outcar_name))
        self._static_dielectric_tensor = np.array(outcar.dielectric_tensor)

    def set_ionic_dielectric_tensor_from_vasp(self, directory_path,
                                              outcar_name="OUTCAR"):
        outcar = Outcar(os.path.join(directory_path, outcar_name))
        self._ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)

    def set_total_dos_from_vasp(self, directory_path,
                                vasprun_name="vasprun.xml"):
        vasprun = Vasprun(os.path.join(directory_path, vasprun_name))
        dos, energies = vasprun.tdos.densities, vasprun.tdos.energies
        # only non-magnetic
        self._total_dos = [list(dos[Spin.up]), list(energies)]

    def set_volume_from_vasp(self, directory_path, contcar_name="CONTCAR"):
        contcar = Poscar.from_file(os.path.join(directory_path, contcar_name))
        self._volume = contcar.structure.volume

    def to_json_file(self, filename="unitcell.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)
