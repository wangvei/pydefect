# -*- coding: utf-8 -*-

import json
import os

import numpy as np
from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn
from pydefect.util.logger import get_logger
from pydefect.util.tools import make_symmetric_matrix
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Outcar, Vasprun, Poscar
from vise.analyzer.band_gap import band_gap_properties

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class UnitcellCalcResults(MSONable):
    """Container class with DFT results for a unitcell. """

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
        self._band_edge = band_edge[:] if band_edge else None
        self._static_dielectric_tensor = static_dielectric_tensor
        self._ionic_dielectric_tensor = ionic_dielectric_tensor
        self._total_dos = total_dos
        self._volume = volume
        self.is_direct = is_direct

    def __repr__(self):

        def to(s):
            return str(s) if s else "None"

        band_edge = self._band_edge if self._band_edge else [None, None]
        is_total_dos = "Exists" if self._total_dos else "None"
        outs = \
            [f"vbm (eV): {to(band_edge[0])}",
             f"cbm (eV): {to(band_edge[1])}",
             f"static dielectric tensor: {to(self._static_dielectric_tensor)}",
             f"ionic dielectric tensor: {to(self._ionic_dielectric_tensor)}",
             f"total dielectric tensor: {to(self.total_dielectric_tensor)}",
             f"volume (A^3): {to(self._volume)}",
             f"Total DOS: {is_total_dos}"]

        return "\n".join(outs)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)

    @staticmethod
    def check_attribute(name, attr):
        if attr is None:
            logger.warning(f"{name}: is None.")
            return
        return attr

    @property
    def band_edge(self):
        return self.check_attribute("band edge", self._band_edge)

    @property
    def static_dielectric_tensor(self):
        return self.check_attribute("static dielectric tensor",
                                    self._static_dielectric_tensor)

    @property
    def ionic_dielectric_tensor(self):
        return self.check_attribute("ionic dielectric tensor",
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
        return self.check_attribute("total dos", self._total_dos)

    @property
    def volume(self):
        return self.check_attribute("volume", self._volume)

    @property
    def is_set_all(self):
        return all([self._band_edge,
                    self._static_dielectric_tensor,
                    self._ionic_dielectric_tensor,
                    self._total_dos,
                    self._volume])

    @band_edge.setter
    def band_edge(self, band_edge: float):
        self._band_edge = band_edge

    @static_dielectric_tensor.setter
    def static_dielectric_tensor(self, d: list):
        self._static_dielectric_tensor = make_symmetric_matrix(d)

    @ionic_dielectric_tensor.setter
    def ionic_dielectric_tensor(self, d: list):
        self._ionic_dielectric_tensor = make_symmetric_matrix(d)

    @total_dos.setter
    def total_dos(self, total_dos):
        self._total_dos = total_dos

    @volume.setter
    def volume(self, volume: float):
        self._volume = volume

    # setter from vasp results
    def set_band_edge_from_vasp(self,
                                directory_path: str,
                                vasprun_name: str = "vasprun.xml") -> None:
        vasprun = Vasprun(os.path.join(directory_path, vasprun_name))

        # 2019/7/13 NEVER USE Vasprun.eigenvalue_band_properties
        # THERE IS A BUG TO ESTIMATE VBM AND CBM of lower band gap materials.
        _, vbm_info, cbm_info = band_gap_properties(vasprun)
        self.is_direct = vbm_info["kpoints"] == cbm_info["kpoints"]
        self._band_edge = [vbm_info["energy"], cbm_info["energy"]]

    def set_static_dielectric_tensor_from_vasp(self,
                                               directory_path: str,
                                               outcar_name: str = "OUTCAR"
                                               ) -> None:
        outcar = Outcar(os.path.join(directory_path, outcar_name))
        outcar.read_lepsilon()
        self._static_dielectric_tensor = np.array(outcar.dielectric_tensor)

    def set_ionic_dielectric_tensor_from_vasp(self,
                                              directory_path: str,
                                              outcar_name: str = "OUTCAR"
                                              ) -> None:
        outcar = Outcar(os.path.join(directory_path, outcar_name))
        outcar.read_lepsilon_ionic()
        self._ionic_dielectric_tensor = np.array(outcar.dielectric_ionic_tensor)

    def set_total_dos_from_vasp(self,
                                directory_path: str,
                                vasprun_name: str = "vasprun.xml") -> None:
        vasprun = Vasprun(os.path.join(directory_path, vasprun_name))
        dos, energies = vasprun.tdos.densities, vasprun.tdos.energies
        # only non-magnetic
        self._total_dos = [list(dos[Spin.up]), list(energies)]

    def set_volume_from_vasp(self,
                             directory_path: str,
                             contcar_name: str = "CONTCAR") -> None:
        contcar = Poscar.from_file(os.path.join(directory_path, contcar_name))
        self._volume = contcar.structure.volume

    def to_json_file(self, filename: str = "unitcell.json") -> None:
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)
