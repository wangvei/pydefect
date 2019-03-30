# -*- coding: utf-8 -*-

import json
import numpy as np
import os
import logging
import warnings

from collections import defaultdict
from monty.json import MontyEncoder, MSONable
from monty.serialization import loadfn

from pymatgen.io.vasp.outputs import Outcar, Vasprun, Poscar
from pymatgen.electronic_structure.core import Spin


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = logging.getLogger(__name__)


class UnitcellDftResults(MSONable):
    """
    DFT results for a unitcell.
    Args:
        band_edge (1x2 list): VBM and CBM.
        static_dielectric_tensor (3x3 numpy array):
        ionic_dielectric_tensor (3x3 numpy array):
        total_dos (2xN numpy array): [[energy1, dos1], [energy2, dos2],...]
        volume (float):
    """

    def __init__(self, band_edge=None, static_dielectric_tensor=None,
                 ionic_dielectric_tensor=None, total_dos=None, volume=None):
        """ """
        self._band_edge = band_edge
        self._static_dielectric_tensor = static_dielectric_tensor
        self._ionic_dielectric_tensor = ionic_dielectric_tensor
        self._total_dos = total_dos
        self._volume = volume

    def __eq__(self, other):
        if other is None or type(self) != type(other):
            raise TypeError
        return self.__dict__ == other.__dict__

    def __str__(self):

        def xstr(s):
            return 'None' if s is None else str(s)

        if self._band_edge is None:
            band_edge = [None, None]
        else:
            band_edge = self._band_edge

        outs = ["vbm (eV): " + xstr(band_edge[0]),
                "cbm (eV): " + xstr(band_edge[1]),
                "static dielectric tensor:",
                xstr(self._static_dielectric_tensor),
                "ionic dielectric tensor:",
                xstr(self._ionic_dielectric_tensor),
                "total dielectric tensor:",
                xstr(self.total_dielectric_tensor),
                "volume (A^3): " + xstr(self._volume)]

        if self._total_dos is None:
            outs.append("Total density of states is not set yet")
        else:
            outs.append("Total density of states is already set.")

        return "\n".join(outs)

    @classmethod
    def load_json(cls, filename):
        """
        Constructs a class object from a json file.
        """
        return loadfn(filename)

    @staticmethod
    def warning_message(name, attr):
        if attr is None:
            logger.warning("Attribute {}: is not set yet.".format(name))
            return None
        return attr

    @property
    def band_edge(self):
        return self.warning_message("band edge", self._band_edge)

    @property
    def static_dielectric_tensor(self):
        return self.warning_message("static dielectric tensor",
                                    self._static_dielectric_tensor)

    @property
    def ionic_dielectric_tensor(self):
        return self.warning_message("ionic dielectric tensor",
                                    self._ionic_dielectric_tensor)

    @property
    def total_dielectric_tensor(self):
        if self.static_dielectric_tensor and self.ionic_dielectric_tensor:
            return self._static_dielectric_tensor + self._ionic_dielectric_tensor
        else:
            return None

    @property
    def total_dos(self):
        return self.warning_message("total dos", self._total_dos)

    @property
    def volume(self):
        return self.warning_message("volume", self._volume)

    def is_set_all(self):
        return not any([self._band_edge,
                        self._static_dielectric_tensor,
                        self._ionic_dielectric_tensor,
                        self._total_dos,
                        self._volume])

    @band_edge.setter
    def band_edge(self, band_edge):
        self._band_edge = band_edge

    @static_dielectric_tensor.setter
    def static_dielectric_tensor(self, d):
        """
        d (list):
            len(d) == 1: Suppose cubic system
            len(d) == 3: Suppose tetragonal or orthorhombic system
            len(d) == 6: Suppose the other system
        """
        if len(d) == 1:
            self._static_dielectric_tensor = \
                np.array([[d[0], 0, 0],
                          [0, d[0], 0],
                          [0, 0, d[0]]])
        elif len(d) == 3:
            self._static_dielectric_tensor = \
                np.array([[d[0], 0, 0],
                          [0, d[1], 0],
                          [0, 0, d[2]]])
        elif len(d) == 6:
            self._static_dielectric_tensor = \
                np.array([[d[0], d[1], d[2]],
                          [d[1], d[3], d[4]],
                          [d[2], d[4], d[5]]])

    @ionic_dielectric_tensor.setter
    def ionic_dielectric_tensor(self, d):
        if len(d) == 1:
            self._ionic_dielectric_tensor = \
                np.array([[d[0], 0, 0],
                          [0, d[0], 0],
                          [0, 0, d[0]]])
        elif len(d) == 3:
            self._ionic_dielectric_tensor = \
                np.array([[d[0], 0, 0],
                          [0, d[1], 0],
                          [0, 0, d[2]]])
        else:
            self._ionic_dielectric_tensor = \
                np.array([[d[0], d[1], d[2]],
                          [d[1], d[3], d[4]],
                          [d[2], d[4], d[5]]])

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
        _, cbm, vbm, _ = vasprun.eigenvalue_band_properties
        self._band_edge = [vbm, cbm]

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
        dos, ene = vasprun.tdos.densities, vasprun.tdos.energies
        # only non-magnetic
        self._total_dos = np.vstack([dos[Spin.up], ene])

    def set_volume_from_vasp(self, directory_path, contcar_name="CONTCAR"):
        contcar = Poscar.from_file(os.path.join(directory_path, contcar_name))
        self._volume = contcar.structure.volume

    def to_json_file(self, filename):
        """
        Returns a json file, named unitcell.json.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)
