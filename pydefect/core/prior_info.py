# -*- coding: utf-8 -*-
import json

from monty.serialization import loadfn
from monty.json import MSONable

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class PriorInfo(MSONable):
    """ Prior information used for input of first-principles calculations

    Args:
        energy_per_atom (float):
        band_gap (float)

    """
    def __init__(self, energy_per_atom=None, band_gap=None,
                 total_magnetization=None, data_source=None,
                 is_molecule=None, mag_threshold=0.01, band_gap_threshold=0.1):
        self.energy_per_atom = energy_per_atom
        self.band_gap = band_gap
        self.total_magnetization = total_magnetization
        self.data_source = data_source
        self.is_molecule = is_molecule
        self.mag_threshold = mag_threshold
        self.band_gap_threshold = band_gap_threshold

    def dump_json(self, filename="prior_info.json"):
        with open(filename, "w") as fw:
            json.dump(self.as_dict(), fw, indent=2)

    @classmethod
    def load_json(cls, filename="prior_info.json"):
        return loadfn(filename)

    @property
    def is_magnetic(self):
        return self.total_magnetization > self.mag_threshold

    @property
    def has_band_gap(self):
        return self.band_gap > self.band_gap_threshold

    @property
    def is_metal(self):
        return self.band_gap < self.band_gap_threshold
