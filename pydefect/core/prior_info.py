# -*- coding: utf-8 -*-
import json

from monty.serialization import loadfn

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


class PriorInfo:

    def __init__(self, energy_per_atom=None, band_gap=None,
                 total_magnetization=None, data_source=None,
                 mag_threshold=0.001, band_gap_threshold=0.1):
        self._energy_per_atom = energy_per_atom
        self._band_gap = band_gap
        self._total_magnetization = total_magnetization
        self._data_source = data_source
        self._mag_threshold = mag_threshold
        self._band_gap_threshold = band_gap_threshold

    def as_dict(self):
        d = {"energy_per_atom": self._energy_per_atom,
             "band_gap": self._band_gap,
             "total_magnetization": self._total_magnetization,
             "data_source": self._data_source}
        return d

    def dump_json(self, filename="prior_info.json"):
        with open(filename, "w") as fw:
            json.dump(self.as_dict(), fw, indent=2)

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a DefectEntry class object from a dictionary.
        """
        # The keys need to be converted to integers.

        return cls(d["energy_per_atom"], d["band_gap"],
                   d["total_magnetization"])

    @classmethod
    def load_json(cls, filename="prior_info.json"):
        """
        Constructs a DefectEntry class object from a json file.
        """
        return cls.from_dict(loadfn(filename))

    @property
    def energy_per_atom(self):
        return self._energy_per_atom

    @property
    def band_gap(self):
        return self._band_gap

    @property
    def total_magnetization(self):
        return self._total_magnetization

    @property
    def is_magnetic(self):
        if self._total_magnetization > self._mag_threshold:
            return True

    @property
    def is_band_gap(self):
        if self._band_gap > self._band_gap_threshold:
            return True