# -*- coding: utf-8 -*-
from abc import ABC, abstractmethod
import json

from monty.json import MSONable, MontyEncoder
from monty.serialization import loadfn


class Correction(ABC, MSONable):
    @property
    @abstractmethod
    def correction_energy(self):
        pass

    def to_json_file(self, filename: str) -> None:
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)


class ManualCorrection(Correction):

    method = "manual_correction"

    def __init__(self,
                 manual_correction_energy: float = 0.0):
        """ Set manual correction energy.

        Used when users want to their own correction values.
        Args:
            manual_correction_energy (float):
        """
        self._manual_correction_energy = manual_correction_energy

    @property
    def correction_energy(self) -> float:
        return self._manual_correction_energy

