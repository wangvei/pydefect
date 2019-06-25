# -*- coding: utf-8 -*-
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.corrections.corrections import Correction

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class Defect:
    def __init__(self,
                 defect_entry: DefectEntry,
                 dft_results: SupercellCalcResults,
                 correction: Correction):

        self.defect_entry = defect_entry
        self.dft_results = dft_results
        self.correction = correction

    @property
    def is_shallow(self):
        for i in self.dft_results.band_edges.values():
            if i.is_shallow:
                return True
        return False
