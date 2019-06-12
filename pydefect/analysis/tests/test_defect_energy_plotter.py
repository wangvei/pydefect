# -*- coding: utf-8 -*-

import os
import unittest

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.analysis.defect_carrier_concentration import DefectConcentration
from pydefect.analysis.defect_energy_plotter import DefectEnergyPlotter
from pydefect.analysis.chempotdiag.chem_pot_diag import ChemPotDiag
from pydefect.corrections.corrections import ExtendedFnvCorrection
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class DefectEnergiesTest(unittest.TestCase):

    def setUp(self):
        """ """
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        unitcell = UnitcellCalcResults.load_json(unitcell_file)
        perfect_file = os.path.join(test_dir,
                                    "MgO/defects/perfect/dft_results.json")
        perfect = SupercellCalcResults.load_json(perfect_file)

        # defect_dirs = ["Mg_O1_0", "Mg_O1_1", "Mg_O1_2", "Mg_O1_3", "Mg_O1_4",
        #                "Mg_i1_0", "Mg_i1_1", "Mg_i1_2", "Va_O1_1", "Va_O1_2",
        #                "Va_O1_0"]
        defect_dirs = ["Mg_O1_0", "Mg_O1_1", "Mg_O1_2", "Mg_O1_3", "Mg_O1_4",
                       "Mg_i1_0", "Mg_i1_1", "Mg_i1_2", "Va_O1_0"]
        defects = []
        for dd in defect_dirs:
            d = os.path.join(test_dir, "MgO/defects", dd)
            defect_entry = \
                DefectEntry.load_json(os.path.join(d, "defect_entry.json"))
            dft_results = \
                SupercellCalcResults.load_json(
                    os.path.join(d, "dft_results.json"))
            correction = \
                ExtendedFnvCorrection.load_json(os.path.join(d, "correction.json"))

            defect = Defect(defect_entry=defect_entry,
                            dft_results=dft_results,
                            correction=correction)

            defects.append(defect)

        # temporary insert values
        chem_pot = ChemPotDiag.load_vertices_yaml(
            os.path.join(test_dir, "MgO/vertices_MgO.yaml"))

        chem_pot_label = "A"

        self.defect_energies = \
            DefectEnergies.from_files(unitcell=unitcell,
                                      perfect=perfect,
                                      defects=defects,
                                      chem_pot=chem_pot,
                                      chem_pot_label=chem_pot_label,
                                      system="MgO")

        # temperature = 10000
        # temperature2 = 1000
        # num_sites_filename = os.path.join(test_dir,
        #                                   "MgO/defects/num_sites.yaml")

        # dc1 = DefectConcentration.from_defect_energies(
        #     defect_energies=self.defect_energies,
        #     temperature=temperature,
        #     unitcell=unitcell,
        #     num_sites_filename=num_sites_filename)

        # self.dc2 = DefectConcentration.from_defect_energies(
        #     defect_energies=self.defect_energies,
        #     temperature=temperature2,
        #     unitcell=unitcell,
        #     num_sites_filename=num_sites_filename,
        #     previous_concentration=dc1,
        #     verbose=False)

    def test(self): pass

    def test_calc_transition_levels(self):
        dp = DefectEnergyPlotter(self.defect_energies)
#        dp = DefectEnergyPlotter(self.defect_energies, self.dc2)
#        plt = dp.plot_energy(filtering_words=["Va_O[0-9]+"],
        plt = dp.plot_energy(x_range=[-1, 5],
                             show_fermi_level=True,
                             show_transition_levels=True,
                             show_all_energies=True)
        plt.show()
#        plt.savefig(fname="energy.eps")


if __name__ == "__main__":
    unittest.main()

