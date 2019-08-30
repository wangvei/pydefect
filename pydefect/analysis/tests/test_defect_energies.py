# -*- coding: utf-8 -*-

import tempfile

from chempotdiag.chem_pot_diag import ChemPotDiag
from pydefect.analysis.defect import Defect
from pydefect.analysis.defect_energies import DefectEnergies
from pydefect.core.defect_entry import DefectEntry
from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.util.testing import PydefectTest
from pydefect.util.tools import sanitize_keys_in_dict

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ConvertStrKeysTest(PydefectTest):
    def setUp(self):
        d = {"Va_O1_0": {"0": {"null": 1.0}, "1": {"inward": 2.0}}}
        print(type(d["Va_O1_0"]))
        print(sanitize_keys_in_dict(d))


class DefectEnergiesTest(PydefectTest):

    def setUp(self):
        """ """
        filename = (self.TEST_FILES_DIR / "defects" / "MgO" / "unitcell.json")
        unitcell = UnitcellCalcResults.load_json(filename)

        filename = ["defects", "MgO", "perfect", "dft_results.json"]
        perfect = self.get_object_by_name(
            SupercellCalcResults.load_json, filename)

        defect_dirs = ["Va_O1_0", "Va_O1_1", "Va_O1_2", "Va_Mg1_-2"]
        defects = []
        for dd in defect_dirs:
            filename = ["defects", "MgO", dd, "defect.json"]
            defects.append(self.get_object_by_name(Defect.load_json, filename))

        # temporary insert values
        filename = ["defects", "MgO", "vertices_MgO.yaml"]
        chem_pot = self.get_object_by_name(
            ChemPotDiag.load_vertices_yaml, filename)
        chem_pot_label = "A"

        self.defect_energies = \
            DefectEnergies.from_objects(unitcell=unitcell,
                                        perfect=perfect,
                                        defects=defects,
                                        chem_pot=chem_pot,
                                        chem_pot_label=chem_pot_label,
                                        system="MgO")

    def test_msonable(self):
        self.assertMSONable(self.defect_energies)

    def test_print(self):
        print(self.defect_energies)

    def test_energies(self):
        d = self.defect_energies.as_dict()
        de = DefectEnergies.from_dict(d)
        dd = de.as_dict()
        self.assertEqual(d, dd)
        print(d)
        print(dd)

    def test_json(self):
        """ round trip test of to_json and from_json """
        tmp_file = tempfile.NamedTemporaryFile()
        self.defect_energies.to_json_file(tmp_file.name)
        defect_entry_from_json = DefectEntry.load_json(tmp_file.name)
        print(defect_entry_from_json.as_dict())
        print(self.defect_energies.as_dict())
        self.assertEqual(defect_entry_from_json.as_dict(),
                         self.defect_energies.as_dict())

    def test_U(self):
        actual, _ = self.defect_energies.u(name="Va_O1", charges=[0, 1, 2])
        expected = 1.072359471189877
        self.assertAlmostEqual(expected, actual, 7)

    def test_calc_transition_levels(self):
        dp = self.defect_energies
#        dp = DefectEnergyPlotter(self.energies, self.dc2)
#        plt = dp.plot_energy(filtering_words=["Va_O[0-9]+"],
        plt = dp.plot_energy(x_range=[-0.15, 4.5],
                             fermi_levels=[[1000, 4.5], [298, 3.3]],
                             show_transition_levels=True,
                             show_all_energies=True)
        plt.show()
#        plt.savefig(fname="energy.eps")

