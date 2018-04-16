# -*- coding: utf-8 -*-

from collections import namedtuple
import numpy as np
import os
import tempfile
import unittest

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

from pydefect.analysis.defect_energies import DefectEnergies, Defect
from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults
from pydefect.core.defect_entry import DefectEntry

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "Feb. 25, 2018"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "core")


class Correction:
    def __init__(self):
        self._total_correction_energy = 0.0

    @property
    def total_correction_energy(self):
        return self._total_correction_energy


class DefectEnergiesTest(unittest.TestCase):

    def setUp(self):
        """ """
        unitcell_file = os.path.join(test_dir, "MgO/defects/unitcell.json")
        unitcell = UnitcellDftResults.json_load(unitcell_file)
        perfect_file = os.path.join(test_dir,
                                    "MgO/defects/perfect/dft_results.json")
        perfect = SupercellDftResults.json_load(perfect_file)

        # defect_dirs = ["Va_O1_1", "Va_O1_2", "Va_O1_0"]
        defect_dirs = ["Mg_O1_0", "Mg_O1_1", "Mg_O1_2", "Mg_O1_3", "Mg_O1_4",
                       "Mg_i1_0", "Mg_i1_1", "Mg_i1_2", "Va_O1_1", "Va_O1_2",
                       "Va_O1_0"]
        defects = []
        for dd in defect_dirs:
            d = os.path.join(test_dir, "MgO/defects", dd)
            defect_entry = DefectEntry.\
                json_load(os.path.join(d, "defect_entry.json"))
            dft_results = SupercellDftResults.\
                json_load(os.path.join(d, "dft_results.json"))
            correction = Correction()

            defect = Defect(defect_entry=defect_entry,
                            dft_results=dft_results,
                            correction=correction)
            defects.append(defect)

        # temporary insert values
        chem_pot = {"A": {"Mg": -2.1, "O": 0}, "B": {"Mg": 0, "O": -2.1}}
        chem_pot_label = "A"

        self.defect_energies = DefectEnergies(unitcell=unitcell,
                                              perfect=perfect,
                                              defects=defects,
                                              chem_pot=chem_pot,
                                              chem_pot_label=chem_pot_label,
                                              filtering_words=["Va_O"],
                                              system_name="MgO")

    def test_energies(self):
        print(self.defect_energies._defect_energies)
        print(self.defect_energies.vbm)
        print(self.defect_energies.cbm)
        print(self.defect_energies.band_gap)

    def test_calc_transition_levels(self):
        self.defect_energies.calc_transition_levels()
        print(self.defect_energies._transition_levels)
        self.defect_energies.plot_energy()
#        self.defect_energies.plot_energy(x_range=[-0.5, 5], y_range=[-5, 20])
#        self.defect_energies.plot_energy(file_name="test.eps")

        # # electrostatic_potential is a property because it is used for
        # # test_relative_potential method.
        # self.electrostatic_potential = \
        #     [-34.69, -35.5244, -35.5244, -35.5244, -35.5244, -35.5244, -35.5244,
        #      -34.59, -70.0739, -70.0739, -70.0739, -70.0739, -70.0739, -70.0739,
        #      -70.4981]

        # self._MgO_Va_O1_2 = SupercellDftResults(
        #     final_structure, total_energy, eigenvalues,
        #     self.electrostatic_potential)

        # self._MgO_Va_O1_2_from_vasp_files = \
        #     SupercellDftResults.from_vasp_files(
        #         os.path.join(test_dir, "MgO/defects/Va_O1_2"))

        # self._MgO_perfect_from_vasp_files = \
        #     SupercellDftResults.from_vasp_files(
        #         os.path.join(test_dir, "MgO/defects/perfect"))

        # self.d = self._MgO_Va_O1_2.as_dict()
        # self.d_from_vasp_files = self._MgO_Va_O1_2_from_vasp_files.as_dict()

        # name = "Va_O1"
        # initial_structure = Structure.from_file(
        #     os.path.join(test_dir, "MgO/defects/Va_O1_2", "POSCAR"))
        # removed_atoms = {8: [0.25, 0.25, 0.25]}
        # inserted_atoms = {}
        # changes_of_num_elements = {"O": -1}
        # charge = 2

        # self._defect_entry_MgO_Va_O1_2 = \
        #     DefectEntry(name, initial_structure, removed_atoms, inserted_atoms,
        #                 changes_of_num_elements, charge)

    # def test_from_vasp_files(self):
    #     # CAUTION: When constructing Structure object from Structure.from_file
    #     #          velocities are not stored.
    #     #          Therefore, equality check of Structure objects returns False.
    #     #          If the structure is converted via poscar file format, it may
    #     #          be solved.
    #     # contcar = Poscar.from_file(os.path.join(DIRNAME_UNITCELL, "CONTCAR"))
    #     # final_structure = contcar.structure

        # # self.assertTrue(self.d["initial_structure"] ==
        # #                 self.d_from_vasp_files["initial_structure"])
        # # self.assertTrue(self.d["final_structure"] ==
        # #                 self.d_from_vasp_files["final_structure"])
        # self.assertEqual(self.d["total_energy"],
        #                  self.d_from_vasp_files["total_energy"])
        # self.assertTrue((self.d["eigenvalues"]["1"] ==
        #                  self.d_from_vasp_files["eigenvalues"]["1"]))
        # self.assertTrue(self.d["electrostatic_potential"] ==
        #                 self.d_from_vasp_files["electrostatic_potential"])

    # def test_dict(self):
    #     MgO_Va_O1_2_fd = SupercellDftResults.from_dict(self.d_from_vasp_files)
    #     np.testing.assert_equal(MgO_Va_O1_2_fd.eigenvalues[Spin.up],
    #                             self._MgO_Va_O1_2.eigenvalues[Spin.up])

    # def test_json(self):
    #     tmp_file = tempfile.NamedTemporaryFile()
    #     self._MgO_Va_O1_2.to_json_file(tmp_file.name)
    #     MgO_Va_O1_2_from_json = SupercellDftResults.json_load(tmp_file.name)
    #     np.testing.assert_equal(MgO_Va_O1_2_from_json.eigenvalues[Spin.up],
    #                             self._MgO_Va_O1_2.eigenvalues[Spin.up])

    # def test_relative_total_energy(self):
    #     actual = self._MgO_Va_O1_2_from_vasp_files.\
    #         relative_total_energy(self._MgO_perfect_from_vasp_files)

        # expected = -93.76904720 - -95.46878101

        # self.assertEqual(actual, expected)

    # def test_relative_potential(self):
    #     actual = self._MgO_Va_O1_2_from_vasp_files.\
    #         relative_potential(self._MgO_perfect_from_vasp_files,
    #                            self._defect_entry_MgO_Va_O1_2)

        # perfect_potential = [-35.2923, -35.2923, -35.2923, -35.2923, -35.2923,
        #                      -35.2923, -35.2923, -35.2923, -69.8160, -69.8160,
        #                      -69.8160, -69.8160, -69.8160, -69.8160, -69.8160]

        # expected = [x - y for x, y in
        #             zip(self.electrostatic_potential, perfect_potential)]

        # self.assertTrue(actual == expected)


# class UnitcellDftResultsTest(unittest.TestCase):

    # def setUp(self):
    #     """ """
    #     self.unitcell = UnitcellDftResults(band_edge=None,
    #                                        band_edge2=None,
    #                                        static_dielectric_tensor=None,
    #                                        ionic_dielectric_tensor=None,
    #                                        total_dos=None)

    # def test_band_edge(self):
    #     print(self.unitcell.band_edge)
    #     self.unitcell.band_edge = [0.0, 1.0]
    #     print(self.unitcell.band_edge)
    #     self.unitcell.set_band_edge_from_vasp(
    #         directory_path=os.path.join(test_dir,
    #                                     "MgO/unitcell/structure_optimization"))
    #     print(self.unitcell.band_edge)

    # def test_band_edge2(self):
    #     print(self.unitcell.band_edge2)
    #     self.unitcell.band_edge2 = [0.0, 2.0]
    #     print(self.unitcell.band_edge2)
    #     self.unitcell.set_band_edge2_from_vasp(
    #         directory_path=os.path.join(test_dir,
    #                                     "MgO/unitcell/structure_optimization"))
    #     print(self.unitcell.band_edge2)

    # def test_dielectric_constant(self):
    #     print(self.unitcell.static_dielectric_tensor)
    #     print(self.unitcell.ionic_dielectric_tensor)
    #     print(self.unitcell.total_dielectric_tensor)
    #     self.unitcell.static_dielectric_tensor = \
    #         np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    #     self.unitcell.ionic_dielectric_tensor = \
    #         np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    #     print(self.unitcell.static_dielectric_tensor)
    #     print(self.unitcell.ionic_dielectric_tensor)
    #     print(self.unitcell.total_dielectric_tensor)
    #     self.unitcell.set_static_dielectric_tensor_from_vasp(
    #         directory_path=os.path.join(test_dir,
    #                                     "MgO/unitcell/dielectric_constants"))
    #     self.unitcell.set_ionic_dielectric_tensor_from_vasp(
    #         directory_path=os.path.join(test_dir,
    #                                     "MgO/unitcell/dielectric_constants"))
    #     print(self.unitcell.static_dielectric_tensor)
    #     print(self.unitcell.ionic_dielectric_tensor)
    #     print(self.unitcell.total_dielectric_tensor)

    # def test_total_dos(self):
    #     print(self.unitcell.total_dos)
    #     self.unitcell.total_dos = \
    #         np.array([[0, 0.5, 1.0], [0.1, 0.2, 0.4]])
    #     print(self.unitcell.total_dos)
    #     self.unitcell.set_total_dos_from_vasp(
    #         directory_path=os.path.join(test_dir,
    #                                     "MgO/unitcell/structure_optimization"))
    #     print(self.unitcell.total_dos)

    # def test_dict_json(self):
    #     self.unitcell.band_edge = [0.0, 1.0]
    #     self.unitcell.band_edge2 = [0.0, 2.0]
    #     self.unitcell.static_dielectric_tensor = \
    #         np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    #     self.unitcell.ionic_dielectric_tensor = \
    #         np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    #     self.unitcell.total_dos = \
    #         np.array([[0, 0.5, 1.0], [0.1, 0.2, 0.4]])
    #     d = self.unitcell.as_dict()
    #     unitcell_from_dict = UnitcellDftResults.from_dict(d)
    #     print(vars(unitcell_from_dict))

        # tmp_file = tempfile.NamedTemporaryFile()
        # self.unitcell.to_json_file(tmp_file.name)
        # unitcell_from_json = UnitcellDftResults.json_load(tmp_file.name)
        # print(vars(unitcell_from_json))


if __name__ == "__main__":
    unittest.main()

