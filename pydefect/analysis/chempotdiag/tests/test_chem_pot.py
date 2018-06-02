#!/usr/bin/env python

import tempfile
import unittest

import ruamel.yaml

from pydefect.analysis.chempotdiag.compound import Compound, CompoundsList
from pydefect.analysis.chempotdiag.vertex import Vertex, VerticesList
from pydefect.analysis.chempotdiag.chem_pot_diag import ChemPotDiag

EXAMPLE_DIR = "../../../../test_files/analysis/chempotdiag/"

FILENAME_1D = EXAMPLE_DIR + "energy_1d.txt"
FILENAME_2D = EXAMPLE_DIR + "energy_2d.txt"
# Output of commit f0d3466e9c752c74618dce1c7acfad870785b35b version.
elements_2d = ["Mg", "O"]
stable_compounds_list_2d \
    = CompoundsList([Compound("Mg", elements_2d, [1, 0], 0.0),
                     Compound("MgO", elements_2d, [0.5, 0.5], -3.08177796),
                     Compound("O", elements_2d, [0, 1], 0.0),
                     Compound("Mg2O", elements_2d, [2/3, 1/3], -2.59114942)])
unstable_compounds_list_2d \
    = CompoundsList([Compound("Mg3O", elements_2d, [3/4, 1/4], 2.0849569025),
                     Compound("MgO2", elements_2d, [1/3, 2/3], 2.1919460833)])
all_compounds_list_2d = stable_compounds_list_2d + unstable_compounds_list_2d
element_energy_2d = {"Mg": -2.11327587, "O": -2.46256238}
vertices_list_2d \
    = VerticesList([Vertex(None, elements_2d, [-6.16356,  0.00000]),
                    Vertex(None, elements_2d, [ 0.00000, -7.77345]),
                    Vertex(None, elements_2d, [-1.60989, -4.55366])])
vertex_near_mgo_2d \
    = VerticesList([Vertex(None, elements_2d, [-6.16356,  0.00000]),
                    Vertex(None, elements_2d, [-1.60989, -4.55366])])
FILENAME_3D = EXAMPLE_DIR + "energy_MP-Ca-Al-O.txt"
FILENAME_4D = EXAMPLE_DIR + "energy_4d.txt"
# For read DFT test. We don't check these files are physically proper.
DFT_DIRECTORIES = [EXAMPLE_DIR + "/dft_data/O2molecule/",
                   EXAMPLE_DIR + "/dft_data/Mg/",
                   EXAMPLE_DIR + "/dft_data/MgO/"]
POSCAR_NAME = "POSCAR-finish"
OUTCAR_NAME = "OUTCAR-finish"
VASPRUN_NAME = "vasprun-finish.xml"
# output of commit f3a7d6d73bcd9d2bb546dd707b64a4fb138a7fbe version.
elements_dft = ["Mg", "O"]
compounds_list_from_dft = \
    CompoundsList([Compound("O2", elements_dft, [0, 1], 0.0),
                   Compound("Mg", elements_dft, [1, 0], 0.0),
                   Compound("MgO", elements_dft, [0.5, 0.5], -2.83816228)])


class TestChemPot(unittest.TestCase):

    #  Simple test to generate class object from text file.
    def test_read_from_file(self):
        _ = ChemPotDiag.from_file(FILENAME_1D)
        _ = ChemPotDiag.from_file(FILENAME_2D)
        _ = ChemPotDiag.from_file(FILENAME_3D)
        _ = ChemPotDiag.from_file(FILENAME_4D)

    def test_read_from_vasp_files(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        vasprun_paths = [d+VASPRUN_NAME for d in DFT_DIRECTORIES]
        print(poscar_paths)
        #  from outcar
        cp = \
            ChemPotDiag.from_vasp_calculations_files(poscar_paths, outcar_paths)
        self.assertTrue(cp.stable_compounds.
                        almost_equal(compounds_list_from_dft))
        #  from vasprun.xml
        cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                      vasprun_paths,
                                                      fmt="vasprun")
        self.assertTrue(cp.stable_compounds.
                        almost_equal(compounds_list_from_dft))
        # When calculation of any simple substance is not found,
        # program should be raise error.
        self.assertRaises(ValueError,
                          lambda: ChemPotDiag.from_vasp_calculations_files(
                              poscar_paths[1:], outcar_paths[1:]))

    #  Test to draw 1d diagram
    def test_draw_diagram_1d(self):
        cp = ChemPotDiag.from_file(FILENAME_1D)
        cp.draw_diagram()

    def test_draw_diagram_1d_twice(self):
        cp = ChemPotDiag.from_file(FILENAME_1D)
        cp.draw_diagram(title="The first drawing of 1d_twice_test")
        cp.draw_diagram(title="The second drawing of 1d_twice_test")

    #  Test to draw 2d diagram

    def test_draw_diagram_2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        cp.draw_diagram()

    def test_draw_diagram_2d_with_vertex_name(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        rc = "MgO"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        cp.draw_diagram(title=title_str, remarked_compound=rc)

        #  Maybe label of previous plot remains if you forget to delete.
        rc = "Mg2O"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        cp.draw_diagram(title=title_str, remarked_compound=rc)

    def test_draw_diagram_2d_without_label(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        cp.draw_diagram(title="diagram_2d_without_label", with_label=False)

    def test_rearrange_and_draw_diagram2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        cp.set_elements(["O", "Mg"])
        cp.draw_diagram(title="rearrange_and_draw_diagram2d")

    def test_draw_range_diagram2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        rc = "MgO"
        draw_range = -20
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        cp.draw_diagram(title=title_str,
                        remarked_compound=rc,
                        draw_range=draw_range)

    def test_wrong_draw_range_diagram2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        draw_range = 100
        rc = "MgO"
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        self.assertRaises(ValueError,
                          lambda: cp.draw_diagram(title=title_str,
                                                  remarked_compound=rc,
                                                  draw_range=draw_range))

    def test_rearrange_and_range_diagram2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        rc = "MgO"
        draw_range = -20
        cp.set_elements(["O", "Mg"])
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        cp.draw_diagram(title=title_str,
                        remarked_compound=rc,
                        draw_range=draw_range)

    #  Test to draw 3d diagram

    def test_draw_diagram_3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        cp.draw_diagram()

    def test_draw_diagram_3d_with_vertex_name(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        rc = "Ca"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        cp.draw_diagram(title=title_str, remarked_compound=rc) # For debug
        #  Maybe label of previous plot remains if you forget to delete.
        rc = "Ca11Al14O32"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        cp.draw_diagram(title=title_str, remarked_compound=rc)
        got_vl = cp.get_neighbor_vertices(rc)
        self.assertTrue(
            got_vl[0].almost_equal(
                Vertex("A", ["Al", "Ca", "O"], [-0.25935, -1.19373, -5.64654])))
        self.assertTrue(
            got_vl[1].almost_equal(
                Vertex("B", ["Al", "Ca", "O"], [-0.70214, -1.19373, -5.45282])))
        self.assertTrue(
            got_vl[2].almost_equal(
                Vertex("C", ["Al", "Ca", "O"], [-8.88137, -6.64655, 0])))
        self.assertTrue(
            got_vl[3].almost_equal(
                Vertex("D", ["Al", "Ca", "O"], [-8.72916, -6.84027,  0])))

    def test_draw_diagram_3d_without_label(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        cp.draw_diagram(title="diagram_3d_without_label",
                        with_label=False)

    def test_rearrange_and_diagram3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        cp.set_elements(["O", "Ca", "Al"])
        cp.draw_diagram()

    def test_draw_range_diagram3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        rc = "Ca11Al14O32"
        draw_range = -20
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        cp.draw_diagram(title=title_str,
                        remarked_compound=rc,
                        draw_range=-20)

    def test_wrong_draw_range_diagram3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        draw_range = 100
        rc = "Ca11Al14O32"
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        self.assertRaises(ValueError,
                          lambda:cp.draw_diagram(title=title_str,
                                                 remarked_compound=rc,
                                                 draw_range=draw_range))

    def test_rearrange_and_range_diagram3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        rc = "Ca11Al14O32"
        draw_range = -20
        cp.set_elements(["O", "Ca", "Al"])
        title_str = "rearranged, {0} is remarked, draw_range = {1}"\
            .format(rc, draw_range)
        cp.draw_diagram(title=title_str,
                        remarked_compound=rc,
                        draw_range=draw_range)
    #  Attribute test

    def test_dim(self):
        cp1 = ChemPotDiag.from_file(FILENAME_1D)
        self.assertEqual(cp1.dim, 1)
        cp2 = ChemPotDiag.from_file(FILENAME_2D)
        self.assertEqual(cp2.dim, 2)
        cp3 = ChemPotDiag.from_file(FILENAME_3D)
        self.assertEqual(cp3.dim, 3)
        cp4 = ChemPotDiag.from_file(FILENAME_4D)
        self.assertEqual(cp4.dim, 4)

    def test_vertices(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        self.assertTrue(cp.equilibrium_points.almost_equal(vertices_list_2d))

    def test_vertices_of_remarked_compound(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        actual_vertices = cp.get_neighbor_vertices("MgO")
        self.assertTrue(actual_vertices.almost_equal(vertex_near_mgo_2d))

    def test_stable_compound_dat(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        self.assertTrue(cp.stable_compounds
                        .almost_equal(stable_compounds_list_2d))

    def test_unstable_compound_dat(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        self.assertTrue(cp.unstable_compounds
                        .almost_equal(unstable_compounds_list_2d))

    def test_all_compound_dat(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        self.assertTrue(cp.all_compounds
                        .almost_equal(all_compounds_list_2d))

    def test_element_energy(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        self.assertAlmostEqual(cp.element_energy[0], element_energy_2d["Mg"])
        self.assertAlmostEqual(cp.element_energy[1], element_energy_2d["O"])

    def test_rearrange(self):
        _ = Compound("rearrange_test", ["A", "B", "C"], [1/6, 2/6, 3/6], 42)
        cp = ChemPotDiag.from_file(FILENAME_3D)
        # initially alphabetical
        self.assertListEqual(cp.elements, ["Al", "Ca", "O"])
        cp.set_elements(["Al", "O", "Ca"])
        self.assertEquals(cp.elements, ["Al", "O", "Ca"])
        cp.set_elements(["O", "Ca", "Al"])
        self.assertEquals(cp.elements, ["O", "Ca", "Al"])

    def test_wrong_rearrange(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        self.assertRaises(ValueError,
                          lambda: cp.set_elements(["He", "Ca", "O"]))
        self.assertRaises(ValueError,
                          lambda: cp.set_elements(["Ca", "O"]))
        self.assertRaises(ValueError,
                          lambda: cp.set_elements(["Ca", "O", "Ca"]))

    # @unittest.skip("You must see many pictures.")
    def test_draw_diagram_3d_all_examples(self):
        for s1 in ("MP", "Oba"):
            for s2 in ("Ca-Al-O", "Mg-Ca-O", "Sr-Bi-N", "Sr-Fe-O", "Sr-Ti-O"):
                path = EXAMPLE_DIR + "energy_" + s1 + "-" + s2 + ".txt"
                cp = ChemPotDiag.from_file(path)
                rc = cp.stable_compounds[1].name
                cp.draw_diagram(remarked_compound=rc)

    def test_draw_diagram_2d_with_temperature_pressure(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        for t in [293.15, 500, 1000]:
            for o2_p in [1e+5, 1e+100]:
                title_str = "T = {0} (K), P_O2 = {1} (Pa)".format(t, o2_p)
                p = {"O2": o2_p}
                cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                              outcar_paths,
                                                              temperature=t,
                                                              pressure=p)
                rc = cp.stable_compounds[1].name
                # print(cp.stable_compounds)
                cp.draw_diagram(
                    title=title_str,
                    remarked_compound=rc)

        title_str = "temperature = 293.15, pressure = None"
        cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                      outcar_paths,
                                                      temperature=293.15,
                                                      pressure=None)
        rc = cp.stable_compounds[1].name
        # print(cp.stable_compounds)
        cp.draw_diagram(
            title=title_str,
            remarked_compound=rc)

    def test_draw_diagram_3d_with_temperature_pressure(self):
        # TODO
        for t in [298.15, 3000, 6000]:
            for o2_p in [1e-10, 1e+5, 1e+20]:
                title_str = "T = {0} (K), P_O2 = {1} (Pa)".format(t, o2_p)
                p = {"O2": o2_p}
                cp = ChemPotDiag.from_file(FILENAME_3D,
                                           temperature=t,
                                           pressure=p)
                rc = cp.stable_compounds[1].name
                # print(cp.stable_compounds)
                cp.draw_diagram(
                    title=title_str,
                    remarked_compound=rc)

    def test_neighbor_vertices_as_dict(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)  # Ca,Al,O
        d = cp.get_neighbor_vertices_as_dict("Ca11Al14O32", comment="comment")
        print(d)
        d = cp.get_neighbor_vertices_as_dict("Ca", comment="comment")
        print(d)

    def test_dump_yaml(self):
        dir_path = "./"

        cp = ChemPotDiag.from_file(FILENAME_2D)  # Mg, O
        comment = "This is from test_output_yaml in test_chem_pot.py"
        cp.dump_vertices_yaml(dir_path, "MgO", comment=comment)
        filename = dir_path + "/vertices_MgO.yaml"
        vertices, standard_energy = ChemPotDiag.load_vertices_yaml(filename)
        vertices.set_elements(cp.elements)
        self.assertTrue(cp.get_neighbor_vertices("MgO").
                        almost_equal(vertices))
        # TODO: element_energy may be better to be represented by dict
        for i, elem in enumerate(cp.elements):
            self.assertAlmostEqual(cp.element_energy[i],
                                   standard_energy[elem])

        cp = ChemPotDiag.from_file(FILENAME_3D)  # Ca, Al, O
        comment = "This is from test_output_yaml in test_chem_pot.py"
        cp.dump_vertices_yaml(dir_path, "Ca11Al14O32", comment=comment)
        filename = dir_path + "/vertices_Ca11Al14O32.yaml"
        vertices, standard_energy = ChemPotDiag.load_vertices_yaml(filename)
        vertices.set_elements(cp.elements)
        self.assertTrue(cp.get_neighbor_vertices("Ca11Al14O32").
                        almost_equal(vertices))
        # TODO: element_energy may be better to be represented by dict
        for i, elem in enumerate(cp.elements):
            self.assertAlmostEqual(cp.element_energy[i],
                                   standard_energy[elem])


if __name__ == "__main__":
    unittest.main()
