import tempfile
import unittest

import ruamel.yaml

from pydefect.analysis.chempotdiag.compound import Compound, CompoundsList
from pydefect.analysis.chempotdiag.vertex import Vertex, VerticesList
from pydefect.analysis.chempotdiag.chem_pot_diag import ChemPotDiag

EXAMPLE_DIR = "../../examples/"

FILENAME_1D = EXAMPLE_DIR + "energy_1d.txt"
FILENAME_2D = EXAMPLE_DIR + "energy_2d.txt"
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
DFT_DIRECTORIES = [EXAMPLE_DIR + "/dft_data/O2/",
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


class TestVertex(unittest.TestCase):

    def test_vertex_dict(self):
        v = Vertex("Label", ["Mg", "O"], [0.00000, -7.77345])
        d = v.as_dict()
        v2 = Vertex.from_dict(d)
        self.assertEqual(v, v2)
        # TODO: test for dummy vertex


if __name__ == "__main__":
    unittest.main()
