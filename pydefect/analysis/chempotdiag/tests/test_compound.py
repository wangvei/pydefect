import unittest

from pydefect.analysis.chempotdiag.compound import Compound, CompoundsList
from pydefect.analysis.chempotdiag.vertex import Vertex, VerticesList

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


class TestCompound(unittest.TestCase):

    def test_set_elements_of_compound(self):
        c = Compound(None, ["H", "He", "Li"], [0.0, 0.3, 0.7], -42)
        #  Can method set_elements order?
        c.set_elements(["He", "H", "Li"])
        expected = Compound(None, ["He", "H", "Li"], [0.3, 0.0, 0.7], -42)
        self.assertTrue(c.almost_equal(expected))
        #  Can method remove element?
        c.set_elements(["Li", "He"])
        expected = Compound(None, ["Li", "He"], [0.7, 0.3], -42)
        self.assertTrue(c.almost_equal(expected))
        #  Can method add element?
        c.set_elements(["He", "Li", "Be"])
        expected = Compound(None, ["He", "Li", "Be"], [0.3, 0.7, 0.0], -42)
        self.assertTrue(c.almost_equal(expected))
        #  Can method raise exception if argument is invalid?
        self.assertRaises(ValueError, lambda: c.set_elements(["U", "Pr"]))
        #  Must be contain He, whose composition is non-zero.
        self.assertRaises(ValueError, lambda: c.set_elements(["Li", "Be"]))

    def test_from_vasp_files(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        vasprun_paths = [d+VASPRUN_NAME for d in DFT_DIRECTORIES]
        #  from outcar
        cl = CompoundsList.from_vasp_calculations_files(poscar_paths,
                                                        outcar_paths)
        print(cl)

        #  from vasprun.xml
        cl = CompoundsList.from_vasp_calculations_files(poscar_paths,
                                                        vasprun_paths,
                                                        fmt="vasprun")

        # When calculation of any simple substance is not found,
        # program should be raise error.
        # TODO: This doesn't raise error
        # self.assertRaises(ValueError,
        #                   lambda: CompoundsList.from_vasp_calculations_files(
        #                       poscar_paths[1:], outcar_paths[1:]))

    def test_gas_shift(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        cl = CompoundsList.from_vasp_calculations_files(poscar_paths,
                                                        outcar_paths)
        print(cl[0].gas_energy_shift(293.15, 1))
        print(cl)
        print(cl.gas_energy_shifts(293.15, {"O2": 1}))
        print("at T=293.15, O2 pressure = 1")
        print(cl.free_energies(293.15, {"O2": 1}))
        print("at T=293.15, O2 pressure = 1e+4")
        print(cl.free_energies(293.15, {"O2": 1e+4}))
        print("at T=293.15, O2 pressure = 1e+5")
        print(cl.free_energies(293.15, {"O2": 1e+5}))
        print("at T=6000, O2 pressure = 1e+5")
        print(cl.free_energies(6000, {"O2": 1e+5}))
        print("at T=6000, O2 pressure = 1e+100")
        print(cl.free_energies(6000, {"O2": 1e+100}))


if __name__ == "__main__":
    unittest.main()
