import unittest

from pymatgen.core.structure import Structure

from pydefect.input_maker.input_maker import *
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.core.site import IrreducibleSite

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

DEFAULT_POTCAR_DIR = "/home/common/default_POTCAR"
FILENAME_PerturbAroundAPointTest_POSCAR = "examples/POSCAR-MgO64atoms"
FILENAME_POSCAR_VA1 = "examples/POSCAR-MgO64atoms-Va_Mg1"
FILENAME_POSCAR_IN1 = "examples/POSCAR-MgO64atoms-O_i1"
FILENAME_POSCAR_IN3 = "examples/POSCAR-MgO64atoms-Al_i1"
FILENAME_POSCAR_AS1 = "examples/POSCAR-MgO64atoms-O_Mg"
FILENAME_POSCAR_SS1 = "examples/POSCAR-MgO64atoms-N_O1"


class NormalizedRandom3dVectorTest(unittest.TestCase):

    def setUp(self): 
        self.v = normalized_random_3d_vector()

    def test_norm(self):
        print("normalized_random_3d_vector: ", self.v)
        print("norm: ", np.linalg.norm(self.v))
        self.assertAlmostEqual(np.linalg.norm(self.v), 1.0) 


class RandomVectorTest(unittest.TestCase):

    def setUp(self): 
        self.distance = 3.0
        normalized_v = normalized_random_3d_vector()
        self.v = random_vector(normalized_v, self.distance)

    def test_norm(self):
        print("random_3d_vector: ", self.v)
        print("distance: ", self.distance)
        print("norm: ", np.linalg.norm(self.v))
        self.assertLessEqual(np.linalg.norm(self.v), self.distance)


class PerturbAroundAPointTest(unittest.TestCase):

    def setUp(self):
        self.structure = \
            Structure.from_file(FILENAME_PerturbAroundAPointTest_POSCAR)
        self.center = [0.0, 0.0, 0.0]
        self.cutoff = 3.0
        self.distance = 0.2

    def test(self):
        # TODO: test the displacement distances
        perturbed_defect_structure, perturbed_sites = \
            perturb_around_a_point(self.structure, self.center, self.cutoff,
                                   self.distance)
        true_perturbed_sites = [0, 40, 44, 48, 50, 56, 57]
        self.assertEqual(perturbed_sites, true_perturbed_sites)


class DefectMakerTest(unittest.TestCase):

    def setUp(self):
        self.structure = \
            Structure.from_file(FILENAME_PerturbAroundAPointTest_POSCAR)
        Mg1 = IrreducibleSite(irreducible_name="Mg1", element="Mg",
                              first_index=1, last_index=32,
                              repr_coords=[0, 0, 0])
        O1 = IrreducibleSite(irreducible_name="O1", element="O",
                             first_index=33, last_index=64,
                             repr_coords=[0.25, 0.25, 0.25])
        self.irreducible_sites = [Mg1, O1]
        self.interstitial_coords = [[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]

        _test_structure = Structure.from_file(FILENAME_POSCAR_VA1)
        _removed_atom_index = 0
        _inserted_atom_index = None
        _defect_coords = [0, 0, 0]
        _in_name = "Va"
        _out_name = "Mg1"
        _charge = -2

        self.test_d = DefectEntry(_test_structure, _removed_atom_index,
                                  _inserted_atom_index, _defect_coords,
                                  _in_name, _out_name, _charge)

    def test_vacancies(self):
        vac1_d = \
            DefectMaker("Va_Mg1_-2", self.structure, self.irreducible_sites,
                        self.interstitial_coords)
        self.assertEqual(vac1_d.defect.as_dict(), self.test_d.as_dict())


if __name__ == "__main__":
    unittest.main()


