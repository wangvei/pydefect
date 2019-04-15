# -*- coding: utf-8 -*-

import numpy as np
import os
import unittest

from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

from pydefect.util.structure import perturb_neighboring_atoms, \
    atomic_distances, create_saturated_interstitial_structure, count_equivalent_clusters

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "input_maker")


class PerturbNeighborsTest(unittest.TestCase):

    def test(self):
        structure = \
            Structure.from_file(os.path.join(test_dir, "POSCAR-MgO64atoms"))
        center = [0.0, 0.0, 0.0]
        cutoff = 3.0
        distance = 0.2

        # TODO: test the displacement distances
        perturbed_defect_structure, perturbed_sites = \
            perturb_neighboring_atoms(structure, center, cutoff, distance)
        true_perturbed_sites = [0, 40, 44, 48, 50, 56, 57]
        self.assertEqual(perturbed_sites, true_perturbed_sites)


        self.structure = Structure.from_file("POSCAR-Sn2Nb2O7")

    # def test_get_coordination_environment(self):
    #     print(get_coordination_environment(self.structure, 0))
    #     # print(get_coordination_environment(self.structure, 4))
    #     # print(get_coordination_environment(self.structure, 8))
    #     # print(get_coordination_environment(self.structure, 16))

    # def test_get_neighbors(self):
    #     a = get_coordination_distances(self.structure, 0)
    #     for k, v in a.items():
    #         print(k + ": " + " ".join([str(round(i, 2)) for i in v]))


class AtomicDistancesTest(PymatgenTest):

    def setUp(self):
        structure = Structure.from_file("POSCAR-atom_distances")
        self.lattice = structure.lattice

    def test(self):
        points = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.15, 0.05]]
        expected = np.array([1.1435112, 1.4864129,  2.123447])
        actual = atomic_distances(self.lattice, points)
        self.assertArrayAlmostEqual(actual, expected)


class CreateSaturatedInterstitialStructureTest(PymatgenTest):
    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgO64atoms")

    def test(self):
        inserted_atom_coords = [[0.125, 0.125, 0.125]]
        saturated_defect_struct, inserted_atom_indices = create_saturated_interstitial_structure(self.structure, inserted_atom_coords)
        print(saturated_defect_struct)
        print(inserted_atom_indices)


class CountEquivalentClustersTest(PymatgenTest):

    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgO64atoms")
        self.inserted_atom_coords = [[0.125, 0.125, 0.125]]
        self.removed_atom_coords = [0, 32]

    def test(self):
        print(count_equivalent_clusters(self.structure, self.inserted_atom_coords, self.removed_atom_coords))


if __name__ == "__main__":
    unittest.main()
