# -*- coding: utf-8 -*-

import numpy as np
from pydefect.util.structure_tools import (
    perturb_neighboring_atoms, get_min_distance, get_displacements,
    defect_center_from_coords, distance_list,
    create_saturated_interstitial_structure, get_neighboring_atom_indices,
    num_equivalent_clusters, first_appearing_index, get_coordination_distances)
from pymatgen.core.structure import Structure
from pydefect.util.testing import PydefectTest


class PerturbNeighborsTest(PydefectTest):

    def setUp(self) -> None:
        self.structure = self.get_structure_by_name("MgO64atoms")
        # Mg at [0, 0, 0] is assumed to be inserted
        center = [0.5, 0.5, 0.5]
        cutoff = self.structure.lattice.a / 4 + 0.01
        self.distance = 0.2
        inserted_atom_indices = [7]

        self.perturbed_structure, self.perturbed_sites = \
            perturb_neighboring_atoms(self.structure, center, cutoff,
                                      self.distance, inserted_atom_indices)

        # Mg-Mg = 2.976
        cutoff_larger = cutoff * 1.42
        _, self.perturbed_sites_2 = \
            perturb_neighboring_atoms(self.structure, center, cutoff_larger,
                                      self.distance, inserted_atom_indices)

    def test_perturbed_sites_smaller(self):
        expected = {35, 39, 54, 55, 61, 63}
        self.assertEqual(sorted(expected), sorted(set(self.perturbed_sites)))

    def test_perturbed_sites_larger(self):
        expected = {12, 13, 14, 15, 18, 19, 22, 23, 25, 27, 29, 31, 35, 39, 54,
                    55, 61, 63}
        self.assertEqual(sorted(expected), sorted(set(self.perturbed_sites_2)))

    def test_displacement_distance(self):
        distances = []
        for s in [43, 47, 53, 55, 62, 63]:
            disp = (self.perturbed_structure[s].frac_coords -
                    self.structure[s].frac_coords)
            distances.append(np.linalg.norm(disp) * self.structure.lattice.a)
        print(distances)
        self.assertLessEqual(max(distances), self.distance)


class GetMinimumDistanceTest(PydefectTest):
    def setUp(self) -> None:
        self.structure = self.get_structure_by_name("KZn4P3")

    def test_minimum_distance(self):
        actual = get_min_distance(self.structure)
        expected = 2.31035
        self.assertAlmostEqual(expected, actual, 5)


class GetDisplacementsTest(PydefectTest):
    def setUp(self):
        contcar = self.get_structure_by_name("MgO64atoms-Va_O_0-relaxed")
        poscar = self.get_structure_by_name("MgO64atoms-Va_O_0-unrelaxed")
        center = [0.25, 0, 0]
        anchor_atom_index = 38

        d = get_displacements(contcar, poscar, center, anchor_atom_index)
        self.disp_vectors = d["displacement_vectors"]
        self.disp_norms = d["displacement_norms"]

    def test(self):
        # disp_vector is in cartesian coordinates.
        vector_expected = [(0.9995962318584262 - 1) * 8.419456, 0, 0]
        norm_expected = abs(0.9995962318584262 - 1) * 8.419456
        self.assertArrayAlmostEqual(vector_expected, self.disp_vectors[0], 4)
        self.assertAlmostEqual(norm_expected, self.disp_norms[0], 4)


class DefectCenterFromCoordsTest(PydefectTest):
    def test_divacancy(self):
        defect_coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
        structure = self.get_structure_by_name("MgO64atoms-Va_Mg+Va_O-unrelax")
        actual = defect_center_from_coords(defect_coords, structure)
        expected = [0.125, 0.125, 0.125]
        self.assertArrayAlmostEqual(expected, actual)


class DistanceListTest(PydefectTest):
    def test_distance_list(self):
        structure = self.get_structure_by_name("MgO")
        coords = [0.125, 0.125, 0.125]
        actual = distance_list(structure, coords)
        expected = [0.9195, 1.7607]
        self.assertArrayAlmostEqual(expected, actual, 4)


class CreateSaturatedInterstitialStructureTest(PydefectTest):
    def setUp(self):
        self.structure = self.get_structure_by_name("MgO64atoms")

    def test(self):
        inserted_atom_coords = [[0.125, 0.125, 0.125]]
        saturated_defect_structure, atom_indices, are_inserted = \
            create_saturated_interstitial_structure(self.structure,
                                                    inserted_atom_coords)
        self.assertEqual(128, len(saturated_defect_structure))
        self.assertArrayEqual(np.array(inserted_atom_coords[0]),
                              saturated_defect_structure.frac_coords[64])
        self.assertEqual([64], atom_indices)
        self.assertEqual(True, are_inserted[0])


class GetNeighboringAtomIndicesTest(PydefectTest):
    def setUp(self) -> None:
        self.structure = self.get_structure_by_name("MgO64atoms")

    def test_only_self(self):
        coords = [0.5, 0.5, 0.5]
        actual = get_neighboring_atom_indices(self.structure,
                                              coords, cutoff=2.1)
        expected = ([7], [0.0])
        self.assertEqual(expected, actual)

    def test_two_neighbors(self):
        coords_2 = [0.125, 0.0, 0.0]
        actual_2 = get_neighboring_atom_indices(self.structure,
                                                coords_2, cutoff=1.07)
        expected_2 = ([0, 32], [1.052432, 1.052432])
        self.assertEqual(expected_2, actual_2)


class CountEquivalentClustersTest(PydefectTest):

    def setUp(self):
        self.structure = self.get_structure_by_name("MgO64atoms")
        self.inserted_atom_coords = [[0.125, 0.125, 0.125]]
        self.removed_atom_indices = [0, 40]

    def test(self):
        num_sites, point_group = \
            num_equivalent_clusters(
                structure=self.structure,
                inserted_atom_coords=self.inserted_atom_coords,
                removed_atom_indices=self.removed_atom_indices)

        self.assertEqual("3m", point_group)


class FirstAppearingIndexTest(PydefectTest):
    def setUp(self):
        self.structure = self.get_structure_by_name("MgO64atoms")

    def test_not_exist_specie(self):
        actual = first_appearing_index(self.structure, "H")
        expected = 0
        self.assertEqual(expected, actual)

    def test_exist_specie(self):
        actual = first_appearing_index(self.structure, "O")
        expected = 32
        self.assertEqual(expected, actual)


class GetCoordinationDistancesTest(PydefectTest):
    def setUp(self) -> None:
        self.structure = self.get_structure_by_name("KZn4P3")

    def test_get_coordination_environment(self):
        actual = get_coordination_distances(self.structure,
                                            atom_index=1, cutoff=2.5)
        expected = {'P': [2.31, 2.5, 2.5, 2.5]}
        self.assertEqual(expected, actual)