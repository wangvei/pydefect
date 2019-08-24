# -*- coding: utf-8 -*-

import os

import numpy as np
from pydefect.util.structure_tools import (
    perturb_neighboring_atoms, get_minimum_distance, get_displacements, defect_center_from_coords,
    distance_list,
    atomic_distances, create_saturated_interstitial_structure,
    count_equivalent_clusters, get_point_group_op_number,
    get_symmetry_multiplicity, get_coordination_distances)
from pymatgen.core.structure import Structure
from vise.util.structure_handler import get_symmetry_dataset
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class PerturbNeighborsTest(PydefectTest):

    def setUp(self) -> None:
        self.structure = self.get_structure_by_name("MgO64atoms")
        # Mg at [0, 0, 0] is assumed to be inserted
        center = [0.5, 0.5, 0.5]
        cutoff = self.structure.lattice.a / 4
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
        expected = [43, 47, 53, 55, 62, 63]
        self.assertEqual(expected, self.perturbed_sites)

    def test_perturbed_sites_larger(self):
        expected = [12, 13, 14, 15, 18, 19, 22, 23, 25, 27, 29, 31, 43, 47, 53,
                    55, 62, 63]
        self.assertEqual(expected, self.perturbed_sites_2)

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
        actual = get_minimum_distance(self.structure)
        expected = 2.3103506255400017
        self.assertEqual(expected, actual)


class GetDisplacementsTest(PydefectTest):
    def setUp(self):
        contcar = self.get_structure_by_name("MgO64atoms-Va_O_0")
        poscar = self.get_structure_by_name("MgO64atoms")
        center = [0.25, 0.25, 0.25]
        anchor_atom_index = 14

        d = get_displacements(contcar, poscar, center, anchor_atom_index)
        self.disp_vectors = d["displacement_vectors"]
        self.disp_norms = d["displacement_norms"]

    def test(self):
        # disp_vector is in cartesian coordinates.
        vector_expected = [0.02123447, -0.10617233,  0.04246893]
        norm_expected = 0.11630595530799753
        self.assertArrayAlmostEqual(self.disp_vectors[1], vector_expected)
        self.assertAlmostEqual(self.disp_norms[1], norm_expected)


class DefectCenterFromCoordsTest(PydefectTest):
    def setUp(self):
        self._defect_coords = [[0.1, 0.2, 0.3], [-0.1, 0, 0.1]]
        self._poscar = Structure.from_file(
            os.path.join(test_dir_core, "Va_O1_2", "POSCAR"))

    def test(self):
        actual = defect_center_from_coords(self._defect_coords, self._poscar)
        print(actual)
        expected = [0.0, 0.1, 0.2]
        self.assertArrayAlmostEqual(actual, expected)


class AtomicDistancesTest(PydefectTest):

    def setUp(self):
        structure = Structure.from_file("POSCAR-atom_distances")
        self.lattice = structure.lattice

    def test(self):
        points = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.15, 0.05]]
        expected = np.array([2.123447, 1.1435112, 1.4864129])
        actual = atomic_distances(self.lattice, points)
        self.assertArrayAlmostEqual(actual, expected)


class CreateSaturatedInterstitialStructureTest(PydefectTest):
    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgO64atoms")

    def test(self):
        inserted_atom_coords = [[0.125, 0.125, 0.125]]
        saturated_defect_struct, atom_indices, are_inserted = \
            create_saturated_interstitial_structure(self.structure,
                                                    inserted_atom_coords)
        print(saturated_defect_struct)
        print(atom_indices)
        print(are_inserted)


class CountEquivalentClustersTest(PydefectTest):

    def setUp(self):
        self.structure = Structure.from_file("POSCAR-MgO64atoms")
        self.inserted_atom_coords = [[0.125, 0.125, 0.125]]
        self.removed_atom_coords = [0, 32]

    def test(self):
        print(count_equivalent_clusters(self.structure, self.inserted_atom_coords, self.removed_atom_coords))


class GetPointGroupOpNumberTest(PydefectTest):
    def setUp(self):
        structure = Structure.from_file("POSCAR-MgO64atoms")
        self.sym_dataset = get_symmetry_dataset(structure)
        self.lattice = structure.lattice.matrix

    def test(self):
        coords = [0.125, 0.125, 0.125]
        print(get_point_group_op_number(self.sym_dataset, coords, self.lattice))
        print(get_symmetry_multiplicity(self.sym_dataset, coords, self.lattice))


class GetCoordinationDistancesTest(PydefectTest):
    def setUp(self) -> None:
        self.structure = self.get_structure_by_name("KZn4P3")

    def test_get_coordination_environment(self):
        actual = get_coordination_distances(self.structure,
                                            atom_index=1, cutoff=2.5)
        expected = {'P': [2.31, 2.5, 2.5, 2.5]}
        self.assertEqual(expected, actual)