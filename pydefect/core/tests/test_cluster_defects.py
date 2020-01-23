# -*- coding: utf-8 -*-

import filecmp
import tempfile
from collections import OrderedDict

from pydefect.core.cluster_defects import ClusterDefect, ClusterDefects
from pydefect.util.testing import PydefectTest

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ClusterDefectTest(PydefectTest):
    def setUp(self):
        # 0.375000 0.375000 0.625000 Cu+
        # 0.625000 0.625000 0.625000 Cu+
        self.cu2o = ClusterDefect(removed_atom_indices=[9, 31],
                                  inserted_atoms=[
                                      {"element": "Cu",
                                       "coords": [0.5, 0.5, 0.625]}],
                                  point_group="mm2",
                                  multiplicity=96,
                                  extreme_charge_state=-1,
                                  annotation="test")

    def test_dict(self):
        d = self.cu2o.as_dict()
        rounded_d = ClusterDefect.from_dict(d).as_dict()
        self.assertEqual(d, rounded_d)

    def test_msonalbe(self):
        self.assertMSONable(self.cu2o)


class ClusterDefectsTest(PydefectTest):
    def setUp(self):
        self.structure = self.get_structure_by_name("Cu2O48atoms")

        self.split = ClusterDefect(removed_atom_indices=[9, 31],
                                   inserted_atoms=[{"element": "Cu",
                                               "coords": [0.5, 0.5, 0.625]}],
                                   point_group="mm2",
                                   multiplicity=96,
                                   extreme_charge_state=-1)

        self.cu2o = \
            ClusterDefects(self.structure, OrderedDict({"split": self.split}))

    def test_dict(self):
        d = self.cu2o.as_dict()
        point_group = d["cluster_defects"]["split"]["point_group"]
        self.assertEqual(point_group,"mm2")
        rounded_d = ClusterDefects.from_dict(d).as_dict()
        self.assertEqual(d, rounded_d)

    def test_yaml(self):
        tmp_file = tempfile.NamedTemporaryFile().name
        self.cu2o.site_set_to_yaml_file(tmp_file)
        with open(tmp_file) as f:
            print(f.read())
        self.assertTrue(filecmp.cmp("expected_cluster_defects.yaml", tmp_file))

    def test_msonable(self):
        self.assertMSONable(self.cu2o)

    def test_from_files(self):
        actual = ClusterDefects.from_files(
            structure=self.structure,
            yaml_filename="expected_cluster_defects.yaml").as_dict()
        expected = self.cu2o.as_dict()
        self.assertEqual(expected, actual)

    def test_add(self):
        self.cu2o_added = ClusterDefects(structure=self.structure)
        self.cu2o_added.add_defect(
            removed_atom_indices=[9, 31],
            inserted_atoms=[{"element": "Cu", "coords": [0.5, 0.5, 0.625]}],
            name="split",
            extreme_charge_state=-1)

        self.assertEqual(self.cu2o.as_dict(), self.cu2o_added.as_dict())



