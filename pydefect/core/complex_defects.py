# -*- coding: utf-8 -*-

from collections import OrderedDict
from typing import Optional

import yaml
from monty.json import MSONable
from pydefect.core.interstitial_site import represent_odict, construct_odict
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

# See document in interstitial_site.py
yaml.add_representer(OrderedDict, represent_odict)
yaml.add_constructor('tag:yaml.org,2002:map', construct_odict)


class ComplexDefect(MSONable):
    """

    Args:
        representative_coords (list):
            Representative coordinates, namely the position of first_index
        wyckoff (str):
            A wyckoff letter.
        site_symmetry (str):
            Site symmetry.
        coordination_distances (dict):
            Coordination environment. An example is
            {"Mg": [1.92, 1.95, 2.01], "Al": [1.82, 1.95]}
        method (str):
            The method name determining the interstitial site.
    """

    def __init__(self,
                 name: str,
                 oxidation_state,
                 annotation,
                 removed_atom_indices: Optional[list] = None,
                 inserted_atoms: Optional[dict] = None):

        self.name = name
        self.oxidation_state = oxidation_state
        self.annotation = annotation
        self.removed_atom_indices = \
            removed_atom_indices if removed_atom_indices else []
        self.inserted_atoms = inserted_atoms if inserted_atoms else {}

    def __repr__(self):
        outs = [f"name: {self.name}",
                f"oxidation_state: {self.oxidation_state}",
                f"annotation: {self.annotation}"]
        return outs
