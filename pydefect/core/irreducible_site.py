# -*- coding: utf-8 -*-
from typing import Optional
from monty.json import MSONable

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class IrreducibleSite(MSONable):
    """ Holds properties related to the symmetrically equivalent atom set.

    Note1: atomic indices need to be sorted, meaning they can be written in a
           sequence, like 17..32
    Note2: first_index atom is assumed to represent the irreducible atoms.

    Args:
        irreducible_name (str):
            Element name with the irreducible index (e.g., Mg1)
        element (str):
            Element name (e.g., Mg)
        first_index (int):
            First index of irreducible_name. Note that the index begins from 1.
        last_index (int):
            Last index of irreducible_name.
        representative_coords (list):
            Representative coordinates, namely the position of first_index
        wyckoff (str):
            A wyckoff letter.
        site_symmetry (str):
            Site symmetry.
        coordination_distances (dict):
            Coordination environment. An example is
            {"Mg": [1.92, 1.95, 2.01], "Al": [1.82, 1.95]}
        magmom (float):
            Local magnetic moment.
    """
    def __init__(self,
                 irreducible_name: str,
                 element: str,
                 first_index: int,
                 last_index: int,
                 representative_coords: list,
                 wyckoff: str,
                 site_symmetry: str,
                 coordination_distances: Optional[dict] = None,
                 magmom: Optional[float] = None):
        self.irreducible_name = irreducible_name
        self.element = element
        self.first_index = first_index
        self.last_index = last_index
        self.representative_coords = representative_coords
        self.wyckoff = wyckoff
        self.site_symmetry = site_symmetry
        self.coordination_distances = coordination_distances
        self.magmom = magmom

    @property
    def num_atoms(self):
        return self.last_index - self.first_index + 1
