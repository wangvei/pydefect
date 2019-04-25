# -*- coding: utf-8 -*-

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import DummySpecie

from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.error_classes import StructureError
from pydefect.util.structure_tools import first_appearance_index, \
    get_symmetry_multiplicity, create_saturated_interstitial_structure
from pydefect.util.logger import get_logger
from pydefect.core.interstitial_site import InterstitialSite, \
    InterstitialSiteSet
from obadb.util.structure_handler import get_coordination_distances

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class InterstitialSiteMaker:
    def __init__(self,
                 interstitial_set: InterstitialSiteSet):

        self.interstitial_set = interstitial_set

    def add_interstitial(self,
                         structure: Structure,
                         specie: str,
                         coords: list,
                         site_name: str,
                         r: float = 0.3,
                         force_add: bool = False,
                         symprec: float = SYMMETRY_TOLERANCE,
                         angle_tolerance: float = ANGLE_TOL,
                         method: str = "manual"):
        """

        1. Check whether other sites exist near the inserted sites.


        """
        # DummySpecie locate at the saturated interstitial sites.
        saturated_structure = \
            create_saturated_interstitial_structure(structure,
                                                    coords,
                                                    dist_tol=symprec)[0]
        cart_coords = saturated_structure.lattice.get_cartesian_coords(coords)
        neighbors = saturated_structure.get_sites_in_sphere(pt=cart_coords, r=r)

        if len(neighbors) > 1:
            for n in neighbors:
                if isinstance(n[0], DummySpecie):
                    specie = "another interstitial site"
                else:
                    specie = n[0]
                distance = n[1]
                if distance < 1e-8:
                    pass
                message = "Inserted position is too close to {}. The distance is" \
                          " {:5.3f}".format(specie, distance)
                if force_add:
                    logger.warning(message)
                else:
                    raise StructureError(message)

        structure = structure.copy()
        inserted_index = first_appearance_index(structure, specie)
        structure.insert(inserted_index, specie, coords)

        sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
        sym_dataset = sga.get_symmetry_dataset()
        wyckoff = sym_dataset["wyckoffs"][inserted_index]
        site_symmetry = sym_dataset["pointgroup"]
        symmetry_multiplicity = \
            get_symmetry_multiplicity(sym_dataset, coords,
                                      structure.lattice.matrix, symprec)
        coordination_distances = \
            get_coordination_distances(structure, inserted_index)

        interstitial_site = \
            InterstitialSite(site_name=site_name,
                             representative_coords=coords,
                             wyckoff=wyckoff,
                             site_symmetry=site_symmetry,
                             symmetry_multiplicity=symmetry_multiplicity,
                             coordination_distances=coordination_distances,
                             method=method)

        self.interstitial_set.append(interstitial_site)
