# -*- coding: utf-8 -*-

from copy import deepcopy
from collections import OrderedDict
from monty.json import MSONable
from typing import Union
import yaml

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer

from obadb.util.structure_handler import get_coordination_distances

from pydefect.core.config \
    import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL, \
    MAX_NUM_INTERSTITIAL_SITES
from pydefect.core.error_classes import StructureError
from pydefect.util.structure_tools \
    import get_symmetry_multiplicity, create_saturated_interstitial_structure
from pydefect.util.logger import get_logger

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class InterstitialSite(MSONable):
    """Holds properties related to the interstitial site.

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
                 representative_coords: list,
                 wyckoff: str = None,
                 site_symmetry: str = None,
                 symmetry_multiplicity: int = None,
                 coordination_distances: dict = None,
                 method: str = None):
        self.representative_coords = representative_coords
        self.wyckoff = wyckoff
        self.site_symmetry = site_symmetry
        self.symmetry_multiplicity = symmetry_multiplicity
        self.coordination_distances = coordination_distances
        self.method = method

    def __repr__(self):
        outs = [f"representative_coords: {self.representative_coords}",
                f"wyckoff: {self.wyckoff}",
                f"site_symmetry: {self.site_symmetry}",
                f"symmetry_multiplicity: {self.symmetry_multiplicity}",
                f"coordination_distances: {self.coordination_distances}",
                f"method: {self.method}"]
        return "\n".join(outs)

    def as_dict(self):
        d = OrderedDict(
            {"representative_coords":  self.representative_coords,
             "wyckoff":                self.wyckoff,
             "site_symmetry":          self.site_symmetry,
             "symmetry_multiplicity":  self.symmetry_multiplicity,
             "coordination_distances": self.coordination_distances,
             "method":                 self.method})

        return d


# The followings are needed to keep order of dictionary for interstitial.yaml.
# https://qiita.com/podhmo/items/aa954ee1dc1747252436
def represent_odict(dumper, instance):
    return dumper.represent_mapping('tag:yaml.org,2002:map', instance.items())


yaml.add_representer(OrderedDict, represent_odict)


def construct_odict(loader, node):
    return OrderedDict(loader.construct_pairs(node))


yaml.add_constructor('tag:yaml.org,2002:map', construct_odict)


class InterstitialSiteSet(MSONable):
    """Holds a set of InterstitialSite objects. """

    def __init__(self,
                 structure: Structure,
                 interstitial_sites: OrderedDict = None):
        """
        Args:
            structure (Structure):
                Structure class object. Supercell used for defect
                calculations.
            interstitial_sites (OrderedDict):
                OrderedDict with keys of site names and values of
                InterstitialSite objects.
        """
        self.structure = structure
        if interstitial_sites is not None:
            self.interstitial_sites = deepcopy(interstitial_sites)
        else:
            self.interstitial_sites = OrderedDict()

    def site_set_as_dict(self):
        d = OrderedDict()
        for k, v in self.interstitial_sites.items():
            d[k] = v.as_dict()

        return d

    def as_dict(self):
        return {"interstitial_site_set": self.site_set_as_dict(),
                "structure":          self.structure.as_dict()}

    def site_set_to_yaml_file(self, filename="interstitials.yaml"):
        with open(filename, "w") as f:
            f.write(yaml.dump(self.site_set_as_dict()))

    @property
    def coords(self):
        """Return list of fractional coordinates of interstitial sites"""
        return [v.representative_coords
                for v in self.interstitial_sites.values()]

    @classmethod
    def from_dict(cls, d):
        # orderedDict disables MSONable.
        structure = d["structure"]
        if isinstance(structure, dict):
            structure = Structure.from_dict(structure)

        interstitial_sites = OrderedDict()
        for k, v in d["interstitial_site_set"].items():
            interstitial_sites[k] = InterstitialSite.from_dict(v)

        return cls(structure=structure,
                   interstitial_sites=interstitial_sites)

    @classmethod
    def from_files(cls,
                   structure: Union[str, Structure] = "DPOSCAR",
                   filename="interstitials.yaml"):
        if isinstance(structure, str):
            d = {"structure": Structure.from_file(structure)}
        else:
            d = {"structure": structure}

        with open(filename, "r") as f:
            d["interstitial_site_set"] = yaml.load(f)

        return cls.from_dict(d)

    def add_site(self,
                 coord: list,
                 site_name: str = None,
                 check_neighbor_radius: float = 0.3,
                 force_add: bool = False,
                 symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                 angle_tolerance: float = ANGLE_TOL,
                 method: str = "manual"):
        """

        Note that the symprec must be the same as that used for defect analysis
        to keep the same point group.
        """
        # Check whether other sites are too close to the inserted sites.
        # Construct saturated structure with existing interstitial sites.
        if self.coords:
            saturated_structure, _, symmetry_dataset = \
                create_saturated_interstitial_structure(
                    self.structure, self.coords, dist_tol=symprec)

            cart_coord = \
                saturated_structure.lattice.get_cartesian_coords(coord)
            neighbors = saturated_structure.get_sites_in_sphere(
                pt=cart_coord, r=check_neighbor_radius)

            for n in neighbors:
                # DummySpecie occupies the saturated interstitial sites.
                if isinstance(n[0].specie, DummySpecie):
                    specie = "another interstitial site"
                else:
                    specie = n[0].frac_coords
                distance = n[1]
                message = f"Inserted position is too close to {specie}.\n " \
                          f"The distance is {distance:5.3f} A."
                if force_add:
                    logger.warning(message)
                else:
                    raise StructureError(message)
        else:
            sga = SpacegroupAnalyzer(self.structure, symprec, angle_tolerance)
            symmetry_dataset = sga.get_symmetry_dataset()

        structure = self.structure.copy()
        structure.insert(0, DummySpecie(), coord)

        # check the symmetry of the newly inserted interstitial site
        sga = SpacegroupAnalyzer(structure, symprec, angle_tolerance)
        sym_dataset = sga.get_symmetry_dataset()
        wyckoff = sym_dataset["wyckoffs"][0]
        site_symmetry = sym_dataset["pointgroup"]
        symmetry_multiplicity = \
            get_symmetry_multiplicity(symmetry_dataset, coord,
                                      structure.lattice.matrix, symprec)
        coordination_distances = get_coordination_distances(structure, 0)

        # Set the default of the site_name
        if site_name is None:
            for i in range(MAX_NUM_INTERSTITIAL_SITES):
                trial_name = "i" + str(i + 1)
                if trial_name not in self.interstitial_sites.keys():
                    site_name = trial_name
                    break
            else:
                raise ValueError("Site name is not assigned automatically.")
        else:
            if site_name in self.interstitial_sites.keys():
                raise ValueError(f"Site {site_name} already exists.")

        interstitial_site = \
            InterstitialSite(representative_coords=coord,
                             wyckoff=wyckoff,
                             site_symmetry=site_symmetry,
                             symmetry_multiplicity=symmetry_multiplicity,
                             coordination_distances=coordination_distances,
                             method=method)

        self.interstitial_sites[site_name] = interstitial_site

    def add_sites_from_charge_density(self,
                                      chgcar_filename: str,
                                      trans_mat: list = None,
                                      threshold_frac: float = None,
                                      threshold_abs: float = None,
                                      min_dist: float = 0.5,
                                      tol: float = 0.2,
                                      radius: float = 0.4,
                                      **kwargs):

        chgcar = Chgcar.from_file(chgcar_filename)
        cda = ChargeDensityAnalyzer(chgcar=chgcar)
        cda.get_local_extrema(threshold_frac=threshold_frac,
                              threshold_abs=threshold_abs)
        cda.sort_sites_by_integrated_chg(r=radius)
        cda.remove_collisions(min_dist)
        cda.cluster_nodes(tol=tol)
        print(cda.extrema_coords)
#        print(cda.charge_distribution_df)

        for c in cda.extrema_coords:
            self.add_site(c)
