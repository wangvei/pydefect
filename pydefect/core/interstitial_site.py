# -*- coding: utf-8 -*-

from collections import OrderedDict
from copy import deepcopy
from typing import Union

import yaml
from monty.json import MSONable
from vise.util.structure_handler import get_coordination_distances
from pydefect.core.config import SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.database.num_symmetry_operation import num_symmetry_operation
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools \
    import create_saturated_interstitial_structure
from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)


class InterstitialSite(MSONable):
    """ Holds properties related to the interstitial site.

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

    def add_sites(self,
                  coords: list,
                  vicinage_radius: float = 0.3,
                  symprec: float = SYMMETRY_TOLERANCE,
                  angle_tol: float = ANGLE_TOL,
                  method: str = "manual"):
        """ Add interstitial sites

        Note that symprec must be the same as that used for
        DefectInitialSetting to keep the symmetry consistency such as
        point group, multiplicity and so on.
        """
        # Don't distinguish old and new interstitial coordinates for simplicity
        total_coords = self.coords + coords
        saturated_structure, atom_indices, are_inserted = \
            create_saturated_interstitial_structure(
              self.structure, total_coords, vicinage_radius, symprec, angle_tol)

        symmetry_dataset = SpacegroupAnalyzer(
            saturated_structure, symprec, angle_tol).get_symmetry_dataset()

        for i, (coord, ai, inserted) \
                in enumerate(zip(coords, atom_indices, are_inserted)):
            if not inserted:
                continue
            site_symmetry = symmetry_dataset["site_symmetry_symbols"][ai]
            wyckoff = symmetry_dataset["wyckoffs"][ai]
            multiplicity = num_symmetry_operation(site_symmetry)
            # Calculate the coordination_distances.
            # Note that saturated_structure includes other interstitial sites.
            defect_str = self.structure.copy()
            defect_str.append(DummySpecie(), coord)
            coordination_distances = \
                get_coordination_distances(defect_str, len(defect_str) - 1)

            interstitial_site = \
                InterstitialSite(representative_coords=coord,
                                 wyckoff=wyckoff,
                                 site_symmetry=site_symmetry,
                                 symmetry_multiplicity=multiplicity,
                                 coordination_distances=coordination_distances,
                                 method=method)

            self.interstitial_sites["i" + str(i)] = interstitial_site


def interstitials_from_charge_density(
        chgcar_filename: str,
        threshold_frac: float = None,
        threshold_abs: float = None,
        min_dist: float = 0.5,
        tol: float = 0.2,
        radius: float = 0.4,
        symprec: float = SYMMETRY_TOLERANCE,
        angle_tol: float = ANGLE_TOL):
    """ Print interstitial sites determined from charge density local minimum

    Note that symprec must be the same as that used for
    DefectInitialSetting to keep the symmetry consistency such as
    point group, multiplicity and so on.
    """

    chgcar = Chgcar.from_file(chgcar_filename)
    cda = ChargeDensityAnalyzer(chgcar=chgcar)
    cda.get_local_extrema(threshold_frac=threshold_frac,
                          threshold_abs=threshold_abs)
    cda.sort_sites_by_integrated_chg(r=radius)
    # Remove sites near host atoms.
    cda.remove_collisions(min_dist)
    # Cluster interstitials that are too close together using a tol.
    cda.cluster_nodes(tol=tol)

    structure = chgcar.structure.copy()
    print(cda.extrema_df)

    start_index = len(structure)
    end_index = len(structure) + len(cda.extrema_coords)
    interstitial_indices = [i for i in range(start_index, end_index)]
    coords = cda.extrema_coords
    for c in coords:
        structure.append(DummySpecie(), c)
    sga = SpacegroupAnalyzer(structure, symprec, angle_tol)
    sym_db = sga.get_symmetry_dataset()
    equiv_atoms = sym_db["equivalent_atoms"]

    print("")
    print(f"Inequivalent indices:")
    for i, ii in enumerate(interstitial_indices):
        if ii == equiv_atoms[ii]:
            print(i, sym_db["site_symmetry_symbols"][ii])

    # Change coords from unitcell to supercell
    # multiply inverse of trans_mat to coords
    # inv_trans_mat = np.linalg.inv(trans_mat)
    # print(inv_trans_mat)
    # supercell_coords = \
    #     [np.dot(inv_trans_mat, c) for c in inequiv_interstitial_coords]
    # print(supercell_coords)

    # self.add_sites(supercell_coords, symprec=symprec, angle_tol=angle_tol,
    #                method="charge", **kwargs)
