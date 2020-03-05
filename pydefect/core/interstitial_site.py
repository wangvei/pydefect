# -*- coding: utf-8 -*-

from collections import OrderedDict
from copy import deepcopy
from typing import Union, Dict, List
import yaml

from monty.json import MSONable

from pydefect.core.config import (
    DEFECT_SYMMETRY_TOLERANCE, INTERSTITIAL_SYMPREC, ANGLE_TOL, CUTOFF_FACTOR)
from pydefect.database.symmetry import num_symmetry_operation
from pydefect.util.logger import get_logger
from pydefect.util.structure_tools import (
    create_saturated_interstitial_structure, get_coordination_distances,
    min_distance_from_coords)

from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


logger = get_logger(__name__)


class InterstitialSite(MSONable):
    """Holds properties related to the interstitial site. """

    def __init__(self,
                 representative_coords: list,
                 wyckoff: str = None,
                 site_symmetry: str = None,
                 multiplicity: int = None,
                 coordination_distances: dict = None,
                 cutoff: float = None,
                 method: str = None):
        """
        Args:
            representative_coords (list):
                Representative coordinates, namely the position of first_index
            wyckoff (str):
                A wyckoff letter.
            site_symmetry (str):
                Site symmetry in Hermann–Mauguin notation.
            multiplicity (int):
                Multiplicity of the interstitial in supercell perfect structure.
            coordination_distances (dict):
                Coordination environment. An example is
                {"Mg": [1.92, 1.95, 2.01], "Al": [1.82, 1.95]}
            method (str):
                The method name determining the interstitial site.
        """

        self.representative_coords = representative_coords
        self.wyckoff = wyckoff
        self.site_symmetry = site_symmetry
        self.multiplicity = multiplicity
        self.coordination_distances = coordination_distances
        self.cutoff = float(cutoff)
        self.method = method

    def __repr__(self):
        outs = [f"representative coords: {self.representative_coords}",
                f"wyckoff: {self.wyckoff}",
                f"site symmetry: {self.site_symmetry}",
                f"multiplicity: {self.multiplicity}",
                f"coordination distances: {self.coordination_distances}",
                f"cutoff: {self.cutoff}",
                f"method: {self.method}"]
        return "\n".join(outs)

    def as_dict(self):
        d = OrderedDict(
            {"representative_coords":  self.representative_coords,
             "wyckoff":                self.wyckoff,
             "site_symmetry":          self.site_symmetry,
             "multiplicity":           self.multiplicity,
             "coordination_distances": self.coordination_distances,
             "cutoff":                 self.cutoff,
             "method":                 self.method})

        return d


# The followings are needed to keep order of dictionary for interstitial.yaml.
# https://qiita.com/podhmo/items/aa954ee1dc1747252436
# If yaml files show weird output, check quantities does not include numpy
# variables (e.g., numpy.float)
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

    def as_dict(self) -> Dict[str, dict]:
        return {"interstitial_site_set": self.site_set_as_dict(),
                "structure":          self.structure.as_dict()}

    def site_set_to_yaml_file(self, yaml_filename="interstitials.yaml") -> None:
        with open(yaml_filename, "w") as f:
            f.write(yaml.dump(self.site_set_as_dict()))

    @property
    def coords(self) -> List[list]:
        """Return list of fractional coordinates of interstitial sites"""
        return \
            [v.representative_coords for v in self.interstitial_sites.values()]

    @classmethod
    def from_dict(cls, d: dict):
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
                   dposcar: Union[str, Structure] = "DPOSCAR",
                   yaml_filename="interstitials.yaml"):
        if isinstance(dposcar, str):
            d = {"structure": Structure.from_file(dposcar)}
        else:
            d = {"structure": dposcar}

        with open(yaml_filename, "r") as f:
            d["interstitial_site_set"] = yaml.load(f, Loader=yaml.FullLoader)

        return cls.from_dict(d)

    def add_sites(self,
                  frac_coords: list,
                  vicinage_radius: float = 0.3,
                  defect_symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                  angle_tolerance: float = ANGLE_TOL,
                  method: str = "manual") -> None:
        """ Add interstitial sites

        Note that symprec must be the same as that used for
        DefectInitialSetting to keep the symmetry consistency such as
        point group, multiplicity and so on.

        Args:
            frac_coords (list):
                Added fractional coords in the self.structure.
            vicinage_radius (float):
                Radius in which atoms are considered too close.
            defect_symprec (float):
                Precision used for symmetry analysis in angstrom.
            angle_tolerance (float):
                Angle tolerance for symmetry analysis in degree
            method (str):
                Name of the method that determine the interstitial sites.
        """
        # Don't distinguish old and new interstitial coordinates for simplicity
        total_coords = self.coords + frac_coords
        saturated_structure, atom_indices, are_inserted = \
            create_saturated_interstitial_structure(
                structure=self.structure,
                inserted_atom_coords=total_coords,
                dist_tol=vicinage_radius,
                symprec=defect_symprec,
                angle_tolerance=angle_tolerance)

        min_dist = min_distance_from_coords(self.structure, total_coords)
        cutoff = round(min_dist * CUTOFF_FACTOR, 2)

        sga = SpacegroupAnalyzer(saturated_structure, defect_symprec,
                                 angle_tolerance)
        symmetry_dataset = sga.get_symmetry_dataset()
        num_symmop = len(sga.get_symmetry_operations())

        i = 1
        for (coord, ai, inserted) in \
                zip(frac_coords, atom_indices, are_inserted):
            if not inserted:
                continue
            while True:
                if "i" + str(i) in self.interstitial_sites:
                    i += 1
                else:
                    site_name = "i" + str(i)
                    break
            site_symmetry = symmetry_dataset["site_symmetry_symbols"][ai]
            wyckoff = symmetry_dataset["wyckoffs"][ai]
            # multiplicity is reduced by the number of symmetry operation
            multiplicity = round(int(num_symmop
                                 / num_symmetry_operation(site_symmetry)))
            # Calculate the coordination_distances.
            # Note that saturated_structure includes other interstitial sites.
            defect_str = self.structure.copy()
            defect_str.append(DummySpecie(), coord)
            coordination_distances = \
                get_coordination_distances(structure=defect_str,
                                           atom_index=len(defect_str) - 1,
                                           cutoff=cutoff)

            interstitial_site = \
                InterstitialSite(representative_coords=coord,
                                 wyckoff=wyckoff,
                                 site_symmetry=site_symmetry,
                                 multiplicity=multiplicity,
                                 coordination_distances=coordination_distances,
                                 cutoff=cutoff,
                                 method=method)

            self.interstitial_sites[site_name] = interstitial_site


def interstitials_from_charge_density(
        chgcar_filename: str,
        threshold_frac: float = None,
        threshold_abs: float = None,
        min_dist: float = 0.5,
        tol: float = 0.2,
        radius: float = 0.4,
        interstitial_symprec: float = INTERSTITIAL_SYMPREC,
        angle_tolerance: float = ANGLE_TOL) -> None:
    """ Print interstitial sites determined from charge density local minimum

    Args:
        chgcar_filename (str):
            CHGCAR-type file name.
        threshold_frac:
            See docstrings of ChargeDensityAnalyzer.
        threshold_abs:
            See docstrings of ChargeDensityAnalyzer.
        min_dist:
            See docstrings of ChargeDensityAnalyzer.
        tol:
            See docstrings of ChargeDensityAnalyzer.
        radius:
            See docstrings of ChargeDensityAnalyzer.
        interstitial_symprec (float):
            Precision used for symmetry analysis in angstrom.
        angle_tolerance (float):
            Angle tolerance for symmetry analysis in degree
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
    sga = SpacegroupAnalyzer(structure, interstitial_symprec, angle_tolerance)
    sym_db = sga.get_symmetry_dataset()
    equiv_atoms = sym_db["equivalent_atoms"]

    sym_struct = sga.get_symmetrized_structure()

    print("")
    print("++ Inequivalent indices and site symmetries ++")
    orig_num_atoms = len(chgcar.structure)
    for i, ii in enumerate(interstitial_indices):
        if ii == equiv_atoms[ii]:
            idx = orig_num_atoms + i
            coords = sym_struct[idx].frac_coords
            idx_coords = \
                f"{i:>3} {coords[0]:8.4f} {coords[1]:8.4f} {coords[2]:8.4f}"

            print(idx_coords, sym_db["site_symmetry_symbols"][ii])

    # Change coords from unitcell to supercell
    # multiply inverse of trans_mat to coords
    # inv_trans_mat = np.linalg.inv(trans_mat)
    # print(inv_trans_mat)
    # supercell_coords = \
    #     [np.dot(inv_trans_mat, c) for c in inequiv_interstitial_coords]
    # print(supercell_coords)

    # self.add_sites(supercell_coords, symprec=symprec, angle_tol=angle_tol,
    #                method="charge", **kwargs)
