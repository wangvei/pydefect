# -*- coding: utf-8 -*-

from collections import Iterable
from monty.json import MSONable

from pydefect.core.error_classes import InvalidFileError
from pydefect.input_maker.defect_initial_setting import get_distances

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class InterstitialSite(MSONable):
    """Holds properties related to the interstitial site.

    Args:
        site_name (str):
            Interstitial site name.
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
                 site_name: str,
                 representative_coords: list,
                 wyckoff: str = None,
                 site_symmetry: str = None,
                 symmetry_multiplicity: int = None,
                 coordination_distances: dict = None,
                 method: str = None):
        self.site_name = site_name
        self.representative_coords = representative_coords
        self.wyckoff = wyckoff
        self.site_symmetry = site_symmetry
        self.symmetry_multiplicity = symmetry_multiplicity
        self.coordination_distances = coordination_distances
        self.method = method


class InterstitialSites(MSONable):
    """Holds a set of InterstitialSite objects.
    """

    def __init__(self,
                 interstitial_sites: Iterable):
        """
        Args:
            interstitial_sites (Iterable):
                List of InterstitialSite objects
        """
        name_set = set()
        for i in interstitial_sites:
            if i.site_name in name_set:
                raise ValueError("Interstitial name {} conflicts.".
                                 format(i.site_name))
            elif i.site_name[0] != "i":
                raise ValueError("Interstitial name {} needs to begin *i*.".
                                 format(i.site_name))
            else:
                name_set.add(i.site_name)

        self.interstitial_sites = list(interstitial_sites)

    @classmethod
    def from_interstitial_in(cls, interstitial_in_file="interstitial.in"):

        interstitial_sites = []
        with open(interstitial_in_file) as i:
            for l in i:

                site_name = ""
                representative_coords = []
                wyckoff = None
                site_symmetry = None
                symmetry_multiplicity = None
                coordination_distances = None
                method = None

                line = l.split()

                # No line also needs to comment out now.
                if line[0][0] == "#":
                    continue
                elif line[0] == "Name":
                    site_name = line[1]
                elif line[0] == "Fractional":
                    representative_coords = [float(x) for x in i[2:]]
                elif line[0] == "Wyckoff":
                    wyckoff = line[2]
                elif line[0] == "Site":
                    site_symmetry = line[2]
                elif line[0] == "Multiplicity:":
                    symmetry_multiplicity = line[1]
                elif line[0] == "Coordination:":
                    coordination_distances = get_distances(line[1:])
                elif line[0] == "Method:":
                    method = line[1]
                elif not line:
                    if not site_name or not representative_coords:
                        raise InvalidFileError(
                            "{} file is invalid!".format(interstitial_in_file))
                    else:
                        interstitial_sites.append(
                            InterstitialSite(
                                site_name=site_name,
                                representative_coords=representative_coords,
                                wyckoff=wyckoff,
                                site_symmetry=site_symmetry,
                                symmetry_multiplicity=symmetry_multiplicity,
                                coordination_distances=coordination_distances,
                                method=method))
                else:
                    raise InvalidFileError(
                        "{} is not supported in interstitials.in!".format(line))

        return cls(interstitial_sites)
