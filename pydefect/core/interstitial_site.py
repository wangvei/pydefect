# -*- coding: utf-8 -*-

from collections import OrderedDict
from monty.json import MSONable
import yaml

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


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

    def as_dict(self):
        d = OrderedDict(
            {"site_name":              self.site_name,
             "representative_coords":  self.representative_coords,
             "wyckoff":                self.wyckoff,
             "site_symmetry":          self.site_symmetry,
             "symmetry_multiplicity":  self.symmetry_multiplicity,
             "coordination_distances": self.coordination_distances,
             "method":                 self.method})

        return d


def represent_odict(dumper, instance):
    return dumper.represent_mapping('tag:yaml.org,2002:map', instance.items())


yaml.add_representer(OrderedDict, represent_odict)


def construct_odict(loader, node):
    return OrderedDict(loader.construct_pairs(node))


yaml.add_constructor('tag:yaml.org,2002:map', construct_odict)


class InterstitialSites(list):
    """Holds a set of InterstitialSite objects.
    """

    def __init__(self, interstitial_sites: list):
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
            else:
                name_set.add(i.site_name)

        super().__init__(interstitial_sites)

    def as_dict(self):
        d = OrderedDict()
        for i, interstitial in enumerate(self):
            d[i] = interstitial.as_dict()
        return d

    def to_yaml_file(self, filename="interstitials.yaml"):
        with open(filename, "w") as f:
            f.write(yaml.dump(self.as_dict()))

    @classmethod
    def from_dict(cls, d):
        return cls([InterstitialSite.from_dict(v) for v in d.values()])

    @classmethod
    def from_yaml_file(cls, filename="interstitials.yaml"):
        with open(filename, "r") as f:
            d = yaml.load(f)

        return cls.from_dict(d)

    # def add_interstitial(self, structure, frac_coords, tolerance):
    #     from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        # # use create_saturated_interstitial_structure??
        # sga = SpacegroupAnalyzer(structure)
        # sg_ops = sga.get_symmetry_operations()
        # primitive = sga.find_primitive()



#        self += new_int



    # @classmethod
    # def from_interstitial_in(cls, interstitial_in_file="interstitials.in"):

        # interstitial_sites = []
        # with open(interstitial_in_file) as i:
        #     for l in i:
        #         line = l.split()

                # # No line also needs to comment out now.
                # if not line:
                #     interstitial_sites.append(
                #         InterstitialSite(
                #             site_name=site_name,
                #             representative_coords=representative_coords,
                #             wyckoff=wyckoff,
                #             site_symmetry=site_symmetry,
                #             symmetry_multiplicity=symmetry_multiplicity,
                #             coordination_distances=coordination_distances,
                #             method=method))

                    # site_name = None
                    # representative_coords = None
                    # wyckoff = None
                    # site_symmetry = None
                    # symmetry_multiplicity = None
                    # coordination_distances = None
                    # method = None

                # elif line[0][0] == "#":
                #     continue
                # elif line[0] == "Name:":
                #     site_name = line[1]
                # elif line[0] == "Fractional":
                #     representative_coords = [float(x) for x in line[2:]]
                # elif line[0] == "Wyckoff":
                #     wyckoff = line[2]
                # elif line[0] == "Site":
                #     site_symmetry = line[2]
                # elif line[0] == "Multiplicity:":
                #     symmetry_multiplicity = line[1]
                # elif line[0] == "Coordination:":
                #     coordination_distances = get_distances(line[1:])
                # elif line[0] == "Method:":
                #     method = line[1]
                # else:
                #     raise InvalidFileError(
                #         "{} is not supported in interstitials.in!".format(line))

        # return cls(interstitial_sites)


