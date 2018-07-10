# -*- coding: utf-8 -*-

import json

from monty.json import MontyEncoder
from monty.serialization import loadfn

from collections import defaultdict, namedtuple
from itertools import combinations

from pydefect.core.supercell_dft_results import SupercellDftResults
from pydefect.core.unitcell_dft_results import UnitcellDftResults

Defect = namedtuple("Defect", ("defect_entry", "dft_results", "correction"))
TransitionLevel = namedtuple("TransitionLevel", ("cross_points", "charges"))


class DefectEnergies:
    """
    A class related to a set of defect formation energies.
    """
    def __init__(self, energies, transition_levels, vbm, cbm,
                 supercell_vbm, supercell_cbm, magnetization, title=None):
        """
        Args:
            energies (defaultdict):
                Defect formation energies. energies[name][charge]
            transition_levels (dict):
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            supercell_vbm (float):
                Valence band maximum in the perfect supercell.
            supercell_cbm (float):
                Conduction band minimum in the perfect supercell.
            magnetization (dict):
                Magnetization in \mu_B. magnetization[defect][charge]
            title (str):
                Title of the system.
        """
        self._energies = energies
        self._transition_levels = transition_levels
        self._vbm = vbm
        self._cbm = cbm
        self._supercell_vbm = supercell_vbm
        self._supercell_cbm = supercell_cbm
        self._magnetization = magnetization
        self._title = title

    @classmethod
    def from_files(cls, unitcell, perfect, defects, chem_pot, chem_pot_label,
                   system=""):
        """
        Calculates defect formation energies.
        Args:
            unitcell (UnitcellDftResults):
            perfect (SupercellDftResults)
            defects (list of namedtuple Defect):
                [[Defect, ...]
                Defect = namedtuple("Defect", "defect_entry", "dft_results",
                                    "correction")
            chem_pot (ChemPot): Chempot class object.
            chem_pot_label (str):
            system (str): A system name used for the title
        """
        # Note: vbm, cbm, perfect_vbm, perfect_cbm are in absolute energy.
        vbm, cbm = unitcell.band_edge
        supercell_cbm, supercell_vbm = perfect.eigenvalue_properties[1:3]

        title = system + " condition " + chem_pot_label

        # Chemical potentials
        relative_chem_pots, standard_e = chem_pot
        relative_chem_pot = relative_chem_pots[chem_pot_label]

        # Calculate defect formation energies at the vbm
        energies = defaultdict(dict)
        magnetization = defaultdict(dict)

        for d in defects:
            name = d.defect_entry.name
            charge = d.defect_entry.charge
            element_diff = d.defect_entry.element_diff

            # calculate four terms for a defect formation energy.
            relative_energy = d.dft_results.relative_total_energy(perfect)
            correction_energy = d.correction.total_correction_energy
            electron_interchange_energy = vbm * charge
            element_interchange_energy = \
                - sum([v * (relative_chem_pot.elem_coords[k] + standard_e[k])
                       for k, v in element_diff.items()])

            energies[name][charge] = \
                relative_energy + correction_energy + \
                electron_interchange_energy + element_interchange_energy

            magnetization[name][charge] = d.dft_results.magnetization

        # Calculates transition levels
        transition_levels = {}

        for name, e_of_c in energies.items():
            points = []
            charge = set()

            for (c1, e1), (c2, e2) in combinations(e_of_c.items(), r=2):
                # Estimate the cross point between two charge states
                x = - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                compared_energy = cls.min_e_at_ef(e_of_c, x)[1] + 0.00001

                # if x_min < x < x_max and y < compared_energy:
                if y < compared_energy:
                    points.append([x, y])
                    charge.add(c1)
                    charge.add(c2)

            transition_levels[name] = \
                TransitionLevel(cross_points=sorted(points, key=lambda x: x[0]),
                                charges=sorted(charge))
        return cls(energies, transition_levels, vbm, cbm, supercell_vbm,
                   supercell_cbm, magnetization, title)

    @staticmethod
    def min_e_at_ef(ec, ef):
        # d[c1] = energy for charge c1
        d = {c: e + c * ef for c, e in ec.items()}
        # return charge for the lowest energy, and its energy value
        return min(d.items(), key=lambda x: x[1])

    @classmethod
    def from_dict(cls, d):
        """
        Constructs a class object from a dictionary.
        """
        # Programmatic access to enumeration members in Enum class.
        return cls(d["defect_energies"], d["transition_levels"], d["vbm"],
                   d["cbm"], d["supercell_vbm"], d["supercell_cbm"],
                   d["title"])

    @classmethod
    def load_json(cls, filename):
        """
        Constructs a class object from a json file.
        """
        d = loadfn(filename)
        return cls.from_dict(d)

    def as_dict(self):
        """
        Dict representation of the class object.
        """

        d = {"defect_energies":   self._defect_energies,
             "transition_levels": self._transition_levels,
             "vbm":               self._vbm,
             "cbm":               self._cbm,
             "supercell_vbm":     self._supercell_vbm,
             "supercell_cbm":     self._supercell_cbm,
             "title":             self._title}
        return d

    def to_json_file(self, filename):
        """
        Returns a json file, named dft_results.json.
        """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def __str__(self):
        pass

    def U_values(self):
        pass

    @property
    def energies(self):
        return self._energies

    @property
    def transition_levels(self):
        return self._transition_levels

    @property
    def vbm(self):
        return self._vbm

    @property
    def cbm(self):
        return self._cbm

    @property
    def supercell_vbm(self):
        return self._supercell_vbm

    @property
    def supercell_cbm(self):
        return self._supercell_cbm

    @property
    def title(self):
        return self._title

    @property
    def band_gap(self):
        return self._cbm - self._vbm

    @property
    def magnetization(self):
        return self._magnetization

