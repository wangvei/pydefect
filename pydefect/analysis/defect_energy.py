# -*- coding: utf-8 -*-

from collections import defaultdict, namedtuple
from itertools import combinations
import json

from monty.json import MontyEncoder, MSONable

from obadb.analyzer.chempotdiag.chem_pot_diag import ChemPotDiag

from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.core.defect import Defect

TransitionLevel = namedtuple("TransitionLevel", ("cross_points", "charges"))


class DefectEnergies(MSONable):
    def __init__(self,
                 energies: dict,
                 transition_levels: dict,
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 magnetization: dict,
                 title: str = None):
        """ A class related to a set of defect formation energies.
        Args:
            energies (dict):
                Defect formation energies. energies[name][charge]
            transition_levels (dict):
                key is defect name and value is TransitionLevel.
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
        self.energies = energies
        self.transition_levels = transition_levels
        self.vbm = vbm
        self.cbm = cbm
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
        self.magnetization = magnetization
        self.title = title

    @classmethod
    def from_files(cls,
                   unitcell: UnitcellCalcResults,
                   perfect: SupercellCalcResults,
                   defects: list,
                   chem_pot: ChemPotDiag,
                   chem_pot_label: str,
                   system: str = ""):
        """ Calculates defect formation energies from several objects.

        Note that all the energies are calculated at 0 eV in the absolute scale.
        Args:
            unitcell (UnitcellCalcResults):
                UnitcellCalcResults object for band edge.
            perfect (SupercellCalcResults):
                SupercellDftResults object of perfect supercell for band edge in
                supercell.
            defects (list of namedtuple Defect):
                List of the Defect namedtuple object.
                Defect = namedtuple(
                    "Defect", "defect_entry", "dft_results", "correction")
            chem_pot (ChemPot):
                Chemical potentials of the competing phases.
            chem_pot_label (str):
                Equilibrium point specified in ChemPot.
            system (str):
                System name used for the title.
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
            element_diff = d.defect_entry.changes_of_num_elements

            # calculate four terms for a defect formation energy.
            relative_energy = d.dft_results.relative_total_energy(perfect)
            correction_energy = d.correction.correction_energy
            element_interchange_energy = \
                - sum([v * (relative_chem_pot.elem_coords[k] + standard_e[k])
                       for k, v in element_diff.items()])

            energies[name][charge] = \
                relative_energy + correction_energy + element_interchange_energy

            magnetization[name][charge] = d.dft_results.total_magnetization

        transition_levels = {}

        # e_of_c means energy as a function of charge: e_of_c[charge] = energy
        for name, e_of_c in energies.items():
            cross_points = []
            charge = []

            for (c1, e1), (c2, e2) in combinations(e_of_c.items(), r=2):
                # The cross point between two charge states.
                x = - (e1 - e2) / (c1 - c2)
                y = (c1 * e2 - c2 * e1) / (c1 - c2)

                # The lowest energy among all the charge states to be compared.
                compared_energy = \
                    min([energy + c * x for c, energy in e_of_c.items()])

                if y < compared_energy + 1e-5:
                    cross_points.append([x, y])
                    charge.append([c1, c2])

            transition_levels[name] = \
                TransitionLevel(
                    cross_points=sorted(cross_points, key=lambda z: z[0]),
                    charges=sorted(charge))

        return cls(energies, transition_levels, vbm, cbm, supercell_vbm,
                   supercell_cbm, magnetization, title)

    def to_json_file(self, filename="defect_energy.json"):
        """ Returns a json file. """
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def __str__(self):
        pass

    def u(self, name, charge):
        """ Return the U value among three charge states.

        Args:
            name (str):
                Name of the defect.
            charge (list):
                1x3 list comprising three charge states.
        """
        if len(charge) != 3:
            assert ValueError("The length of charge states must be 3.")
        elif charge[0] + charge[2] - 2 * charge[1] != 0:
            assert ValueError("The charge states {} {} {} are not balanced."
                              .format(*charge))

        energies = [self.energies[name][c] for c in charge]
        return energies[0] + energies[2] - 2 * energies[1]

    @property
    def band_gap(self):
        return self.cbm - self.vbm
