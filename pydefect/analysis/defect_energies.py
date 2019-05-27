# -*- coding: utf-8 -*-

from collections import defaultdict, namedtuple
import json

from monty.json import MontyEncoder, MSONable

from obadb.analyzer.chempotdiag.chem_pot_diag import ChemPotDiag

from pydefect.core.supercell_calc_results import SupercellCalcResults
from pydefect.core.unitcell_calc_results import UnitcellCalcResults
from pydefect.database.num_symmetry_operation import num_symmetry_operation

TransitionLevel = namedtuple("TransitionLevel", ("cross_points", "charges"))
Defect = namedtuple("Defect", ("defect_entry", "dft_results", "correction"))


class DefectEnergies(MSONable):
    def __init__(self,
                 energies: dict,
                 multiplicity: dict,
                 magnetization: dict,
                 convergence: dict,
                 shallow: dict,
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 title: str = None):
        """ A class related to a set of defect formation energies.

        Args:
            energies (dict):
                Defect formation energies. energies[name][charge]
            multiplicity (dict):
                Multiplicity. multiplicity[name][charge]
            magnetization (dict):
                Magnetization in mu_B. magnetization[name][charge]
            convergence (bool):
                Whether to remove the unconverged defects.
            shallow (bool):
                Whether to remove the shallow defects.
            vbm (float):
                Valence band maximum in the unitcell in the absolute scale.
            cbm (float):
                Conduction band minimum in the unitcell in the absolute scale.
            supercell_vbm (float):
                Valence band maximum in the perfect supercell.
            supercell_cbm (float):
                Conduction band minimum in the perfect supercell.
            title (str):
                Title of the system.
        """
        self.energies = energies
        self.multiplicity = multiplicity
        self.magnetization = magnetization
        self.convergence = convergence
        self.shallow = shallow
        self.vbm = vbm
        self.cbm = cbm
        self.supercell_vbm = supercell_vbm
        self.supercell_cbm = supercell_cbm
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
        multiplicity = defaultdict(dict)
        magnetization = defaultdict(dict)
        convergence = defaultdict(dict)
        shallow = defaultdict(dict)

        for d in defects:
            name = d.defect_entry.name
            charge = d.defect_entry.charge

            element_diff = d.defect_entry.changes_of_num_elements

            # calculate four terms for a defect formation energy.
            relative_energy = d.dft_results.relative_total_energy
            correction_energy = d.correction.correction_energy
            element_interchange_energy = \
                - sum([v * (relative_chem_pot.elem_coords[k] + standard_e[k])
                       for k, v in element_diff.items()])

            energies[name][charge] = \
                relative_energy + correction_energy + element_interchange_energy

            initial_num_symmops = \
                num_symmetry_operation(d.defect_entry.initial_site_symmetry)
            final_num_symmops = \
                num_symmetry_operation(d.dft_results.site_symmetry)

            multiplicity[name][charge] = \
                int(d.defect_entry.num_equiv_sites / final_num_symmops
                    * initial_num_symmops)

            magnetization[name][charge] = d.dft_results.total_magnetization
            convergence[name][charge] = d.dft_results.is_converged
            shallow[name][charge] = d.dft_results.shallow

        return cls(dict(energies), dict(multiplicity), dict(magnetization),
                   dict(convergence), dict(shallow), vbm, cbm,
                   supercell_vbm, supercell_cbm, title)

    def to_json_file(self, filename="defect_energy.json"):
        with open(filename, 'w') as fw:
            json.dump(self.as_dict(), fw, indent=2, cls=MontyEncoder)

    def __str__(self):
        outs = []
        for n in self.energies.keys():
            for c in self.energies[n].keys():
                outs.extend(
                    ["- name: {} - charge: {}".format(n, c),
                     "Energy at e_f = 0 (eV): {}".format(
                         round(self.energies[n][c], 4)),
                     "Multiplicity: {}".format(self.multiplicity[n][c]),
                     "Magnetization: {}".format(
                         round(self.magnetization[n][c], 2)),
                     "Convergence: {}".format(self.convergence[n][c]),
                     "Shallow: {}".format(self.shallow[n][c])])
            outs.append("")

        return "\n".join(outs)

    def u(self, name, charge):
        """ Return the U value among three charge states.

        Args:
            name (str):
                Name of the defect.
            charge (list):
                1x3 list comprising three charge states.

        Return:
            U value.
        """
        if len(charge) != 3:
            assert ValueError("The length of charge states must be 3.")
        elif charge[2] - charge[1] != 1 or charge[1] - charge[0] != 1:
            assert ValueError("The charge states {} {} {} are not sequential."
                              .format(*charge))

        energies = [self.energies[name][c] for c in charge]
        return energies[0] + energies[2] - 2 * energies[1]

    @property
    def band_gap(self):
        return self.cbm - self.vbm
