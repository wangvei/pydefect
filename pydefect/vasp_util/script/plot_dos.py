# -*- coding: utf-8 -*-

from collections import OrderedDict
from pymatgen import Structure
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import Vasprun

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "May 15, 2018"


def get_dos_plot(vasprun_file, site=False, element=False, orbital=False):
    v = Vasprun(vasprun_file)
    dos = v.complete_dos

    all_dos = OrderedDict()
    all_dos["Total"] = dos

    structure = v.final_structure

    if site:
        for i in range(len(structure)):
            site = structure[i]
            all_dos["Site " + str(i) + " " + site.specie.symbol] = \
                dos.get_site_dos(site)

    if element:
        syms = [tok.strip() for tok in element[0].split(",")]
        all_dos = {}
        for el, dos in dos.get_element_dos().items():
            if el.symbol in syms:
                all_dos[el] = dos
    if orbital:
        all_dos = dos.get_spd_dos()

    plotter = DosPlotter()
    plotter.add_dos_dict(all_dos)
    return plotter.get_plot()