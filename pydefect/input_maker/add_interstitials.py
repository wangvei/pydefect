from typing import List

import numpy as np

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting

from pymatgen import Structure


def add_interstitials(coords_in_unitcell: List[float],
                      vicinage_radius: float,
                      uposcar: str = "UPOSCAR",
                      interstitials_yaml: str = "interstitials.yaml",
                      defect_symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                      angle_tolerance: float = ANGLE_TOL) -> None:
    if len(coords_in_unitcell) == 3:
        coords_in_unitcell = [coords_in_unitcell]
    elif len(coords_in_unitcell) % 3 == 0:
        n_coords = int(len(coords_in_unitcell) / 3)
        coords_in_unitcell = [[coords_in_unitcell[3 * i + j]
                               for j in range(3)] for i in range(n_coords)]
    else:
        raise ValueError(f"Interstitial coordinates {coords_in_unitcell} "
                         f"in the unit cell are invalid.")

    try:
        interstitial_set = \
            InterstitialSiteSet.from_files(uposcar=uposcar,
                                           yaml_filename=interstitials_yaml)
    except FileNotFoundError:
        structure = Structure.from_file(uposcar)
        interstitial_set = InterstitialSiteSet(structure=structure)

    interstitial_set.add_sites(added_coords=coords_in_unitcell,
                               vicinage_radius=vicinage_radius,
                               defect_symprec=defect_symprec,
                               angle_tolerance=angle_tolerance)
    interstitial_set.site_set_to_yaml_file(yaml_filename=interstitials_yaml)
