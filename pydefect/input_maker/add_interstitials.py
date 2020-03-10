from typing import List

import numpy as np

from pydefect.core.config import DEFECT_SYMMETRY_TOLERANCE, ANGLE_TOL
from pydefect.core.interstitial_site import InterstitialSiteSet
from pydefect.input_maker.defect_initial_setting import DefectInitialSetting

from pymatgen import Structure


def add_interstitials(uc_coords: List[float],
                      vicinage_radius: float,
                      dposcar: str = "DPOSCAR",
                      interstitials_yaml: str = "interstitials.yaml",
                      defect_in_file: str = "defect.in",
                      defect_symprec: float = DEFECT_SYMMETRY_TOLERANCE,
                      angle_tolerance: float = ANGLE_TOL) -> None:
    if len(uc_coords) == 3:
        uc_coords = [uc_coords]
    elif len(uc_coords) % 3 == 0:
        n_coords = int(len(uc_coords) / 3)
        uc_coords = \
            [[uc_coords[3 * i + j] for j in range(3)] for i in range(n_coords)]
    else:
        raise ValueError(f"Interstitial coordinates {uc_coords} in the unit "
                         f"cell are invalid.")

    defect_initial_setting = \
        DefectInitialSetting.from_defect_in(
            poscar=dposcar,
            interstitials_yaml=interstitials_yaml,
            defect_in_file=defect_in_file)

    tm_list = defect_initial_setting.transformation_matrix
    trans_mat = [[tm_list[3 * i + j] for j in range(3)] for i in range(3)]
    inv_trans_mat = np.linalg.inv(trans_mat)
    # inv_trans_mat must be multiplied with coords from the right as the
    # trans_mat is multiplied to the unitcell lattice vector from the left.
    # see __mul__ of IStructure in pymatgen.
    # x_u, x_s means the frac coordinates in unitcell and supercell,
    # while a, b, c are the unitcell lattice vector.
    # (a_u, b_u, c_u) . (a, b, c) = (a_s, b_s, c_s) . trans_mat . (a, b, c)
    # (a_u, b_u, c_u) = (a_s, b_s, c_s) . trans_mat
    # so, (a_s, b_s, c_s) = (a_u, b_u, c_u) . inv_trans_mat
    supercell_coords = [np.dot(c, inv_trans_mat).tolist() for c in uc_coords]

    try:
        interstitial_set = \
            InterstitialSiteSet.from_files(dposcar=dposcar,
                                           yaml_filename=interstitials_yaml)
    except FileNotFoundError:
        structure = Structure.from_file(dposcar)
        interstitial_set = InterstitialSiteSet(structure=structure)

    interstitial_set.add_sites(frac_coords=supercell_coords,
                               vicinage_radius=vicinage_radius,
                               defect_symprec=defect_symprec,
                               angle_tolerance=angle_tolerance)
    interstitial_set.site_set_to_yaml_file(yaml_filename=interstitials_yaml)