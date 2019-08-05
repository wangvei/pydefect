from collections import defaultdict
from xml.etree.ElementTree import ParseError

import numpy as np
from pydefect.util.logger import get_logger
from pymatgen import Spin


logger = get_logger(__name__)


def spin_key_to_str(arg: dict, value_to_str=False):
    if arg is not None:
        if value_to_str:
            return {str(spin): str(v) for spin, v in arg.items()}
        else:
            return {str(spin): v for spin, v in arg.items()}
    else:
        return


def str_key_to_spin(arg: dict, method_from_str_for_value=None):
    if arg is not None:
        x = {}
        for spin, value in arg.items():
            x[Spin(int(spin))] = value
            if method_from_str_for_value:
                x[Spin(int(spin))] = method_from_str_for_value(value)
            else:
                x[Spin(int(spin))] = value
        return x
    else:
        return


def parse_file(classmethod_name, parsed_filename):
    try:
        logger.info("Parsing {}...".format(parsed_filename))
        return classmethod_name(parsed_filename)
    except ParseError:
        logger.warning("Parsing {} failed.".format(parsed_filename))
        raise ParseError
    except FileNotFoundError:
        logger.warning("File {} doesn't exist.".format(parsed_filename))
        raise FileNotFoundError


def defaultdict_to_dict(d):
    """Recursively change defaultdict to dict"""
    if isinstance(d, defaultdict):
        d = dict(d)
    if isinstance(d, dict):
        for key, value in d.items():
            d[key] = defaultdict_to_dict(value)

    return d


def make_symmetric_matrix(d):
    """
    d (list or float):
        len(d) == 1: Suppose cubic system
        len(d) == 3: Suppose tetragonal or orthorhombic system
        len(d) == 6: Suppose the other system
    """
    if isinstance(d, float):
        tensor = np.array([[d, 0, 0],
                           [0, d, 0],
                           [0, 0, d]])
    elif len(d) == 1:
        tensor = np.array([[d[0], 0,  0],
                           [0,  d[0], 0],
                           [0,  0,  d[0]]])
    elif len(d) == 3:
        tensor = np.array([[d[0], 0, 0],
                           [0, d[1], 0],
                           [0, 0, d[2]]])
    elif len(d) == 6:
        from pymatgen.util.num import make_symmetric_matrix_from_upper_tri
        """ 
        Given a symmetric matrix in upper triangular matrix form as flat array 
        indexes as:
        [A_xx, A_yy, A_zz, A_xy, A_xz, A_yz]
        This will generate the full matrix:
        [[A_xx, A_xy, A_xz], [A_xy, A_yy, A_yz], [A_xz, A_yz, A_zz]
        """
        tensor = make_symmetric_matrix_from_upper_tri(d)
    else:
        raise ValueError("{} is not valid to make symmetric matrix".format(d))

    return tensor