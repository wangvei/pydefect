from collections import defaultdict
from typing import Optional, Callable, Union
from xml.etree.ElementTree import ParseError

import numpy as np
from pydefect.util.logger import get_logger
from pymatgen import Spin

logger = get_logger(__name__)


def is_str_digit(n: str) -> bool:
    try:
        float(n)
        return True
    except ValueError:
        return False


def is_str_int(n: str) -> bool:
    try:
        if int(n) - float(n) < 1e-5:
            return True
        else:
            return False
    except ValueError:
        return False


def spin_key_to_str(arg: Optional[dict],
                    value_to_str=False) -> Optional[dict]:
    """ Change Spin object in key to int.

    Optionally, value is also changed to str.
    """
    if arg is not None:
        if value_to_str:
            return {str(spin): str(v) for spin, v in arg.items()}
        else:
            return {str(spin): v for spin, v in arg.items()}
    else:
        return


def str_key_to_spin(arg: Optional[dict],
                    method_from_str_for_value=None) -> Optional[dict]:
    """ Change int key to Spin object. """
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


def parse_file(classmethod_name: Callable, parsed_filename: str):
    """Check filename and parse and return cls via __init__ or classmethod """
    try:
        logger.info("Parsing {}...".format(parsed_filename))
        return classmethod_name(parsed_filename)
    except ParseError:
        logger.warning("Parsing {} failed.".format(parsed_filename))
        raise ParseError
    except FileNotFoundError:
        logger.warning("File {} doesn't exist.".format(parsed_filename))
        raise FileNotFoundError


def defaultdict_to_dict(d: dict) -> dict:
    """Recursively change defaultdict to dict"""
    if isinstance(d, defaultdict):
        d = dict(d)
    if isinstance(d, dict):
        for key, value in d.items():
            d[key] = defaultdict_to_dict(value)

    return d


def make_symmetric_matrix(d: Union[list, float]) -> np.ndarray:
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


def sanitize_keys_in_dict(d: dict) -> dict:
    """ Recursively sanitize keys in dict from str to int, float and None.
    Args
        d (dict):
            d[name][charge][annotation]
    """
    if not isinstance(d, dict):
        return d
    else:
        new_d = dict()
        for key, value in d.items():
            try:
                key = int(key)
            except (ValueError, TypeError):
                try:
                    key = float(key)
                except (ValueError, TypeError):
                    if key == "null":
                        key = None
            value = None if value == "null" else sanitize_keys_in_dict(value)
            new_d[key] = value
        return new_d


def construct_obj_in_dict(d: dict, cls: Callable) -> dict:
    """ Recursively sanitize keys in dict from str to int, float and None.
    Args
        d (dict):
            d[name][charge][annotation]
    """
    from copy import deepcopy
    if not isinstance(d, dict):
        return d
    else:
        new_d = deepcopy(d)
        for key, value in d.items():
            if value.get("@class", "") == cls.__name__:
                new_d[key] = cls.from_dict(value)
            else:
                new_d[key] = construct_obj_in_dict(value, cls)
        return new_d


def flatten_dict(d: dict, depth: int = None) -> list:
    """ Flatten keys and values in dic

    d[a][b][c] = x -> [a, b, c, x]
    """
    flattened_list = list()
    for key, value in d.items():
        if isinstance(value, dict) and depth != 1:
            depth = depth - 1 if depth else None
            flattened_list.extend(
                [[key] + v for v in flatten_dict(value, depth)])
        else:
            flattened_list.append([key, value])

    return flattened_list
