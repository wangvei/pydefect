# -*- coding: utf-8 -*-
import re
from inspect import signature, _empty
from os.path import join
from pathlib import Path
from typing import Callable, Optional, List

import yaml
from pydefect.util.logger import get_logger
from pydefect.util.tools import is_str_digit, is_str_int

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

logger = get_logger(__name__)

setting_keys = ["symprec",
                "defect_symprec",
                "angle_tolerance",
                "kpt_density",
                "defect_kpt_density",
                "perfect_incar_setting",
                "defect_incar_setting",
                "perfect_vos_kwargs",
                "defect_vos_kwargs",
                "ldauu",
                "ldaul",
                "xc",
                "no_wavecar",
                "potcar_set",
                "outcar",
                "contcar",
                "vicinage_radius",
                "cutoff",
                "displacement_distance",
                "volume_dir",
                "static_diele_dir",
                "ionic_diele_dir",
                "band_edge_dir",
                "dos_dir",
                "unitcell_json",
                "perfect_json",
                "chem_pot_yaml",
                "competing_phases_incar_setting"]


def get_user_settings(yaml_filename: str = "pydefect.yaml") -> dict:
    """Get the user specifying settings written in yaml_filename

    Note1: The yaml_filename is explored in the parent folders up to home
           or root directory until it's found. If it does not exist, empty
           dictionary is returned.
    Note2: When the key includes "/", the absolute path is added as a prefix.
           E.g., unitcell/unitcell.json -> /something/../unitcell/unitcell.json
    Note3: The value of "potcar_set: Mg_pv O_h" is "Mg_pv O_h" string, which
           is suited used for main default value.

    Args:
        yaml_filename (str): User setting yaml filename.

    Return:
        Dictionary of configs.
    """

    config_path = Path.cwd()
    home = Path.home()

    while True:
        if config_path == home or config_path == Path("/"):
            return {}

        f = config_path / yaml_filename
        if f.exists():
            with open(f, "r") as f:
                user_settings = yaml.load(f)
            break

        else:
            config_path = config_path.parent

    # Add full path
    for k, v in user_settings.items():
        if k not in setting_keys:
            raise ValueError(f"Keys in pydefect.yaml is invalid."
                             f"The candidate keys are {setting_keys}")
        if isinstance(v, str) and re.match(r'\S*/\S*', v):
            user_settings[k] = str(config_path / v)

    return user_settings


def get_default_args(function: Callable) -> dict:
    """Get the default values of the arguments in the method/function.

    Args:
        function (Callable):
            Method or function. when class is inserted, cls.__init__ is called.

    Return:
        default dict
    """
    defaults = {}
    signature_obj = signature(function)
    for name, param in signature_obj.parameters.items():
        if param.default != _empty:
            defaults[name] = param.default

    return defaults


def potcar_str2dict(potcar_list: Optional[str]) -> dict:
    """Sanitize the string type potcar setting to dict.

    An exmaple is "Mg_pv O_h" -> {"Mg": "Mg_pv", "O": "O_h"}
    If potcar_list is None, {} is returned.

    Args:
         potcar_list (str/None)

    Return:
         Dictionary of potcar types.
    """

    potcar_list = potcar_list.split() if potcar_list else []
    d = {}
    for p in potcar_list:
        element = p.split("_")[0]
        d[element] = p
    return d


def dict2list(d: dict) -> list:
    """Sanitize the string type potcar setting to dict.

    The string is also separated by space. An example is
    dict2list({"a": 1, "b": "2 3 4", "c": True}) =
                                 ["a", "1", "b", "2", "3", "4", "c", "True"]

    Args:
         d (dict)

    Return:
         list of flattened dict
    """

    d = d if d else {}
    flattened_list = []
    for k, v in d.items():
        flattened_list.append(k)
        if isinstance(v, str):
            flattened_list.extend(v.split())
        else:
            flattened_list.append(str(v))

    return flattened_list


def list2dict(flattened_list: Optional[list], key_candidates: list) -> dict:
    """Sanitize the list to dict with keys in the flags

    If a string in l does not exist in key_candidates, raise ValueError.

    key_candidates = ["ENCUT", "MAGMOM", "LWAVE", ...]
    list2dict(["ENCUT", "500", "MAGMOM", "4", "4", "LWAVE", "F"]) =
                        {"ENCUT": 500, "MAGMOM": [4, 4], "LWAVE": False}

    arg_list = ["ENCUT", "500", "MAGMAM", "4", "4"]
    raise ValueError

    Args:
        flattened_list (list): Input list
        key_candidates (list): List of key candidates, e.g., INCAR flags.
    Return:
        Sanitized dict
    """
    flattened_list = flattened_list if flattened_list else []

    d = {}
    key = None
    value_list = []

    def insert():
        if not value_list:
            raise ValueError(f"Invalid input: {flattened_list}.")
        if len(value_list) == 1:
            d[key] = value_list[0]
        else:
            d[key] = value_list

    for string in flattened_list:
        if key is None and string not in key_candidates:
            raise ValueError(f"Keys are invalid: {flattened_list}.")
        elif string in key_candidates:
            if key:
                insert()
                key = None
                value_list = []
            key = string
        else:
            if string.lower() == "true" or string == "T":
                value = True
            elif string.lower() == "false" or string == "F":
                value = False
            elif is_str_int(string):
                value = int(string)
            elif is_str_digit(string):
                value = float(string)
            else:
                value = string

            value_list.append(value)
    else:
        if key:
            insert()

    return d


def generate_objects_from_json_files(
        directory: str,
        json_filenames: List[str],
        classes: list,
        raise_error: bool = True) -> Optional[list]:
    """Generate class objects from directory and file names

    Note1: Filenames and classes correspond to each other. E.g.,
           filenames = ["defect_entry.json", "correction.json"]
           classes = [DefectEntry, ExtendedFnvCorrection]
    Note2: All the classes need to have load_json classmethod.

    Args:
        directory (str): Directory name where files locate
        json_filenames (list): List of file names 
        classes (list): List of classes
        raise_error (bool): Whether to raise error when parsing failed.

    Return:
        List of class objects.
    """

    objects = []
    for filename, class_obj in zip(json_filenames, classes):
        try:
            objects.append(class_obj.load_json(join(directory, filename)))
        except IOError:
            logger.warning(f"Parsing {filename} in {directory} failed.")
            if raise_error:
                raise
            else:
                logger.warning(f"Parsing {filename} in {directory} failed.")
                return

    return objects
