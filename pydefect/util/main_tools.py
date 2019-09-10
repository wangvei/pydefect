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
                "vasprun",
                "procar",
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
                "chem_pot_yaml"]


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
            raise ValueError(f"Key in pydefect.yaml {k} is invalid."
                             f"The candidate keys are {setting_keys}")
        if isinstance(v, str) and re.match(r'\S*/\S*', v):
            user_settings[k] = str(config_path / v)

    return user_settings


def get_default_args(function: Callable) -> dict:
    """Get the default values of the arguments in the method/function.

    inspect._empty means no default.

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
