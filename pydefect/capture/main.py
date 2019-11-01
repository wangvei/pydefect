#!/usr/bin/env python3

""" Command line interface for the capture sub-package within pydefect.
"""

import argparse
from typing import Union

from pydefect.capture.main_functions import identify_capture_processes
from vise.util.main_tools import get_user_settings, dict2list


def main():

    # copied from pydefect/main - to refactor
    def simple_override(d: dict, keys: Union[list, str]) -> None:
        """Override dict if keys exist in user_settings.

        When the value in the user_settings is a dict, it will be changed to
        list using dict2list.
        """
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            if key in user_settings:
                v = user_settings[key]
                if isinstance(v, dict):
                    v = dict2list(v)
                d[key] = v

    setting_keys = ["defect_dir"]

    user_settings = get_user_settings(yaml_filename="capture_settings.yaml",
                                      setting_keys=setting_keys)

    parser = argparse.ArgumentParser(
        description="""
        capture is a sub-package that helps pydefect users to do 
        carrier capture calculations with the VASP code.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    parser_identify_processes = subparsers.add_parser(
        name="identify_processes",
        description="Tools for identifying possible carrier capture "
                    "processes from defect calculation data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ip'])

    gs_defaults = {"defect_dir": "../defects/"}

    simple_override(gs_defaults, list(gs_defaults.keys()))  # redundant arg?

    parser_identify_processes.add_argument(
        "--defect_dir", dest="defect_dir",
        default=gs_defaults["defect_dir"], type=str,
        help="Identify capture processes from a directory containing "
             "defects.json files (recursive search)")

    del gs_defaults  # is this for memory considerations only?

    parser_identify_processes.set_defaults(func=identify_capture_processes)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
