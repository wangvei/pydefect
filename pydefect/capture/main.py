#!/usr/bin/env python3

""" Command line interface for the capture sub-package within pydefect.
"""

import argparse

from pydefect.capture.main_functions import identify_capture_processes
from vise.util.main_tools import get_user_settings
from pydefect.util.main_tools import simple_override


def main():
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

    gs_defaults = simple_override(gs_defaults,
                                  user_settings,
                                  list(gs_defaults.keys())) # redundant argument?

    parser_identify_processes.add_argument(
        "--defect_dir", dest="defect_dir",
        default=gs_defaults["defect_dir"], type=str,
        help="Identify capture processes from a directory containing "
             "defects.json files (recursive search)")

    del gs_defaults # is this for memory conisderations only?

    parser_identify_processes.set_defaults(func=identify_capture_processes)

if __name__ == "__main__":
    main()