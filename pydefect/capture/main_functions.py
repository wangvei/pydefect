""" Module containing functions used in main.py
"""

from pydefect.capture.capture_processes import (locate_defect_jsons,
                                                identify_processes)


def identify_capture_processes(args):
    """ Identify the carrier capture processes which are possible from a given
    set of defect data"""

    try:
        defect_json_files = locate_defect_jsons(args.defect_dir)
        identify_processes(defect_json_files)

    except IOError:
        raise FileNotFoundError(args.defect_dir, 'does not contain multiple '
                                                 'defects.json')
