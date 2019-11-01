""" Module containing functions used in main.py
"""

from pydefect.capture.capture_processes import CaptureProcesses


def identify_capture_processes(args):
    """ Identify the carrier capture processes which are possible from a given
    set of defect data"""

    try:
        capture_processes = CaptureProcesses.from_directory(args.defect_dir)
        capture_processes.to_yaml_file()
        capture_processes.print_to_terminal()

    except IOError:
        raise FileNotFoundError(args.defect_dir, 'does not contain multiple '
                                                 'defects.json')
