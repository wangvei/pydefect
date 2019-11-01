""" Module to identify the possible carrier capture processes which can be
considered from a given set of defect data
"""

import json
import itertools

from typing import Generator
from pathlib import Path


def locate_defect_jsons(defect_dir: str) -> Generator[Path, None, None]:
    """Returns a list of all the defect.json filepaths contained in a given
    directory (recursive search).

    Parameters:
        defect_dir (str): Search directory.

    Returns:
        Generator[Path, None, None]: Paths to the defect.json files contained
        in defect_dir.
    """

    return Path(defect_dir).rglob('defect.json')


def identify_processes(path_list: Generator[Path, None, None]) -> None:
    """Identifies possible carrier capture processes which can be considered
    from a given set of defect.json files.

    Parameters:
        path_list (Generator[Path, None, None]): Paths to defect.json files.

    Returns:
        None.
    """

    valid_defects = {}
    for path in path_list:
        with path.open() as file:
            data = json.load(file)
            if data['@class'] == "Defect" and data['is_converged'] == True:
                valid_defects[data['name'] + "_" + str(data['charge'])] = {
                    'name': data['name'], 'charge': data['charge'],
                    'path': path.resolve(),
                    'band_edge_states': data['band_edge_states']}

    capture_processes = {}
    for d1, d2 in itertools.combinations(valid_defects,2):
        if d1['name']==d2['name'] and abs(d1['charge']-d2['charge'])==1:
            capture_processes[d1]

    # loop through and record key (name, charge, localised, path) info if
    # 1) defect class 3) converged into a dict
    # in new dict, create dictionary entries where charge state diff by one with path/name/initial->final
    # option to calc delta Q and delta E (use displacements)
    # option to do unconverged
    # option to include non-local states
    # print out to screen / to file
    # error handling and tests
