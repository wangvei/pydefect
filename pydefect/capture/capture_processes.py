""" Module to identify the possible carrier capture processes which can be
considered from a given set of defect data
"""

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

    print("I ")
    return