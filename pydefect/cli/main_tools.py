# -*- coding: utf-8 -*-

from os.path import join
from typing import Optional, List

from pydefect.util.logger import get_logger

logger = get_logger(__name__)


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
                objects.append(None)

    return objects
