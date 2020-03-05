# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Callable, Union

from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest


class PydefectTest(PymatgenTest):
    """ Extends PymatgenTest with some path modification. """
    MODULE_DIR = Path(__file__).absolute().parent
    TEST_FILES_DIR = MODULE_DIR / ".." / "test_files"
    DEFECTS_MGO_DIR = TEST_FILES_DIR / "defects" / "MgO"
    POSCARS_DIR = TEST_FILES_DIR / "poscars"
    CORE_DIR = TEST_FILES_DIR

    @classmethod
    def get_structure_by_name(cls, name: str):
        filename = cls.POSCARS_DIR / ("POSCAR-" + name)
        return Structure.from_file(filename)

    @classmethod
    def get_object_by_name(cls, method: Callable, names: Union[list, str]):
        """ return cls by passing classmethod returning cls """
        if isinstance(names, str):
            names = [names]

        filename = cls.TEST_FILES_DIR
        for name in names:
            filename = filename / name

        return method(filename)

    @classmethod
    def get_filename(cls, name):
        return cls.TEST_FILES_DIR / name


