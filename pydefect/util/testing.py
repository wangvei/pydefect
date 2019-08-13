# -*- coding: utf-8 -*-
from typing import Union
from pathlib import Path
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
import tempfile

class PydefectTest(PymatgenTest):
    """
    Extends PymatgenTest with some path modification.
    """
    MODULE_DIR = Path(__file__).absolute().parent
    TEST_FILES_DIR = MODULE_DIR / ".." / ".." / "test_files"
    POSCARS_DIR = TEST_FILES_DIR / "poscars"
    CORE_DIR = TEST_FILES_DIR

    @classmethod
    def get_structure_by_name(cls, name):
        filename = cls.POSCARS_DIR / ("POSCAR-" + name)
        return Structure.from_file(filename)

    @classmethod
    def get_object_by_name(cls, classmethod, name):
        filename = cls.TEST_FILES_DIR / name
        return classmethod(filename)

    @classmethod
    def get_filename(cls, name):
        return cls.TEST_FILES_DIR / name


