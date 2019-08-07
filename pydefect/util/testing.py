# -*- coding: utf-8 -*-
from pathlib import Path
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure


class PydefectTest(PymatgenTest):
    """
    Extends PymatgenTest with some path modification.
    """
    MODULE_DIR = Path(__file__).absolute().parent
    TEST_FILES_DIR = MODULE_DIR / ".." / ".." / "test_files"
    POSCARS_DIR = TEST_FILES_DIR / "poscars"

    @classmethod
    def get_structure_by_name(cls, name):
        filename = cls.POSCARS_DIR / ("POSCAR-" + name)
        return Structure.from_file(filename)