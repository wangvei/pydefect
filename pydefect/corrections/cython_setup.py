# -*- coding: utf-8 -*-
from Cython.Build import cythonize
from distutils.core import setup

setup(ext_modules=cythonize("pydefect.corrections.calc_ewald_sum.pyx"))

