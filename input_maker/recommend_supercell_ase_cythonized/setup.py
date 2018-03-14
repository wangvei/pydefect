from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

setup(name="ase_cython",
      ext_modules=cythonize("ase_cython.pyx"),
      include_dirs = [numpy.get_include()])
