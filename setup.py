from setuptools import setup, find_packages
from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = { }
ext_modules = [ ]

if use_cython:
    ext_modules += [
    Extension("pydefect.", [ "pydefect/analysis/recommend_supercell_ase_cythonized/ase_cython.pyx"]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
    Extension("pydefect.", [ "pydefect/analysis/recommend_supercell_ase_cythonized/ase_cython.c"]),
    ]

setup(
    name='PyDefect',
    version='0.1.dev0',
    author='Yu Kumagai',
    author_email='yuuukuma@gmail.com',
    url='https://scholar.google.co.jp/citations?user=xST4MSEAAAAJ&hl=ja',
    packages=['bin','test',],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ],
    install_requires=['numpy', 'pymatgen', 'monty', 'matplotlib'],
    cmdclass = cmdclass,
    ext_modules=ext_modules,
)
