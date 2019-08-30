from setuptools import setup, find_packages
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = {}
ext_modules = []

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pydefect-yuuukuma',
    version='0.1.dev0',
    author='Yu Kumagai',
    author_email='yuuukuma@gmail.com',
    url='https://github.com/oba-group/pydefect"',
    packages=find_packages(),
    license='MIT license',
    description="Integrated enveironment for first-principles point-defect "
                "calculations using vasp",
    long_description=long_description,
    classifiers=[
        'Programming Language :: Python :: 3.6',
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=['numpy', 'pymatgen', 'monty', 'matplotlib', 'argcomplete',
                      'seekpath', 'spglib', 'scipy', 'ase', 'tqdm', 'yaml'],
    cmdclass = cmdclass,
    ext_modules=ext_modules,
)
