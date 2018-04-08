#!/usr/bin/zsh

python3 setup.py build_ext -i && python3 ase_generation_supercell.py POSCAR_yk
