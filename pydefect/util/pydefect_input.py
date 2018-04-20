#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
A master convenience script with many tools for vasp and structure analysis.
"""

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "April 19, 2018"


def input_pydefect(args):
    if args.potcar_dirs:
        setup_potcars(args)
    elif args.install:
        install_software(args)
    elif args.var_spec:
        add_config_var(args)