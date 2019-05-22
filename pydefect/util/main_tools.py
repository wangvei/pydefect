# -*- coding: utf-8 -*-
from inspect import signature, _empty

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def overwrite_default_args(class_method, main_args):
    """ Use the defaults in class.classmethod

    Args:
        class_method (classmethod):
            classmethod. When using __init__, class is fine.
        main_args (dict):
            Args set by main

    Return:
        args (dict): Overwritten args by options
    """

    args_with_default = []
    sig = signature(class_method)

    for name, param in sig.parameters.items():
        if param.default != _empty:
            args_with_default.append(name)

    args = {}
    for a in args_with_default:
        if hasattr(main_args, a):
            if getattr(main_args, a) is not None:
                args[a] = getattr(main_args, a)

    return args


def get_default_args(class_method):
    defaults = {}
    sig = signature(class_method)
    for name, param in sig.parameters.items():
        if param.default != _empty:
            defaults[name] = param.default

    return defaults
