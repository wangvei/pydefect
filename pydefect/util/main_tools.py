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


def list2dict(arg_list, flags):
    """
    flags: incar flags

    arg_list = ["ENCUT", "500", "MAGMOM", "4", "4", "LWAVE": "F"]
    return {"ENCUT": 500, "MAGMOM": [4, 4], "LWAVE": False}

    arg_list = ["ENCUT", "500", "MAGMAM", "4", "4"]
    raise ValueError

    :param arg_list:
    :param flags:
    :return:
    """
    flag_indices = []
    arg_list = [] if arg_list is None else arg_list

    for i, element in enumerate(arg_list):
        if element in flags:
            flag_indices.append(i)

    num_flag = len(flag_indices)
    flag_indices.append(len(arg_list))

    d = {}
    for i in range(num_flag):
        if flag_indices[i + 1] - flag_indices[i] == 2:
            value_str = arg_list[flag_indices[i] + 1]
            if value_str.lower() == "true" or value_str == "T":
                value = True
            elif value_str.lower() == "false" or value_str == "F":
                value = False
            elif value_str.isdigit():
                value = int(value_str)
            else:
                try:
                    value = float(value_str)
                except ValueError:
                    value = value_str

        else:
            try:
                value = [float(arg_list[v]) for v in range(flag_indices[i] + 1,
                                                           flag_indices[i + 1])]
            except ValueError:
                raise ValueError("Invalid input {}".format(arg_list))

        d[arg_list[flag_indices[i]]] = value

    return d
