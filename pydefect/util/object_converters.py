from pymatgen import Spin


def spin_key_to_str(arg: dict, value_to_str=False):
    if arg is not None:
        if value_to_str:
            return {str(spin): str(v) for spin, v in arg.items()}
        else:
            return {str(spin): v for spin, v in arg.items()}
    else:
        return


def str_key_to_spin(arg: dict, method_from_str_for_value=None):
    if arg is not None:
        x = {}
        for spin, value in arg.items():
            x[Spin(int(spin))] = value
            if method_from_str_for_value:
                x[Spin(int(spin))] = method_from_str_for_value(value)
            else:
                x[Spin(int(spin))] = value
        return x
    else:
        return