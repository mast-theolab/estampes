"""Provide tools to facilitate some output operations.

The module provides some simple tools for common output operations.
"""
import re
from collections.abc import Sequence

from estampes.base import ArgumentError


def fortran_fmt_D(number_s: float | Sequence[float], fmt: str) -> str:
    """Print simulating the double precision format.

    Simulate the printing a value `num` with a double precision format
    specification.
    Fortran double precision format has the form: Dw.d[Ee] with:

    - w: total width of the field.
    - d: number of decimal digits.
    - e: number of digits for the exponentiation.

    with the integer part always 0.

    If multiple values are provided, the resulting string is a
    concatenation of each one.

    Parameters
    ----------
    number_s
        Number to represent.
    fmt
        Fortran format.

    Returns
    -------
    str
        Representation of the content `number_s` with the proper format.
    """
    def num_to_repr(num: float) -> str:
        n_sign = 1 if num < 0 else 0
        if comps['e'] is not None:
            # sign + 0 + dot + decimal_part + E + sign_E + length_e
            w_min = n_sign + 2 + len_d + 2 + len_e
        else:
            # sign + 0 + dot + decimal_part + 4 (E+dd or +ddd)
            w_min = n_sign + 2 + len_d + 4
        if len_w < w_min:
            return len_w*'*'

        # Python represents number with integer part > 0 so we must correct d
        # by reducing it by -1
        py_fmt = f'{{:.{len_d-1}e}}'
        py_res = py_fmt.format(abs(num))
        p_num, p_exp = py_res.split('e')
        if num > 0:
            p_num = '0.' + p_num.replace('.', '')
        else:
            p_num = '-0.' + p_num.replace('.', '')
        if comps['e'] is not None:
            if len_e + 1 < len(p_exp):
                return len_w*'*'
            else:
                py_fmt = f'{{:+{len_e}d}}'
                p_exp = 'D' + py_fmt.format(int(p_exp))
        else:
            if len(p_exp) == 3:
                p_exp = 'D' + p_exp
            elif len(p_exp) == 4:
                p_exp = py_fmt.format(int(p_exp))
            else:
                return len_w*'*'

        py_fmt = f'{{:>{len_w}s}}'
        return py_fmt.format(p_num+p_exp)

    parse = re.match(r'D(?P<w>\d+)\.(?P<d>\d+)(?:E(?P<e>\d+))?', fmt, re.I)
    if parse is None:
        raise ArgumentError('fmt', 'Incorrect Fortran format')

    comps = parse.groupdict()
    len_w = int(comps['w'])
    len_d = int(comps['d'])
    len_e = int(comps['e'] if comps['e'] is not None else 0)
    if 0 < len_e < 2:
        raise NotImplementedError('Small specs for D not implemented.')

    if isinstance(number_s, float):
        return num_to_repr(number_s)
    elif isinstance(number_s, (tuple, list)):
        repr_s = []
        for item in number_s:
            repr_s.append(num_to_repr(item))
        return ''.join(repr_s)
    else:
        raise TypeError('Unsupported representation of number_s')
