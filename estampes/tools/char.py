"""Module providing basic character strings-related functions

A basic module providing methods related to character string operations,
for ESTAMPES tools.

"""

import re
import typing as tp

from estampes.base.errors import ArgumentError


# ==============
# Module Methods
# ==============

def convert_expr(expr: str,
                 variable: tp.Optional[str] = None,
                 natural: bool = True) -> str:
    """Convert mathematical expression to correct Python form.

    A very basic function to convert a mathematical expression given in
    a more natural language to Python-compatible expressions.
    It is intended primarily for expressions like scaling or "scaled"
    units, like `"10^-40 esu^2 cm^2"`, not to provide a complete
    mathematical framework.

    Parameters
    ----------
    expr
        String expression to interpret.
    variable
        Variable accepted in the expression.
    natural
        Assumes more natural function names (e.g., ln, log).

    Returns
    -------
    string
        Python-compatible expression.

    Raises
    ------
    ValueError
        String contains invalid elements
    """
    # Very simple test on possible safety issues.
    if 'exec' in expr or 'lambda' in expr:
        raise ValueError('Pure mathematical expressions expected')

    if natural:
        _expr = expr.strip().replace('log', 'log10').replace('ln', 'log')
    else:
        _expr = expr

    if variable is not None:
        if not re.search(r'\b'+variable+r'\b', _expr):
            char = _expr[-1]
            if char in '+-*/^':
                _expr += variable
            else:
                _expr += '*' + variable
        pattern = r'\^([-\+]\d*' + variable + r'?)'
    else:
        pattern = r'\^([-\+]\d*)'

    if '^' in _expr:
        # Replace case without protecting parentheses (ex: 10^-3)
        _expr = re.sub(pattern, r'**(\1)', _expr)
        # Replace other cases
        _expr = _expr.replace('^', '**')

    return _expr


def unit_to_tex(unit: str,
                only_dot: bool = False,
                use_cdot: bool = True) -> str:
    """Convert a unit to TeX format.

    Converts a unit written in compact format to a proper LaTeX format.

    Parameters
    ----------
    unit
        Unit label.
    only_dot
        If True, slashes are converted to dot and the exponent corrected.
        Example: /cm -> .cm^{-1}, /cm2 -> .cm^{-2}
    use_cdot
        If True, cdot is used as "multiplier", otherwise the simple dot.

    Notes
    -----
    The function is not very robust and assume the string is correctly
    built following the conventions used inside ESTAMPES.
    """
    if use_cdot:
        dot = r'$\cdot$'
    else:
        dot = '.'

    new_unit = ''
    # Check first if there is a factor.
    # The factor is assumed separated from the rest by one or more blanks.
    items = unit.split(' ', maxsplit=1)
    if len(items) == 2:
        num = re.fullmatch(
            r'(?P<int>\d+)\b\^?(?P<exp>(?:[-+]?\d+|\([-+]?\d+\)))?',
            items[0])
        if num:
            res = num.groupdict()
            if res['exp'] is None:
                new_unit += f'{res["int"]}'
            else:
                new_unit += f'{res["int"]}^{{{res["exp"]}}}'
            unit_ = items[1]
        else:
            unit_ = unit
    else:
        unit_ = unit

    items = re.split(r'(?=[./])', unit_)
    # Since unit may start with '/', the first block may be empty, bypass
    if unit_.startswith('/') and not items[0]:
        items = items[1:]
    for item in items:
        res = re.fullmatch(
            r'(?P<op>[./]?)(?P<unit>[^0-9+-]+)(?P<sign>[-+]?)(?P<exp>\d*)',
            item.replace(' ', ''))
        if not res:
            raise ArgumentError(f'Wrong format of unit, item: {item}')
        oper, unit, sign, exp = res.groups()

        if sign and not exp:
            raise ArgumentError(f'Wrong format of unit, item: {item}')

        if oper == '/':
            if only_dot:
                if not exp:
                    sign = '-'
                    exp = '1'
                else:
                    if not sign:
                        sign = '-'
                    elif sign == '+':
                        sign = '-'
                oper = dot
        elif oper == '.':
            oper = dot
        # We remove multiplication symbol if nothing before.
        if not new_unit and oper == dot:
            oper = ''
        if sign+exp:
            new_unit += f'{oper}{unit}$^{{{sign+exp}}}$'
        else:
            new_unit += f'{oper}{unit}'

    return new_unit
