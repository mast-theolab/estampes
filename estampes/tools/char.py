"""Module providing basic character strings-related functions

A basic module providing methods related to character string operations,
  for ESTAMPES tools.

Methods
-------
convert_expr
    Converts mathematical expression to correct Python form.
"""

import re
import typing as tp


# ==============
# Module Methods
# ==============

def convert_expr(expr: str,
                 variable: tp.Optional[str] = None,
                 natural: bool = True) -> str:
    """Converts mathematical expression to correct Python form.

    A very basic function to convert a mathematical expression given in
      a more natural language to Python-compatible expressions.
    It is intended primarily for expressions like scaling or "scaled"
      units, like "10^-40 esu^2 cm^2", not to provide a complete
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
