"""Toolbox related to anharmonic treatment.

Module providing tools to parse/analyse anharmonic data.
"""

import typing as tp

from estampes.base import QLabel
from estampes.base.errors import ArgumentError, ParsingError
from estampes.parser import DataFile


def variational_notation(dfile: DataFile | None = None,
                         assign: list[tp.Any] | None = None) -> bool:
    """Check if variational notation used for vibrational states.

    Returns True if the variational notation has been used to represent
    the vibrational states.
    The function accepts in input a Datafile or the assignment
    information, the latter taking precedence.

    Parameters
    ----------
    dfile
        Datafile.
    assign
        Assignment information, as provided by the internal parser.

    Returns
    -------
    bool
        True if the variational notation is used.

    Raises
    ------
    ValueError
        The notation cannot be recognized or is missing.
    """
    if assign is not None:
        try:
            answer = assign[1][1][0][1] == 0
        except IndexError as err:
            raise ValueError('Incorrect assignment information') from err
    elif dfile is not None:
        res = dfile.get_data(key=QLabel(quantity='vtrans', level='A'))
        answer = res['key'].data[1][1][0][1] == 0
    else:
        raise ArgumentError('Missing input data')

    return answer
