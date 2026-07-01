"""Submodule of the parser related to Gaussian files for ESTAMPES

A module providing the submodules to be used by ESTAMPES:

Sub-modules
-----------
fchk
    Formatted checkpoint files.
glog
    Gaussian output files.

See submodules for details.

The module also provides some conversion functions specific to
Gaussian, which apply to different types, like the structure of electric
quadrupole moment
"""
from collections.abc import Sequence

from estampes.base import InternalError


def g_elquad_LT_to_2D(vec: Sequence[str] | Sequence[int] | Sequence[float],
                      factor: float = 1.0
                      ) -> list[list[float]]:
    """Convert electric quadrupole from LT form to 2D tensor.

    Gaussian has a particular way to store internally the electric
    quadrupole, with the following sequence:
    XX, YY, ZZ, XY, XZ, YZ.
    The function builds a list of lists with the correct order.

    Parameters
    ----------
    vec
        Vector of data corresponding to the sequence stored in memory.
        Note that the list can have a longer size, only the first
        elements are used.
    factor
        Factor to be applied to the components, primarily intended to
        convert the elements in Gaussian dipole-quadrupole A tensor.
    """
    if len(vec) < 6:
        raise InternalError('Missing elements to build symmetric tensor')
    return [
        [factor*float(vec[0]), factor*float(vec[3]), factor*float(vec[4])],
        [factor*float(vec[3]), factor*float(vec[1]), factor*float(vec[5])],
        [factor*float(vec[4]), factor*float(vec[5]), factor*float(vec[2])]
    ]


def gfloat(number: str) -> float:
    """Convert a string in Fortran notation to float.

    Converts a string in number to float, taking care of D->E conversion.

    Parameters
    ----------
    number
        Number to convert.

    Returns
    -------
    float
        Converted number.
    """
    return float(number.replace('D', 'e'))
