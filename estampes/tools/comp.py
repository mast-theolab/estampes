"""Module providing basic computer-related functions

A basic module providing methods related to basic computer operations,
for ESTAMPES tools.

"""

import typing as tp

from estampes.base import ArgumentError


# =================
# Module Attributes
# =================

__BIT_POWERS = (' ', 'k', 'm', 'g', 't', 'p', 'e', 'z', 'y')


# ==============
# Module Methods
# ==============

def convert_storage(label: str) -> int:
    """Converts storage string to number of bytes.

    Given a storage specification with unit, converts it to a number of
    bytes.

    Parameters
    ----------
    label
        Storage specification (ex: 32GB).

    Returns
    -------
    int
        Number of bytes corresponding to the storage specification.
    """

    _label = label.strip().replace(' ', '').lower()
    if _label.endswith('ib'):
        metric = 1024
        offset = -3
        magnitude = _label[offset]
        str_num = _label[:offset]
    elif _label.endswith('b'):
        metric = 1000
        offset = -2
        magnitude = _label[offset]
        if magnitude not in __BIT_POWERS:
            magnitude = None
            offset = -1
        str_num = _label[:offset]
    elif _label.endswith('w'):  # Gaussian word unit (8 bytes / 64 bits)
        metric = 8192  # 8*1024
        offset = -2
        magnitude = _label[offset]
        if magnitude not in __BIT_POWERS:
            magnitude = None
            offset = -1
        str_num = _label[:offset]
    else:
        magnitude = None
        str_num = _label
    try:
        value = int(str_num)
    except ValueError as err:
        raise ValueError('Unsupported storage format.') from err
    if magnitude is not None:
        if magnitude not in __BIT_POWERS:
            raise ValueError('Unsupported byte unit.')
        unit = metric**(__BIT_POWERS.index(magnitude))
    else:
        unit = 1
    return value*unit


def bytes_units(num_bytes: int,
                prec: int = 0,
                binary: bool = False,
                power: tp.Optional[str] = None) -> str:
    """Converts a number of bytes to a human readable unit.

    Given a number of bytes (integer), returns a string with the
    storage expressed in a human readable form or in a specific
    unit.

    Parameters
    ----------
    num_bytes
        Total number of bytes
    prec
        Precision. Number of decimal digits.
    binary
        Uses binary metric instead of SI metric.
    power
        Select specific power of bytes (must be in `BIT_POWERS`).

    Returns
    -------
    str
        Storage in a specific or most adapted unit.
    """
    if binary:
        metric = 1024
        end_unit = 'iB'
    else:
        metric = 1000
        end_unit = 'B'
    if prec < 0:
        raise ValueError('Precision must be positive or null')
    fmt = '{{:.{:d}f}}{{}}'.format(prec)
    if power is None:
        magnitude = len(__BIT_POWERS)
        while magnitude > 0:
            value = num_bytes/metric**magnitude
            if value >= 1.:
                break
            magnitude -= 1
    else:
        try:
            magnitude = __BIT_POWERS.index(power.lower()[0])
            value = num_bytes/metric**magnitude
        except ValueError as err:
            raise ValueError('Unknown byte power.') from err
    return fmt.format(value, __BIT_POWERS[magnitude].upper()+end_unit)


def lin_storage(shape: tp.Sequence[int],
                *indx: int,
                dfmt='C') -> int:
    """Converts an list of indexes to linear storage.

    Converts a list of indexes to a linear storage based on shape.
    The shape can contain:
    * positive or negative values, treated the same way.
    * 0, for a dim. symm. with previous one, stored as lower triangular.

    Parameters
    ----------
    shape
        Shape of the array.
    indx
        Indexes along each coordinate.
    dfmt
        Storage format ("C"-like or "F"ortran-like).

    Returns
    -------
    int
        Linear storage.

    Notes
    -----
    * For C storage format, the indexes should start at 0; 1 for Fortran.
    * Fortran index is shifted to match Python indexing start at 0.
    """
    if len(indx) > len(shape):
        raise ArgumentError('Too many indexes')
    num = len(indx)
    res = 0
    offset = 1
    lprv = 1
    if dfmt.upper() == 'C':
        for i in range(num, 0, -1):
            res += indx[i-1]*offset
            lcur = abs(shape[i-1])
            if lcur > 0:
                offset *= lcur
                lprv = lcur
            else:
                offset *= lprv
    else:
        for i in range(num):
            res += indx[i]
            lcur = abs(shape[i-1])
            if lcur > 0:
                offset *= lcur
                lprv = lcur
            else:
                offset *= lprv
        offset -= 1
    return res
