"""Toolbox for molecular vibrations

Module providing tools to manipulate data relative to vibrations.

Methods
-------
orient_modes
    Sets a unique orientation for normal coordinates.
"""

import numpy as np

# ==============
# Module Methods
# ==============


def orient_modes(Lmat: np.ndarray) -> np.ndarray:
    """Sets a unique orientation for normal coordinates.

    Sets the sign for each column of `Lmat` so that the highest element
      is positive.  This way, the vectors can adopt a consistent
      orientation, especially for the numerical differentiation, where
      it is critical.
    A small shift is given to this highest value to avoid uncertainty
      if two elements have very close or equal magnitudes but
      potentially opposite signs.

    Parameters
    ----------
    Lmat
        Transformation matrix from cart to normal coord. (nvib,3natoms).

    Returns
    -------
    :obj:`np.ndarray`
        `Lmat` with corrected columns.
    """
    if not isinstance(Lmat, np.ndarray):
        raise IndexError('Wrong parameter, Lmat should be an array')
    if Lmat.ndim != 2:
        raise IndexError('Lmat is expected to have dimension 2')
    shift = 1.0e-6
    nvib, nat3 = Lmat.shape
    for i in range(nvib):
        xmax = abs(Lmat[i, 0]) + shift
        imax = 0
        for j in range(1, nat3):
            if abs(Lmat[i, j]) > xmax:
                xmax = abs(Lmat[i, j]) + shift
                imax = j
        if Lmat[i, imax] < 0.:
            Lmat[i, :] *= -1.

    return Lmat
