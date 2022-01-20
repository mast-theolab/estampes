"""Toolbox for molecules

Module providing tools to manipulate data relative to molecules.

Methods
-------
list_bonds
    Finds and lists bonds between atoms.
"""

import numpy as np
import numpy.typing as npt

from estampes.base import TypeAtLab, TypeBonds
from estampes.data.atom import atomic_data


# ================
# Module Constants
# ================

TypeAtCrd = npt.ArrayLike


# ==============
# Module Methods
# ==============

def list_bonds(at_lab: TypeAtLab,
               at_crd: TypeAtCrd,
               rtol: float = 1.1) -> TypeBonds:
    """Finds and lists bonds between atoms.

    Parameters
    ----------
    at_lab
        List of atoms labels (as string).
    at_crd
        Atomic coordinates as XYZ vectors, in Ang.
    rtol
        Radius tolerance, i.e. scaling factor applied to Rcov for bond
          identification.

    Returns
    -------
    list
        List of bonds as `(atom1, atom2)`.
    """
    bonds = []
    natoms = len(at_lab)
    atdat = atomic_data(*at_lab)
    for i in range(natoms-1):
        xyz_i = at_crd[i]
        rad_i = max(atdat[at_lab[i]]['rcov'], key=lambda x: x or 0.)
        if rad_i is None:
            raise ValueError('Missing rcov for atom {}'.format(i))
        for j in range(i+1, natoms):
            xyz_j = at_crd[j]
            rad_j = max(atdat[at_lab[j]]['rcov'], key=lambda x: x or 0.)
            if rad_j is None:
                raise ValueError('Missing rcov for atom {}'.format(j))
            if np.linalg.norm(xyz_j-xyz_i) < (rad_i + rad_j)*rtol/100.:
                bonds.append((i, j))

    return bonds
