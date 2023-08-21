"""Toolbox for molecules

Module providing tools to manipulate data relative to molecules.

"""

import typing as tp

import numpy as np

from estampes.base import TypeAtCrd, TypeAtLab, TypeAtMas, TypeBonds
from estampes.data.atom import atomic_data


# ================
# Module Constants
# ================

# Typeat_crd = npt.ArrayLike
# Typeat_mass = npt.ArrayLike


# ==============
# Module Methods
# ==============

def center_of_mass(at_crd: TypeAtCrd,
                   at_mass: TypeAtMas) -> np.ndarray:
    """Computes the center of a mass.

    Computes and returns the center of mass.

    Parameters
    ----------
    at_crd
        Atomic coordinates as XYZ vectors, as a 2D numpy array.
    at_mass
        Atomic masses, as a 1D numpy array.

    Returns
    -------
    np.ndarray
        (3) Coordinates of the center of mass.

    Raises
    ------
    IndexError
        Size inconsistency between masses and coordinates.
    """
    if at_crd.shape[0] != at_mass.shape[0]:
        raise IndexError('Size inconsistency between masses and coordinates.')
    return np.einsum('ij,i->j', at_crd, at_mass)/np.sum(at_mass)


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
