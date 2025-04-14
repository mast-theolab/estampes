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
        List of atoms labels (as string or integer).
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


def eckart_orient(at_crd: TypeAtCrd,
                  at_mass: TypeAtMas) -> np.ndarray:
    """Return Eckart orientation.

    Given a set of atomic coordinates and the masses, return a new
    orientation suitable to meet the Eckart conditions.

    Parameters
    ----------
    at_crd
        Atomic coordinates as XYZ vectors, as a 2D array.
    at_mass
        Atomic masses, as a 1D array.

    Returns
    -------
    np.ndarray
        (n_atoms,3) New oriented structure.
    """
    c_new = at_crd - center_of_mass(at_crd, at_mass)
    _, pmom_evec = inertia_mom(at_crd, at_mass)
    c_new = c_new @ pmom_evec

    return c_new


def inertia_mom(at_crd: TypeAtCrd,
                at_mass: TypeAtMas) -> tp.Tuple[np.ndarray, np.ndarray]:
    """Build inertia moments

    Computes the inertia moments and the returns the principal moments
    of inertia and the corresponding eigenvectors.

    Parameters
    ----------
    at_crd
        Atomic coordinates as XYZ vectors, as a 2D array.
    at_mass
        Atomic masses, as a 1D array.

    Returns
    -------
    np.ndarray
        (3) Principal moments of inertia.
    np.ndarray
        (3,3) Eigenvectors from the diagonalization.
    """
    if at_crd.shape[0] != at_mass.shape[0]:
        raise IndexError('Size inconsistency between masses and coordinates.')

    Imom = np.zeros((3, 3))
    for i, mi in enumerate(at_mass):
        Imom[0, 0] += mi*(at_crd[i, 1]**2 + at_crd[i, 2]**2)
        Imom[1, 1] += mi*(at_crd[i, 0]**2 + at_crd[i, 2]**2)
        Imom[2, 2] += mi*(at_crd[i, 0]**2 + at_crd[i, 1]**2)
        Imom[0, 1] -= mi * at_crd[i, 0] * at_crd[i, 1]
        Imom[0, 2] -= mi * at_crd[i, 0] * at_crd[i, 2]
        Imom[1, 2] -= mi * at_crd[i, 1] * at_crd[i, 2]
    Imom[1, 0] = Imom[0, 1]
    Imom[2, 0] = Imom[0, 2]
    Imom[2, 1] = Imom[1, 2]

    if abs(Imom[0, 1]) + abs(Imom[0, 2]) + abs(Imom[1, 2]) > 1.0e-6:
        eigval, eigvec = np.linalg.eigh(Imom, UPLO='u')
    else:
        eigval = np.diag(Imom)
        eigvec = np.eye(3)

    return eigval, eigvec
