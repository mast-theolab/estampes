"""Toolbox for molecules

Module providing tools to manipulate data relative to molecules.

"""

import numpy as np
import numpy.typing as npt

from estampes.base.types_npt import AtsMasType
from estampes.base import AtsCrdType, AtsLabType, BondsType
from estampes.data.atom import atomic_data


# ==============
# Module Methods
# ==============

def center_of_mass(at_crd: AtsCrdType,
                   at_mass: AtsMasType) -> np.ndarray:
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
    at_crd = np.asarray(at_crd)
    at_mass = np.asarray(at_mass)
    if at_crd.shape[0] != at_mass.shape[0]:
        raise IndexError('Size inconsistency between masses and coordinates.')
    return np.einsum('ij,i->j', at_crd, at_mass)/np.sum(at_mass)


def list_bonds(at_lab: AtsLabType,
               at_crd: AtsCrdType,
               rtol: float = 1.1) -> BondsType:
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
        xyz_i = np.asarray(at_crd[i])
        rad_i = max(atdat[at_lab[i]]['rcov'], key=lambda x: x or 0.)
        if rad_i is None:
            raise ValueError(f'Missing rcov for atom {i}')
        for j in range(i+1, natoms):
            xyz_j = np.asarray(at_crd[j])
            rad_j = max(atdat[at_lab[j]]['rcov'], key=lambda x: x or 0.)
            if rad_j is None:
                raise ValueError(f'Missing rcov for atom {j}')
            if np.linalg.norm(xyz_j-xyz_i) < (rad_i + rad_j)*rtol/100.:
                bonds.append((i, j))

    return bonds


def eckart_orient(at_crd: AtsCrdType,
                  at_mass: AtsMasType,
                  get_rotation: bool = False,
                  get_translation: bool = False
                  ) -> dict[str, npt.NDArray]:
    """Return Eckart orientation.

    Given a set of atomic coordinates and the masses, return a new
    orientation suitable to meet the Eckart conditions, and optionally
    additional parameters.

    Parameters
    ----------
    at_crd
        Atomic coordinates as XYZ vectors, as a 2D array.
    at_mass
        Atomic masses, as a 1D array.
    get_rotation
        Return the rotation matrix.
    get_translation
        Return the translation vector to origin as well.

    Returns
    -------
    dict
        Dictionary containing the following NumPy arrays, depending
        on options:

        'coords': (n_atoms,3) New oriented structure.
        'rotmat': (3,3) rotation matrix (not by default).
        'transl': (3) translation vector (not by default).
    """
    tvec = - center_of_mass(at_crd, at_mass)
    c_new = np.asarray(at_crd) + tvec
    _, pmom_evec = inertia_mom(at_crd, at_mass)
    # Check if chirality preserved:
    x = (
        pmom_evec[0, 0] * (
            pmom_evec[1, 1]*pmom_evec[2, 2]
            - pmom_evec[2, 1]*pmom_evec[1, 2])
        + pmom_evec[0, 1] * (
            pmom_evec[1, 2]*pmom_evec[2, 0]
            - pmom_evec[1, 0]*pmom_evec[2, 2])
        + pmom_evec[0, 2] * (
            pmom_evec[1, 0]*pmom_evec[2, 1]
            - pmom_evec[1, 1]*pmom_evec[2, 0]))
    if x < 0.0:
        pmom_evec[:, 0] *= -1

    c_new = c_new @ pmom_evec

    res = {'coords': c_new}
    if get_rotation:
        res['rotmat'] = pmom_evec
    if get_translation:
        res['transl'] = tvec

    return res


def inertia_mom(at_crd: AtsCrdType,
                at_mass: AtsMasType) -> tuple[npt.NDArray, npt.NDArray]:
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
    at_crd = np.asarray(at_crd)
    at_mass = np.asarray(at_mass)
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
