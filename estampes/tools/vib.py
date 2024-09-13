"""Toolbox for molecular vibrations

Module providing tools to manipulate data relative to vibrations.

"""

import typing as tp

import numpy as np

from estampes.base import TypeAtCrd, TypeAtMas, ArgumentError


# ==============
# Module Methods
# ==============

def build_dusch_J(Lmat_A: np.ndarray,
                  Lmat_B: np.ndarray,
                  orth_L: bool = True,
                  at_mass_A: tp.Optional[np.ndarray] = None,
                  at_mass_B: tp.Optional[np.ndarray] = None) -> np.ndarray:
    """Build Duschinsky matrix J.

    Build the Duschinsky matrix J corresponding to:

    :math:`Q_A = J Q_B + K`

    with :math:`J = L_A M_A^{1/2} M_B^{-1/2} L_B^{-1}`.

    By default, the L matrices are assumed pseudo-orthogonal and the
    masses are the same, so :math:`J = L_A L_B^{T}`

    Parameters
    ----------
    Lmat_A
        Eigenvectors of the Hessian matrix of system A, as (N,3Na) array.
    Lmat_B
        Eigenvectors of the Hessian matrix of system B, as (N,3Na) array.
    orth_L
        If True, `Lmat_A` and `Lmat_B` are assumed pseudo-orthogonal.
    at_mass_A
        Vector of atomic masses of system A.
    at_mass_B
        Vector of atomic masses of system B.

    Returns
    -------
    np.ndarray
        (N,N) Duschinsky matrix J.
    """
    if orth_L:
        if at_mass_A is None and at_mass_B is None:
            Jmat = np.einsum('ij,kj->ik', Lmat_A, Lmat_B)
        else:
            if not (at_mass_A is not None and at_mass_B is not None):
                raise ArgumentError('Atomic masses of both systems are needed')
            massA = np.sqrt(np.repeat(at_mass_A, 3))
            massB = np.repeat(at_mass_B, 3)**(-1/2)
            Jmat = np.einsum('ij,j,j,kj->ij', Lmat_A, massA, massB, Lmat_B)
    else:
        if at_mass_A is None and at_mass_B is None:
            Jmat = np.einsum('ij,jk->ik', Lmat_A, np.linalg.pinv(Lmat_B))
        else:
            if not (at_mass_A is not None and at_mass_B is not None):
                raise ArgumentError('Atomic masses of both systems are needed')
            massA = np.sqrt(np.repeat(at_mass_A, 3))
            massB = np.repeat(at_mass_B, 3)**(-1/2)
            Jmat = np.einsum('ij,j,j,jk->ij', Lmat_A, massA, massB,
                             np.linalg.pinv(Lmat_B))

    return Jmat


def build_dusch_K(Lmat: np.ndarray,
                  at_mass: TypeAtMas,
                  at_deltaR: tp.Optional[TypeAtCrd] = None) -> np.ndarray:
    r"""Build Duschinsky-transformation shift vector K.

    Build the shift vector associated to the Duschinsky transformation,

    :math:`Q_A = J Q_B + K`

    The form of K depends on the availability of two equilibrium
    geometries.  In electronic transitions, this is related to the
    adiabatic or vertical models adopted to describe the transition.

    Parameters
    ----------
    Lmat
        Transfo. matrix from mass-weighted cart. to normal coord. (N,3 Na).
    at_mass
        Atomic masses, as a 1D array.
    at_deltaR
        Difference between equilibrium geometries (:math:`\Delta R`).
    """
    if at_deltaR is not None:
        mass = np.sqrt(np.repeat(at_mass, 3))
        Kvec = np.einsum('ij,j,j->i', Lmat, mass, np.reshape(at_deltaR, (-1,)))
    else:
        raise NotImplementedError('Vertical K not yet implemented.')
    return Kvec


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
        Transformation matrix from cart to normal coord. (N,3Na).

    Returns
    -------
    :obj:`np.ndarray`
        `Lmat` with corrected columns.
    """
    if not isinstance(Lmat, np.ndarray):
        raise IndexError('Wrong parameter, Lmat should be an array')
    shift = 1.0e-6
    if Lmat.ndim == 2:
        nvib, nat3 = Lmat.shape
        for i in range(nvib):
            xmax = abs(Lmat[i, 0]) + shift
            jmax = 0
            for j in range(1, nat3):
                if abs(Lmat[i, j]) > xmax:
                    xmax = abs(Lmat[i, j]) + shift
                    jmax = j
            if Lmat[i, jmax] < 0.:
                Lmat[i, :] *= -1.
    elif Lmat.ndim == 3:
        nvib, nat, xyz = Lmat.shape
        if xyz != 3:
            raise IndexError('Expected 3rd dimension to have 3 elements.')
        for i in range(nvib):
            xmax = abs(Lmat[i, 0, 0]) + shift
            jmax = 0
            kmax = 0
            for j in range(nat):
                for k in range(3):
                    if abs(Lmat[i, j, k]) > xmax:
                        xmax = abs(Lmat[i, j, k]) + shift
                        jmax = j
                        kmax = k
            if Lmat[i, jmax, kmax] < 0.:
                Lmat[i, ...] *= -1.
    else:
        raise IndexError('Lmat is expected to have dimension 2 or 3')

    return Lmat
