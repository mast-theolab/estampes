"""Module providing mathematical functions

A basic module providing methods for some useful mathematical
  operations, for ESTAMPES tools.

Methods
-------
f_gauss
    Gaussian distribution function.
f_lorentz
    Lorentzian distribution function.
square_ltmat
    Squares a lower-triangular matrix.
superpose
    Returns the transformation matrices to superpose 2 structures.
vrotate_3D
    Rotates a vector in a 3D space.
"""

from math import exp, log, pi, sqrt
import typing as tp

import numpy as np

from estampes.base import ArgumentError


# ==============
# Module Methods
# ==============

def f_gauss(x: float,
            hwhm: float = 2.0,
            x0: float = 0.0,
            y0: float = 1.0) -> float:
    """Value of a Gaussian broadening function at x.

    Computes the value of Gaussian distribution function with an
        integrated area of `y0` centered on position `x0` at `x`.

    Parameters
    ----------
    x
        Position at which the function must be computed.
    hwhm
        Half-width at half maximum.
    x0
        Center of the function.
    y0
        Total, integrated intensity

    Returns
    -------
    float
        Value of the Gaussian function in position x.

    Raises
    ------
    ValueError
        Wrong value for the HWHM.
    """
    if hwhm <= 0.0:
        raise ValueError('Wrong HWHM.')
    sigma2 = hwhm**2/log(2)  # technically, sigma2 = 2*sigma
    # Normalization factor
    norm = 1./sqrt(sigma2*pi)
    return y0*norm*exp(-(x-x0)**2/sigma2)


def f_lorentz(x: float,
              hwhm: float = 2.0,
              x0: float = 0.0,
              y0: float = 1.0) -> float:
    """Value of a Lorentzian broadening function at x.

    Computes the value of Lorentzian distribution function with an
        integrated area of `y0` centered on position `x0` at `x`.

    Parameters
    ----------
    x
        Position at which the function must be computed.
    hwhm
        Half-width at half maximum.
    x0
        Center of the function.
    y0
        Total, integrated intensity

    Returns
    -------
    float
        Value of the Lorentzian function in position x.

    Raises
    ------
    ValueError
        Wrong value for the HWHM.
    """
    norm = hwhm/pi
    return y0*norm/((x-x0)**2+hwhm**2)


def square_ltmat(ltmat: np.ndarray, what: str = 'symm') -> np.ndarray:
    """Squares a lower-triangular matrix.

    Takes a lower-triangular matrix and returns the 2D symmetric or
      antisymmetric matrix.

    Parameters
    ----------
    ltmat
        Lower-triangular matrix as 1D Numpy array.
    what
        Operation to do:
        * `symm`: symmetrizes
        * `anti`: antisymmetrizes

    Returns
    -------
    :obj:np.ndarray
        2D matrix.

    Raises
    ------
    ValueError
        Input matrix does not seem to have lower-triangular size.
    """
    n = int((-.5 + sqrt(.5**2 + 2*ltmat.size))/1)
    if ltmat.size != n*(n+1)//2:
        raise ValueError('Inconsistency size in LT matrix')
    sqmat = np.zeros((n, n))
    ij = 0
    if what == 'symm':
        for i in range(n):
            for j in range(i):
                sqmat[i, j] = ltmat[ij]
                sqmat[j, i] = ltmat[ij]
                ij += 1
            sqmat[i, i] = ltmat[ij]
            ij += 1
    elif what == 'anti':
        for i in range(n):
            for j in range(i):
                sqmat[i, j] = ltmat[ij]
                sqmat[j, i] = -ltmat[ij]
                ij += 1
            sqmat[i, i] = ltmat[ij]
            ij += 1
    else:
        raise ArgumentError('Unrecognized symmetrization option')
    return sqmat


def superpose(cref: np.ndarray,
              cnew: np.ndarray,
              atmass: tp.Optional[np.ndarray] = None,
              get_ctrans: bool = False
              ) -> tp.Union[tp.Tuple[np.ndarray, np.ndarray],
                            tp.Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Returns the transformation matrices to superpose cnew onto cref.

    Returns the rotation matrix and transition vector to maximize the
      superposition of `cnew` onto `cref`.  The translated and rotated
      coordinates can be returned on request.
    The superposition can be mass-weighted if requested.
    The coordinates should have the form (N,3).

    Parameters
    ----------
    cref
        Reference structure, as a Numpy 2D array (N,3).
    cnew
        New structure to superpose, as a Numpy 2D array (N,3).
    atmass
        Atomic masses, as a Numpy 1D array (N).
    get_ctrans
        If True, return the transformed structure.

    Returns
    -------
    np.ndarray
        Rotation matrix (3,3).
    np.ndarray
        Transition vector (3).
    np.ndarray : optional
        Transformed structure, on request (N,3).

    Raises
    ------
    IndexError
        Inconsistency in structure shapes.

    References
    ----------
    .. [1] G.R. Kneller, Mol. Sim. 1991, 7, 113-119
    .. [2] G.R. Kneller, J. Chim. Phys. 1991, 88, 2709-2715
    """
    DEBUG = False
    com_ref = [0., 0., 0.]
    com_new = [0., 0., 0.]
    size_ref = cref.shape
    size_new = cnew.shape
    if size_ref != size_new:
        raise IndexError('Inconsistency in shapes of the structure')
    natoms = size_ref[0]
    if atmass is not None:
        use_m = True
        if atmass.shape[0] != natoms:
            raise IndexError('Atomic masses inconsistent with structures')
    else:
        use_m = False
    rotmat = np.identity(3)
    qmat = np.zeros((4, 4))

    # Calculate the center of mass and move the structure
    if use_m:
        totwt = np.sum(atmass)
        com_ref = np.einsum('ij,i->j', cref, atmass)
        com_new = np.einsum('ij,i->j', cnew, atmass)
    else:
        totwt = 1.0
        com_ref = np.einsum('ij->j', cref)
        com_new = np.einsum('ij->j', cnew)
    cnew_ = cnew - com_new/totwt
    cref_ = cref - com_ref/totwt

    # Computes the rotation matrix
    for ia in range(natoms):
        weight = atmass[ia] if use_m else 1.0
        xnew = cnew_[ia, 0]
        ynew = cnew_[ia, 1]
        znew = cnew_[ia, 2]
        xref = cref_[ia, 0]
        yref = cref_[ia, 1]
        zref = cref_[ia, 2]
        xy = xnew*yref
        xz = xnew*zref
        yx = ynew*xref
        yz = ynew*zref
        zx = znew*xref
        zy = znew*yref
        diag0 = xnew*xnew + ynew*ynew + znew*znew + \
            xref*xref + yref*yref + zref*zref
        diag1 = 2.0*xnew*xref
        diag2 = 2.0*ynew*yref
        diag3 = 2.0*znew*zref
        qmat[0, 0] += weight*(diag0 - diag1 - diag2 - diag3)
        qmat[1, 1] += weight*(diag0 - diag1 + diag2 + diag3)
        qmat[2, 2] += weight*(diag0 + diag1 - diag2 + diag3)
        qmat[3, 3] += weight*(diag0 + diag1 + diag2 - diag3)
        qmat[1, 0] += weight*2.0*(yz - zy)
        qmat[2, 0] += weight*2.0*(zx - xz)
        qmat[3, 0] += weight*2.0*(xy - yx)
        qmat[2, 1] += weight*2.0*(-(xy + yx))
        qmat[3, 1] += weight*2.0*(-(xz + zx))
        qmat[3, 2] += weight*2.0*(-(yz + zy))
    _, qmvec = np.linalg.eigh(qmat, UPLO='L')
    q0q0 = qmvec[0, 0] * qmvec[0, 0]
    q1q1 = qmvec[1, 0] * qmvec[1, 0]
    q2q2 = qmvec[2, 0] * qmvec[2, 0]
    q3q3 = qmvec[3, 0] * qmvec[3, 0]
    rotmat[0, 0] = q0q0 + q1q1 - q2q2 - q3q3
    rotmat[1, 1] = q0q0 - q1q1 + q2q2 - q3q3
    rotmat[2, 2] = q0q0 - q1q1 - q2q2 + q3q3
    q0q1 = qmvec[0, 0] * qmvec[1, 0]
    q0q2 = qmvec[0, 0] * qmvec[2, 0]
    q0q3 = qmvec[0, 0] * qmvec[3, 0]
    q1q2 = qmvec[1, 0] * qmvec[2, 0]
    q1q3 = qmvec[1, 0] * qmvec[3, 0]
    q2q3 = qmvec[2, 0] * qmvec[3, 0]
    rotmat[0, 1] = 2.0 * (q1q2 - q0q3)
    rotmat[1, 0] = 2.0 * (q1q2 + q0q3)
    rotmat[0, 2] = 2.0 * (q1q3 + q0q2)
    rotmat[2, 0] = 2.0 * (q1q3 - q0q2)
    rotmat[1, 2] = 2.0 * (q2q3 - q0q1)
    rotmat[2, 1] = 2.0 * (q2q3 + q0q1)
    if DEBUG:
        print('ROTATION MATRIX')
        for i in range(3):
            print('{0[0]:8.5f}{0[1]:8.5f}{0[2]:8.5f}'.format(rotmat[i, :]))
    if get_ctrans:
        return rotmat, com_new-com_ref, cnew_ @ rotmat
    else:
        return rotmat, com_new-com_ref


def vrotate_3D(vec: np.ndarray,
               ref: np.ndarray) -> np.ndarray:
    """Rotates a vector in a 3D space.

    Returns the rotation matrix for `vec` to match the orientation of a
      reference vector `ref`.
    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311

    Parameters
    ----------
    vec
        Vector to rotate, as a numpy 1D array
    ref
        Reference vector, as a numpy 1D array

    Returns
    -------
    np.ndarray
        (3,3) rotation matrix, as a numpy 2D array
    """
    def norm(A):
        return sqrt(np.dot(A, A))
    # G = np.matrix([
    #     [np.dot(A, B), -norm(np.cross(A, B)), 0.0],
    #     [norm(np.cross(A, B)), np.dot(A, B),  0.0],
    #     [0.0, 0.0, 1.0]
    # ])
    # F = np.matrix([
    #     A,
    #     (B-np.dot(A, B)*A)/norm(B-np.dot(A, B)*A),
    #     np.cross(B, A)/norm(np.cross(B, A))
    # ])
    # return F.I*G*F
    V = np.cross(vec, ref)
    S = norm(V)
    if abs(S) < 1.0e-6:
        # Already collinear, nothing to do
        return np.eye(3)
    else:
        C = np.dot(vec, ref)
        Vx = np.matrix([[0.0, -V[2], V[1]],
                        [V[2], 0.0, -V[0]],
                        [-V[1], V[0], 0.0]])
        return np.eye(3) + Vx + Vx**2*(1.0-C)/S**2
