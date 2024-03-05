"""Module providing mathematical functions.

A basic module providing methods for some useful mathematical
operations, for ESTAMPES tools.

"""

from math import acos, atan2, degrees, exp, log, pi, sqrt
import typing as tp

import numpy as np
import numpy.typing as npt

from estampes.base import TypeAtCrd, ArgumentError
from estampes.data.physics import PHYSFACT


# ==============
# Module Methods
# ==============


def angle(r1: npt.ArrayLike, r2: npt.ArrayLike, r3: npt.ArrayLike,
          radians: bool = False) -> float:
    """Compute the angle formed by the points 1, 2, 3.

    Computes the angle between points 1, 2, 3 at positions r_i.
    It is assumed that the angle is formed by taking the points
    sequentially.

    Parameters
    ----------
    r1
        Cartesian position of point 1.
    r2
        Cartesian position of point 2.
    r3
        Cartesian position of point 3.
    radians
        Returns angle in radians.

    Returns
    -------
    float
        angle.
    """
    v1 = np.array(r1) - np.array(r2)
    v2 = np.array(r3) - np.array(r2)

    ang = acos(v1@v2/(np.linalg.norm(v1)*np.linalg.norm(v2)))

    if not radians:
        return degrees(ang)
    else:
        return ang


def bond(r1: npt.ArrayLike, r2: npt.ArrayLike, in_au: bool = False) -> float:
    """Compute the bond length between points 1 and 2.

    Computes the distance between points 1 and 2 at positions r_i.
    The positions are assumed in atomic units and returned in Angstroms
    by default for Z-matrix specification for instance.

    Parameters
    ----------
    r1
        Cartesian position of point 1.
    r2
        Cartesian position of point 2.
    in_au
        Keep in atomic units.

    Returns
    -------
    float
        bond length.
    """
    if in_au:
        return np.linalg.norm(np.array(r2) - np.array(r1))
    else:
        return np.linalg.norm(np.array(r2) - np.array(r1))*PHYSFACT.bohr2ang


def dihedral(r1: npt.ArrayLike, r2: npt.ArrayLike, r3: npt.ArrayLike,
             r4: npt.ArrayLike, radians: bool = False) -> float:
    """Compute the dihedral formed between points 1, 2, 3, 4.

    Computes the dihedral angle between points 1, 2, 3, 4 at positions
    ri.  It is assumed that the dihedral is formed by taking the points
    sequentially, with 2, 3 shared between the two planes.

    Parameters
    ----------
    r1
        Cartesian position of point 1.
    r2
        Cartesian position of point 2.
    r3
        Cartesian position of point 3.
    r4
        Cartesian position of point 4.
    radians
        Returns angle in radians.

    Returns
    -------
    float
        Dihedral angle.
    """
    v1 = np.array(r2) - np.array(r1)
    v2 = np.array(r3) - np.array(r2)
    v3 = np.array(r4) - np.array(r3)

    ang = atan2(np.linalg.norm(v2)*v1@np.cross(v2, v3),
                np.cross(v1, v2)@np.cross(v2, v3))
    if not radians:
        return degrees(ang)
    else:
        return ang


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


def levi_civita_tens() -> np.ndarray:
    """Return the Levi-Civita tensor.

    Builds and returns the Levi-Civita tensor.

    Returns
    -------
    :obj:np.ndarray
        3D tensor (3,3,3)
    """
    tensor = np.zeros((3, 3, 3))
    tensor[0, 1, 2] = tensor[1, 2, 0] = tensor[2, 0, 1] = 1
    tensor[0, 2, 1] = tensor[2, 1, 0] = tensor[1, 0, 2] = -1

    return tensor


def square_ltmat(ltmat: tp.Union[tp.Sequence[float], np.ndarray],
                 what: str = 'symm') -> np.ndarray:
    """Square a lower-triangular matrix.

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
    n = int(-.5 + sqrt(.5**2 + 2*np.size(ltmat)))
    if np.size(ltmat) != n*(n+1)//2:
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


def superpose(c_ref: TypeAtCrd,
              c_new: TypeAtCrd,
              at_mass: tp.Optional[np.ndarray] = None,
              get_ctrans: bool = False,
              at_mask: tp.Optional[np.ndarray] = None
              ) -> tp.Union[tp.Tuple[np.ndarray, np.ndarray],
                            tp.Tuple[np.ndarray, np.ndarray, TypeAtCrd]]:
    """Return the transformation matrices to superpose c_new onto c_ref.

    Returns the rotation matrix and transition vector to maximize the
    superposition of `c_new` onto `c_ref`.  The translated and rotated
    coordinates can be returned on request.
    The superposition can be mass-weighted if requested.
    The coordinates should have the form (N,3).

    The algorithm is described in Refs. [1]_ [2]_

    Parameters
    ----------
    c_ref
        Reference structure, as a Numpy 2D array (N,3).
    c_new
        New structure to superpose, as a Numpy 2D array (N,3).
    at_mass
        Atomic masses, as a Numpy 1D array (N).
    get_ctrans
        If True, return the transformed structure.
    at_mask
        Atoms to consider for the superposition scheme.
        The translation still refers to the full system.
        Mask should be a Numpy 1D array(N).

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

    Notes
    -----
    See Ref. [1,2] for details on the algorithms.

    References
    ----------
    .. [1] G.R. Kneller, Mol. Sim. 1991, 7, 113-119
    .. [2] G.R. Kneller, J. Chim. Phys. 1991, 88, 2709-2715
    """
    DEBUG = False
    com_ref = [0., 0., 0.]
    com_new = [0., 0., 0.]
    size_ref = c_ref.shape
    size_new = c_new.shape
    if size_ref != size_new:
        raise IndexError('Inconsistency in shapes of the structure')
    natoms = size_ref[0]
    if at_mass is not None:
        use_m = True
        if at_mass.shape[0] != natoms:
            raise IndexError('Atomic masses inconsistent with structures')
    else:
        use_m = False
    if at_mask is not None:
        mask = at_mask
        if mask.shape[0] != natoms:
            raise IndexError('Mask inconsistent with structures.')
    else:
        mask = np.full(natoms, True)
    rotmat = np.identity(3)
    qmat = np.zeros((4, 4))

    # Calculate the center of mass and move the structure
    if use_m:
        totwt = np.sum(at_mass[mask])
        com_ref = np.einsum('ij,i->j', c_ref[mask, :], at_mass[mask])
        com_new = np.einsum('ij,i->j', c_new[mask, :], at_mass[mask])
    else:
        totwt = np.sum(at_mass[mask])
        com_ref = np.einsum('ij->j', c_ref[mask, :])
        com_new = np.einsum('ij->j', c_new[mask, :])
    c_new_ = c_new - com_new/totwt
    c_ref_ = c_ref - com_ref/totwt

    # Computes the rotation matrix
    for ia in range(natoms):
        if not mask[ia]:
            continue
        weight = at_mass[ia] if use_m else 1.0
        xnew = c_new_[ia, 0]
        ynew = c_new_[ia, 1]
        znew = c_new_[ia, 2]
        xref = c_ref_[ia, 0]
        yref = c_ref_[ia, 1]
        zref = c_ref_[ia, 2]
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
        return rotmat, (com_new-com_ref)/totwt, c_new_ @ rotmat
    else:
        return rotmat, (com_new-com_ref)/totwt


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
