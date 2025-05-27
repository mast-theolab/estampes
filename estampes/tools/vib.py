"""Toolbox for molecular vibrations

Module providing tools to manipulate data relative to vibrations.

"""

import sys
from math import exp, sqrt, copysign
import typing as tp

import numpy as np
import numpy.typing as npt

from estampes.base import TypeAtCrd, TypeAtMas
from estampes.base.errors import ArgumentError, QuantityError
from estampes.data.physics import phys_fact
from estampes.tools.mol import eckart_orient


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
            Jmat = np.einsum('ij,j,j,kj->ik', Lmat_A, massA, massB, Lmat_B)
    else:
        if at_mass_A is None and at_mass_B is None:
            Jmat = np.einsum('ij,jk->ik', Lmat_A, np.linalg.pinv(Lmat_B))
        else:
            if not (at_mass_A is not None and at_mass_B is not None):
                raise ArgumentError('Atomic masses of both systems are needed')
            massA = np.sqrt(np.repeat(at_mass_A, 3))
            massB = np.repeat(at_mass_B, 3)**(-1/2)
            Jmat = np.einsum('ij,j,j,jk->ik', Lmat_A, massA, massB,
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
        Atomic masses in unified atomic mass, as a 1D array.
    at_deltaR
        Difference between equilibrium geometries (:math:`\Delta R`).
    """
    if at_deltaR is not None:
        mass = np.sqrt(np.repeat(at_mass/phys_fact('au2amu'), 3))
        Kvec = np.einsum('ij,j,j->i', Lmat, mass, np.reshape(at_deltaR, (-1,)))
    else:
        raise NotImplementedError('Vertical K not yet implemented.')
    return Kvec


def build_vibrations(fc_cart: npt.ArrayLike,
                     at_mass: npt.ArrayLike,
                     at_crd: npt.ArrayLike,
                     fc_not_weighted: bool = True,
                     get_evec: bool = True,
                     get_eval: bool = True,
                     get_rmas: bool = True,
                     get_lweigh: bool = True,
                     nvib: bool = None,
                     remove_rottrans: bool = True
                     ) -> tp.Dict[str, tp.Optional[npt.NDArray]]:
    """Build vibrational data from force constants matrix.

    Diagonalizes the force constants matrix and build the frequencies
    (wavenumbers) and eigenvectors.
    The function takes also care of projecting out the rotations and
    translations.

    Parameters
    ----------
    fc_cart
        Force constants matrix, as a symmetric matrix.
    at_mass
        Atomic masses, as a vector.
    at_crd
        Atomic coordinates as a 2D array-like.
    fc_not_weighted
        The force constants are not weighted by the atomic masses.
    get_evec
        Return the eigenvectors from the diagonalization.
    get_eval
        Return the eigenvalues (wavenumbers) from the diagonalization.
    get_rmas
        Return the reduced masses.
    get_lweigh
        Return the mass-weighted eigenvectors matrix.
    nvib
        If provided, the system checks that the number of modes found
        matches the reference number.  Raises an error otherwise.
    remove_rottrans
        Remove any residual rotation and translation by projecting them
        out.

    Returns
    -------
    dict
        Dictionary with the quantities of interest filled.
    """
    def hessval_to_freq(hessval):
        """Convert the Hessian eigenvalue to wavenumbers.

        Given an eigenvalue from the Hessian matrix, convert to
        wavenumbers, preserving their sign in the final energies.
        """
        return copysign(sqrt(abs(hessval)*phys_fact('fac2au')), hessval)

    result = {'evec': None, 'freq': None, 'redmas': None, 'lmweigh': None}
    # First, quick size check
    n_at = len(at_mass)
    if not isinstance(at_mass, np.ndarray):
        at_mass = np.array(at_mass)
    if not isinstance(fc_cart, np.ndarray):
        fc_cart = np.array(fc_cart)
    if not isinstance(at_crd, np.ndarray):
        at_crd = np.array(at_crd)
    n_at3 = fc_cart.shape[0]
    if n_at*3 != n_at3:
        raise ArgumentError('Size',
                            'Arrays in arguments have inconsistent size')
    # Weigh force constants matrix if needed.
    sqmas_nat3 = np.tile(np.sqrt(at_mass)[:, np.newaxis], (1, 3))
    inv_sqmas = 1.0/np.reshape(sqmas_nat3, (-1, ))
    if fc_not_weighted:
        fc_mweig = np.einsum('i,ij,j->ij', inv_sqmas, fc_cart, inv_sqmas)
    else:
        fc_mweig = fc_cart

    # Compute the eigenvalues and eigenvectors
    hessval, hessvec = np.linalg.eigh(fc_mweig)
    if not remove_rottrans:
        thresh = 8.0
        thresh_incr = 2.0
        max_thresh = 50.0
        lmweig = np.einsum('i,ij->ij', inv_sqmas, hessvec)
        red_mas = 1/(lmweig**2).sum(axis=0)
        while True:
            vibs = np.full(n_at3, True)
            freqs = []
            for i, val in enumerate(hessval):
                freq = hessval_to_freq(val)
                if abs(freq) < thresh:
                    vibs[i] = False
                else:
                    freqs.append(freq)
            if nvib is not None:
                if np.count_nonzero(vibs) != nvib:
                    if thresh < max_thresh:
                        thresh += thresh_incr
                    else:
                        msg = 'Unable to identify vibrations from rot/trans'
                        raise QuantityError(msg)
                else:
                    break
            elif np.count_nonzero(vibs) in (5, 6):
                break
        if get_eval:
            result['freq'] = freqs
        if get_evec:
            result['evec'] = np.transpose(hessvec[:, vibs])
        if get_lweigh:
            result['lmweigh'] = np.transpose(
                lmweig[:, vibs]*np.sqrt(red_mas[vibs]))
        if get_rmas:
            result['redmas'] = red_mas[vibs]
    else:
        # Let us build the internal coordinates FC matrix
        # following the white paper from Gaussian
        # First, we build the reference rotations and translations
        # Compute the center of mass and translates the molecule
        c_eck, p_evec = eckart_orient(at_crd, at_mass, True).values()
        # Second, build the rotation and translations
        new_evec = np.zeros((n_at3, n_at3))
        # Note: the matrix will be built putting eigenvectors as columns
        #       for consistency with eigenvectors given by eigh (in cols.)
        # -- we duplicate p_evec over each atom to facilitate operations
        p_evec_at = np.repeat(p_evec[np.newaxis, :, :], n_at, axis=0)
        # -- Construct translations
        new_evec[:, :3] = np.reshape(
            np.broadcast_to(sqmas_nat3[:, :, np.newaxis],
                            (n_at, 3, 3))*p_evec_at,
            (n_at3, 3))
        new_evec[:, 3:6] = np.reshape(np.cross(
            np.tile((sqmas_nat3*c_eck)[:, np.newaxis, :], (1, 3, 1)),
            p_evec_at),
            (n_at3, 3))
        # -- we check if we have any non-null eigenvalues
        trrot_norm2 = np.einsum('ij,ij->j', new_evec[:, :6], new_evec[:, :6])
        valid_trrot = trrot_norm2 > 1.0e-9
        n_trrot_non0 = sum(valid_trrot)
        if n_trrot_non0 == 5:
            n_trrot = 5
        elif n_trrot_non0 < 5:
            raise ValueError(
                'Too many negligible eigenvalues from translation-rotation')
        else:
            n_trrot = 6
        nvib0 = n_at3 - n_trrot
        if nvib is not None:
            if nvib0 != nvib:
                msg = 'Unable to identify rot/trans'
                raise QuantityError(msg)
        for i in range(6):
            if valid_trrot[i]:
                new_evec[:, i] /= trrot_norm2[i]
        # Third, build the remaining matrix
        # -- 1. look for the overlap between original eigenvectors
        #       and the "clean" translational/rotational vectors
        overlap = ((hessvec.T @ new_evec[:, :6])**2).sum(axis=1)
        # We now sort the result of the overlap
        max_overlap = np.sort(overlap)[-n_trrot]
        new_evec[:, n_trrot:] = hessvec[:, overlap < max_overlap]
        # -- 2. construct the Gram-Schmidt orthogonalization through QR
        evec_orth, _ = np.linalg.qr(new_evec)

        if (np.linalg.norm(evec_orth, axis=0) < 1.0e-6).any():
            raise ValueError(
                'Null vectors found in Gram-Schmidt orthogonalization')
        fc_int = (evec_orth.T @ fc_mweig) @ evec_orth
        hessval, hessvec = np.linalg.eigh(fc_int[n_trrot:, n_trrot:])
        lmat = np.einsum('ij,jk->ik', evec_orth[:, n_trrot:], hessvec)
        lmweig = np.einsum('i,ij->ij', inv_sqmas, lmat)
        red_mas = 1/(lmweig**2).sum(axis=0)
        freqs = []
        for i, val in enumerate(hessval):
            freqs.append(hessval_to_freq(val))
        if get_eval:
            result['freq'] = freqs
        if get_evec:
            result['evec'] = np.transpose(lmat)
        if get_lweigh:
            result['lmweigh'] = np.transpose(lmweig*np.sqrt(red_mas))
        if get_rmas:
            result['redmas'] = red_mas

    return result


def convert_hess_evec(evec: npt.ArrayLike,
                      at_mass: tp.Optional[npt.ArrayLike] = None,
                      natoms: tp.Optional[int] = None,
                      form: str = 'L.M^{-1/2}',
                      do_norm: bool = True) -> np.ndarray:
    """Convert Hessian eigenvectors matrix.

    Converts the eigenvectors matrix from the diagonalization of the
    Hessian matrix, correcting mass weights and fixing the shape if
    necessary.

    Parameters
    ----------
    evec
        Original eigenvectors matrix.
    at_mass
        Atomic masses (in u).
    natoms
        Number of atoms.
    form
        Current form of the eigenvectors, as a ESTAMPES datatype.
    do_norm
        Normalize the final eigenvectors.

    Returns
    -------
    np.ndarray
        L matrix.
    """
    eigvec = None
    if form in ('L.M^-1/2', 'L/M^1/2', 'L.M^{-1/2}', 'L/M^{1/2}'):
        if at_mass is not None:
            nat3 = at_mass.size
            eigvec = np.einsum(
                'ij,j->ij',
                np.reshape(evec, (-1, nat3)),
                np.sqrt(at_mass)
            )
            if do_norm:
                eigvec = norm_evec(eigvec)
        else:
            raise QuantityError('Missing atomic masses to correct evec')
    else:
        if natoms is not None:
            eigvec = np.reshape(evec, (-1, 3*natoms))
        else:
            # Compute nat3 based on: 3*nat*(3nat-ntrro) = size(evec)
            # ntrro = 6 for non-linear, 5 otherwise
            # The positive root should be last one (ascending order)
            # assume first most common case: non-linear
            N = evec.size
            val = np.polynomial.polynomial.polyroots((-N, -6, 1))[-1]
            if val.is_integer():
                nat3 = int(val)
            else:
                nat3 = int(np.polynomial.polynomial.polyroots(
                    (-N, -5, 1))[-1])
            eigvec = np.reshape(evec, (-1, nat3))
    return eigvec


def get_vib_trans(trans_data: tp.Sequence[tp.Any]
                  ) -> tp.Dict[int, int]:
    """Return vibrational transition information.

    Taking a vibrational transition information, returns the transition
    specification.

    Parameters
    ----------
    trans_data
        Transition information as built by ESTAMPES' parsers.

    Returns
    -------
    dict
        Number of quanta for each excited mode/state.
    """
    state_i, state_f, *_ = trans_data
    trans = {}
    is_var = False
    for i, ni in state_i:
        if i > 0:
            if ni == 0:
                if is_var:
                    raise ValueError(
                        'Wrong variational state specifications for starting',
                        ' state')
                is_var = True
                trans[i] = -1
            else:
                trans[i] = -ni
    is_var = False
    for i, ni in state_f:
        if i > 0:
            if ni == 0:
                if is_var:
                    raise ValueError(
                        'Wrong variational state specifications for starting',
                        ' state')
                if i not in trans:
                    trans[i] = 1
                else:
                    trans[i] += 1
            else:
                if i not in trans:
                    trans[i] = ni
                else:
                    trans[i] += ni

    return trans


def norm_evec(evec: npt.ArrayLike) -> np.ndarray:
    """Normalize Hessian eigenvectors matrix.

    Normalizes the matrix of eigenvectors from the Hessian/force
    constants matrix.

    Parameters
    ----------
    evec
        Eigenvectors matrix, assumed to have shape (nvib, nat3).

    Returns
    -------
    np.ndarray
        Matrix of normalized eigenvectors.
    """
    res = np.empty(evec.shape)
    norms = np.sum(evec**2, axis=1)
    for i in range(evec.shape[0]):
        if norms[i] > sys.float_info.epsilon:
            res[i, :] = evec[i, :] / sqrt(norms[i])
        else:
            res[i, :] = 0.0
    return res


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


def weigh_trans_progress(trans_spec: tp.Dict[int, int],
                         energy: tp.Union[tp.Sequence[float],
                                          tp.Dict[int, float]],
                         temp: float = 298.15,
                         i_offset: int = 0) -> float:
    r"""Weigh transition progressions considering all populations.

    Weighs transitions considering all possible transitions involving
    the same modes and the same difference in quanta.
    Note that this method is only true for the harmonic approximation.
    Considering the following:

    .. math::
        Z = \sum_s
            e^{-\sum_{i=1}^N \frac{\hbar \omega_i}{k_B T}
               (n^s_i + \frac{1}{2})}
        = \sum_s \prod_{i=1}^N e^{-\frac{\hbar \omega_i}{k_B T}
                                  (n^s_i + \frac{1}{2})}
        = \sum_{n_1=0}^\infty \sum_{n_2=0}^\infty \ldots
          \sum_{n_N=0}^\infty
          \prod_{i=1}^N e^{-\frac{\hbar \omega_i}{k_B T}
                           (n_i + \frac{1}{2})}

        Z = \prod_{i=1}^N \sum_{n_i=0}^\infty
            e^{-\frac{\hbar \omega_i}{k_B T}(n_i + \frac{1}{2})}
        = \prod_{i=1}^N Z_i

        Z_i = \sum_{n_i=0}^\infty
            e^{-\frac{\hbar \omega_i}{k_B T}(n_i + \frac{1}{2})}
        = \frac{e^{-\frac{\hbar \omega_i}{2 k_B T}}}
               {1 - e^{-\frac{\hbar \omega_i}{k_B T}}}

    Within the harmonic approximation, the probability
    :math:`p_{n_i}` for a mode (:math:`i`) to be excited at a given
    number of quanta :math:`n_i` can be generalized to:

    .. math::
        p_{n_i} = \frac{e^{-\frac{\hbar \omega_i}{k_B T}
                           (n_i + \frac{1}{2})}
                        \prod_{j \neq i}^N \sum_{n_i=0}^\infty
                            e^{-\frac{\hbar \omega_j}{k_B T}
                               (n_j + \frac{1}{2})}
                        }{Z}
        = \frac{e^{-\frac{\hbar \omega_i}{k_B T}(n_i + \frac{1}{2})}}
               {Z_i}

    Now considering all possible excitations from this mode, we would
    have, for an arbitrary intensity, at the harmonic level,

    .. math::
        I = \sum_{n_i=0}^\infty p_{n_i}
            \frac{\partial P_a}{\partial q_i}
            \frac{\partial {P_b}^\ast}{\partial q_i}
            | \langle n_i \mid q_i | n_i + \delta n \rangle |^2

        I = e^{\frac{\hbar \omega_i}{2 k_B T}} \biggl(
            1 - e^{-\frac{\hbar \omega_i}{k_B T}}
        \biggr)
        \frac{\partial P_a}{\partial q_i}
        \frac{\partial {P_b}^\ast}{\partial q_i}
        \sum_{n_i=0}^\infty
            e^{-\frac{\hbar \omega_i}{k_B T}(n_i + \frac{1}{2})}
            | \langle n_i \mid q_i | n_i + \delta n \rangle |^2

        I = \biggl(
            1 - e^{-\frac{\hbar \omega_i}{k_B T}}
        \biggr)
        \frac{\partial P_a}{\partial q_i}
        \frac{\partial {P_b}^\ast}{\partial q_i}
        \sum_{n_i=0}^\infty
            e^{-\frac{\hbar n_i \omega_i}{k_B T}}
            | \langle n_i \mid q_i | n_i + \delta n \rangle |^2

    Standard case for harmonic level, :math:`\delta n = 1`,

    .. math::
        \sum_{n_i=0}^\infty
            e^{-\frac{\hbar n_i \omega_i}{k_B T}}
            (n_i + 1)
        = \frac{e^{\frac{2 \omega_i}{k_B T}}}
               {\biggl(e^{\frac{\omega_i}{k_B T}} - 1\biggr)^2}
        = \frac{1}
               {\biggl(1 - e^{-\frac{\omega_i}{k_B T}}\biggr)^2}

    This can be extended to :math:`\delta n = 2`,

    .. math::
        \sum_{n_i=0}^\infty
            e^{-\frac{\hbar n_i \omega_i}{k_B T}}
            (n_i + 1)(n_i + 2)
        = \frac{2 e^{\frac{3 \omega_i}{k_B T}}}
               {\biggl(e^{\frac{\omega_i}{k_B T}} - 1\biggr)^3}
        = \frac{2}
               {\biggl(1 - e^{-\frac{\omega_i}{k_B T}}\biggr)^3}

    and :math:`\delta n = 3`,

    .. math::
        \sum_{n_i=0}^\infty
            e^{-\frac{\hbar n_i \omega_i}{k_B T}}
            (n_i + 1)(n_i + 2)(n_i + 3)
        = \frac{6 e^{\frac{4 \omega_i}{k_B T}}}
               {\biggl(e^{\frac{\omega_i}{k_B T}} - 1\biggr)^4}
        = \frac{6}
               {\biggl(1 - e^{-\frac{\omega_i}{k_B T}}\biggr)^4}

    Transitions involving multiple modes can be also included
    considering the probability of each involved mode having a set
    number of quanta.

    Parameters
    ----------
    trans_spec
        Transition specification, as: (mode, quanta difference).
    energy
        Mode energies, in cm\ :sup:`-1`.
    temp
        Temperature, in K.
    i_offset
        Offset to apply between mode index and index in `energy`.

    Returns
    -------
    float
        weigh to be applied to transition.
    """
    weight = 1.0
    # kbTm1: 1/(k_B.T) in cm
    kbTm1 = 2.0 * phys_fact('factb') / temp
    for i, ni in trans_spec.items():
        wi = energy[i + i_offset]
        wikbT = wi * kbTm1
        if ni == 1:
            weight *= 1.0 / (1 - exp(-wikbT))
        elif ni == 2:
            weight *= 2.0 / (1 - exp(-wikbT))**2
        elif ni == 3:
            weight *= 6.0 / (1 - exp(-wikbT))**3
        elif ni == -1:
            weight *= exp(-wikbT) / (1 - exp(-wikbT))
        elif ni == -2:
            weight *= 2.0 * exp(-2*wikbT) / (1 - exp(-wikbT))**2
        elif ni == 3:
            weight *= 6.0 * exp(-3* wikbT) / (1 - exp(-wikbT))**3
        else:
            raise ValueError(
                f'Unsupported quanta difference for mode {i}')

    return weight
