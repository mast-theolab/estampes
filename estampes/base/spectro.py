"""Module providing instruments to handle/generate spectral observables

A basic module providing classes and methods to facilitate the handling
    of spectral observables

Attributes
----------
RAMAN_SETUPS
    Provides a list of accepted setups to be used in `raman_intensities`
    for instance.

Classes
-------
RamanInvariants
    Build the Raman/ROA invariants.

Methods
-------
raman_intensities
    Compute the Raman/ROA intensities or activities.

TODO
!! Add check of consistent types between numbers
"""

from math import pi
import sys
import re
import typing as tp

import numpy as np

from estampes.tools.math import levi_civita_tens
from estampes.data.physics import PHYSCNST, PHYSFACT
from estampes.base.types import TypeRlCx


# =================
# Module Attributes
# =================

# Speed of light, in atomic units
_c_au = 1.0/PHYSCNST.finestruct

TypeTensor2D = tp.Union[np.ndarray, tp.Sequence[tp.Sequence[TypeRlCx]]]
TypeTensor3D = tp.Union[np.ndarray,
                        tp.Sequence[tp.Sequence[tp.Sequence[TypeRlCx]]]]

RAMAN_SETUPS = {
    'full': (  # full version
        'ICP(0)', 'ICP(90)x', 'ICP(90)z', 'ICP(90)u', 'ICP(90)*', 'ICP(180)',
        'SCP(0)', 'SCP(90)x', 'SCP(90)z', 'SCP(90)u', 'SCP(90)*', 'SCP(180)',
        'DCP_I(0)', 'DCP_I(90)x', 'DCP_I(90)z', 'DCP_I(90)u', 'DCP_I(90)*',
        'DCP_I(180)',
        'DCP_II(0)', 'DCP_II(90)x', 'DCP_II(90)z', 'DCP_II(90)u',
        'DCP_II(90)*', 'DCP_II(180)'
    ),
    'short': (  # Shortened version
        'ICP0', 'ICP90x', 'ICP90z', 'ICP90u', 'ICP90*', 'ICP180',
        'SCP0', 'SCP90x', 'SCP90z', 'SCP90u', 'SCP90*', 'SCP180',
        'DCPI0', 'DCPI90x', 'DCPI90z', 'DCPI90u', 'DCPI90*', 'DCPI180',
        'DCPII0', 'DCPII90x', 'DCPII90z', 'DCPII90u', 'DCPII90*', 'DCPII180'
    )
}


# ==============
# Module Classes
# ==============

# Exception classes
class InvariantNA(Exception):
    """Error for unavailable invariant.

    Error generated were data are missing to compute a particular
      invariant.

    Parameters
    ----------
    name : str
        Name of the invariant
    """
    def __init__(self, name: str) -> None:
        msg = f'Invariant {name} is not available.'
        super().__init__(msg)


# Main classes

class RamanInvariants():
    """Build Raman/ROA invariants

    Builds the Raman/ROA invariants, based on available quantities.
    The labels of the invariants are based on Ref. [1]_.
    The class supports the far-from-resonance case, as well as the more
      general case, in complex.

    Parameters
    ----------
    alpha
        alpha or electric dipole-electric dipole tensor.
        Assumed complete.  If stored lower-triangular, use `set_alpha`.
    indedip_mdip
        electric dipole-magnetic dipole tensor (induced electric dipole).
    G_tensor
        Alias for `indedip_mdip`
    indedip_equa
        electric dipole-electric quadrupole (induced electric dipole).
    A_tensor
        Alias for `indedip_equa`
    edip_indmdip
        electric dipole-magnetic dipole tensor (induced magnetic dipole).
    edip_indequa
        electric dipole-electric quadrupole (induced electric quadrupole).
    E_inc
        Incident energy, in Hartree.
    E_scat
        Scattered energy, in Hartree.
    assume_FFR
        Assume that the far-from-resonance approximation holds.
        Only unique terms are computed.
    is_complex
        If True, properties are assumed complexes.
        By default, this is inferred from the input quantity.
            It cannot be changed later.

    Attributes
    ----------
    alpha
        Electric dipole-electric dipole tensor
    G_indED
        Induced electric dipole-magnetic dipole tensor
    Gtensor
        Alias for `G_indED`
    A_indED
        Induced electric dipole-electric quadrupole tensor
    Atensor
        Alias for "A_indED"
    G_indMD
        Electric dipole-induced magnetic dipole tensor
    A_indEQ
        Electric dipole-induced electric quadrupole tensor
    inc_frq
        Incident frequency (s^-1)
    scat_frq
        Scattered frequency (s^-1)
    alpha2
        Invariant alpha^2
    b_s_alph
        Invariant beta_s(alpha)^2
    b_a_alph
        Invariant beta_a(alpha)^2
    alp_romG
        Invariant alpha.G{roman}
    b_s_romG
        Invariant beta_s(G{roman})
    b_a_romG
        Invariant beta_a(G{roman})
    alp_calG
        Invariant alpha.G{script}
    b_s_calG
        Invariant beta_s(G{script})
    b_a_calG
        Invariant beta_a(G{script})
    b_s_romA
        Invariant beta_s(A{roman})
    b_a_romA
        Invariant beta_a(A{roman})
    b_s_calA
        Invariant beta_s(A{script})
    b_a_calA
        Invariant beta_a(A{script})

    Raises
    ------
    InvariantNA
        Missing quantities to compute the invariant.

    Notes
    -----
    * All units should in atomic units.
    * The polarizability is expected to be volumic (bohr^3).
    * All tensors should be stored with the right structure.

    References
    ----------
    .. [1] L.Hecht, L. Nafie, Mol. Phys. 1991, 72, 441.
    .. [2] L. Nafie, Chem. Phys. 1996, 205, 309.
    """
    def __init__(self,
                 alpha: TypeTensor2D,
                 *,
                 indedip_mdip: tp.Optional[TypeTensor2D] = None,
                 G_tensor: tp.Optional[TypeTensor2D] = None,
                 indedip_equa: tp.Optional[TypeTensor3D] = None,
                 A_tensor: tp.Optional[TypeTensor3D] = None,
                 edip_indmdip: tp.Optional[TypeTensor2D] = None,
                 edip_indequa: tp.Optional[TypeTensor3D] = None,
                 E_inc: tp.Optional[TypeRlCx] = None,
                 E_scat: tp.Optional[TypeRlCx] = None,
                 is_complex: tp.Optional[bool] = None):
        self.__alph_s = None
        self.__alph_a = None
        self.__romG_s = None
        self.__romG_a = None
        self.__romA1_s = None
        self.__romA1_a = None
        self.__romA2_a = None
        self.__calG_s = None
        self.__calG_a = None
        self.__calA1_s = None
        self.__calA1_a = None
        self.__calA2_a = None

        self.alpha = alpha
        if is_complex is None:
            self.__cmplx = np.iscomplex(self.alpha)
        else:
            self.__cmplx = is_complex

        if indedip_mdip is not None and G_tensor is not None:
            raise ValueError('Overlapping definition of the G tensor')
        elif indedip_mdip is not None:
            self.G_indED = indedip_mdip
        elif G_tensor is not None:
            self.G_indED = G_tensor
        else:
            self.__G_indED = None
        if indedip_equa is not None and A_tensor is not None:
            raise ValueError('Overlapping definition of the A tensor')
        elif indedip_equa is not None:
            self.A_indED = indedip_equa
        elif A_tensor is not None:
            self.A_indED = A_tensor
        else:
            self.__A_indED = None
        if edip_indmdip is not None:
            self.G_indMD = edip_indmdip
        else:
            self.__G_indMD = None
        if edip_indequa is not None:
            self.A_indEQ = edip_indequa
        else:
            self.__A_indEQ = None

        self.inc_frq = E_inc
        self.scat_frq = E_scat

    # Invariants
    @property
    def alpha2(self):
        """Compute/return alpha^2"""
        if self.__alpha2 is None:
            self.__alpha2 = np.real(
                np.einsum('ii,jj->', self.__alph_s, self.__alph_s.conj())
            )/9.0
        return self.__alpha2

    @property
    def b_s_alph(self):
        """Compute/return beta_s(alpha)^2"""
        if self.__b_s_alph is None:
            self.__b_s_alph = np.real(
                3*np.einsum('ij,ij->', self.__alph_s, self.__alph_s.conj())
                - np.einsum('ii,jj->', self.__alph_s, self.__alph_s.conj())
            )/2.0
        return self.__b_s_alph

    @property
    def b_a_alph(self):
        """Compute/return beta_a(alpha)^2"""
        if self.__b_a_alph is None:
            self.__b_a_alph = 3.0*np.real(
                np.einsum('ij,ij->', self.__alph_a, self.__alph_a.conj())
            )/2.0
        return self.__b_a_alph

    @property
    def alp_romG(self):
        """Compute/return alpha.G{roman}"""
        if self.__alp_romG is None:
            if self.__G_indED is None:
                raise InvariantNA('alpha.romG')
            self.__alp_romG = np.imag(
                np.einsum('ii,jj->', self.__alph_s, self.__romG_s.conj())
            )/(9.0*_c_au)
        return self.__alp_romG

    @property
    def b_s_romG(self):
        """Compute/return beta_s(G{roman})"""
        if self.__b_s_romG is None:
            if self.__G_indED is None:
                raise InvariantNA('beta_s(romG)')
            self.__b_s_romG = np.imag(
                3*np.einsum('ij,ij->', self.__alph_s, self.__romG_s.conj())
                - np.einsum('ii,jj->', self.__alph_s, self.__romG_s.conj())
            )/(2.0*_c_au)
        return self.__b_s_romG

    @property
    def b_a_romG(self):
        """Compute/return beta_a(G{roman})"""
        if self.__b_a_romG is None:
            if self.__G_indED is None:
                raise InvariantNA('beta_a(romG)')
            self.__b_a_romG = 3.0*np.imag(
                    np.einsum('ij,ij->', self.__alph_a, self.__romG_a.conj())
                )/(2.0*_c_au)
        return self.__b_a_romG

    @property
    def alp_calG(self):
        """Compute/return alpha.G{script}"""
        if self.__alp_calG is None:
            if self.__G_indMD is None:
                raise InvariantNA('alpha.calG')
            self.__alp_calG = np.imag(
                np.einsum('ii,jj->', self.__alph_s, self.__calG_s.conj())
            )/(9.0*_c_au)
        return self.__alp_calG

    @property
    def b_s_calG(self):
        """Compute/return beta_s(G{script})"""
        if self.__b_s_calG is None:
            if self.__G_indMD is None:
                raise InvariantNA('beta_s(calG)')
            self.__b_s_calG = np.imag(
                3*np.einsum('ij,ij->', self.__alph_s, self.__calG_s.conj())
                - np.einsum('ii,jj->', self.__alph_s, self.__calG_s.conj())
            )/(2.0*_c_au)
        return self.__b_s_calG

    @property
    def b_a_calG(self):
        """Compute/return beta_a(G{script})"""
        if self.__b_a_calG is None:
            if self.__G_indMD is None:
                raise InvariantNA('beta_a(calG)')
            self.__b_a_calG = 3.0*np.imag(
                    np.einsum('ij,ij->', self.__alph_a, self.__calG_a.conj())
                )/(2.0*_c_au)
        return self.__b_a_calG

    @property
    def b_s_romA(self):
        """Compute/return beta_s(A{roman})"""
        if self.__b_s_romA is None:
            if self.__A_indED is None:
                raise InvariantNA('beta_s(romA)')
            self.__b_s_romA = self.inc_frq*np.imag(
                1j*np.einsum('ij,ij->', self.__alph_s, self.__romA1_s.conj())
            )/(2.0*_c_au)
        return self.__b_s_romA

    @property
    def b_a_romA(self):
        """Compute/return beta_a(A{roman})"""
        if self.__b_a_romA is None:
            if self.__A_indED is None:
                raise InvariantNA('beta_a(romA)')
            self.__b_a_romA = self.inc_frq*np.imag(
                1j*np.einsum('ij,ij->', self.__alph_a,
                             self.__romA1_a.conj() + self.__romA2_a.conj())
            )/(2.0*_c_au)
        return self.__b_a_romA

    @property
    def b_s_calA(self):
        """Compute/return beta_s(A{script})"""
        if self.__b_s_calA is None:
            if self.__A_indEQ is None:
                raise InvariantNA('beta_s(calA)')
            self.__b_s_calA = self.scat_frq*np.imag(
                1j*np.einsum('ij,ij->', self.__alph_s, self.__calA1_s.conj())
            )/(2.0*_c_au)
        return self.__b_s_calA

    @property
    def b_a_calA(self):
        """Compute/return beta_a(A{script})"""
        if self.__b_a_calA is None:
            if self.__A_indEQ is None:
                raise InvariantNA('beta_a(calA)')
            self.__b_a_calA = self.inc_frq*np.imag(
                1j*np.einsum('ij,ij->', self.__alph_a,
                             (self.__calA1_a.conj() + self.__calA2_a.conj()))
            )/(2.0*_c_au)
        return self.__b_a_calA

    # Encapsulation methods
    def __set_alpha(self, tens: TypeTensor2D):
        """Set electric dipole-electric dipole tensor."""
        if isinstance(tens, np.ndarray):
            self.__alpha = tens
        elif isinstance(tens, (tuple, list)):
            if len(tens) != 3:
                raise NotImplementedError('Non 2D tensors NYI')
            else:
                self.__alpha = np.array(tens)
        else:
            raise NotImplementedError('Non-sequence tensor NYI.')
        self.__alph_s = (self.alpha + self.alpha.T)/2.0
        self.__alph_a = (self.alpha - self.alpha.T)/2.0
        self.__alpha2 = None
        self.__b_s_alph = None
        self.__b_a_alph = None
        self.__alp_romG = None
        self.__b_s_romG = None
        self.__b_a_romG = None
        self.__alp_calG = None
        self.__b_s_calG = None
        self.__b_a_calG = None
        self.__b_s_romA = None
        self.__b_a_romA = None
        self.__b_s_calA = None
        self.__b_a_calA = None

    def __get_alpha(self) -> np.ndarray:
        """Get electric dipole-electric dipole tensor."""
        return self.__alpha

    alpha = property(__get_alpha, __set_alpha)

    def __set_G_indED(self, tens: TypeTensor2D):
        """Set G tensor (induced electric dipole).

        Sets the electric dipole-magnetic dipole (induced electric
            dipole).
        This would be the default for G.
        """
        if isinstance(tens, np.ndarray):
            self.__G_indED = tens
        elif isinstance(tens, (tuple, list)):
            if len(tens) != 3:
                raise NotImplementedError('Non 2D tensors NYI')
            else:
                self.__G_indED = np.array(tens)
        else:
            raise NotImplementedError('Non-sequence tensor NYI.')
        self.__romG_s = (self.__G_indED + self.__G_indED.T)/2.0
        self.__romG_a = (self.__G_indED - self.__G_indED.T)/2.0
        self.__alp_romG = None
        self.__b_s_romG = None
        self.__b_a_romG = None

    def __get_G_indED(self) -> tp.Optional[np.ndarray]:
        """Get G tensor (induced electric dipole)."""
        return self.__G_indED

    G_indED = property(__get_G_indED, __set_G_indED)
    Gtensor = property(__get_G_indED, __set_G_indED)

    def __set_A_indED(self, tens: TypeTensor3D):
        """Set A tensor (induced electric dipole).

        Sets the electric dipole-electric quadrupole (induced electric
            dipole).
        This would be the default for A.
        """
        if isinstance(tens, np.ndarray):
            self.__A_indED = tens
        elif isinstance(tens, (tuple, list)):
            if len(tens) != 3:
                raise NotImplementedError('Non 3D tensors NYI')
            else:
                self.__A_indED = np.array(tens)
        else:
            raise NotImplementedError('Non-sequence tensor NYI.')
        eps = levi_civita_tens()
        epsA = np.einsum('ikl,klj->ij', eps, self.__A_indED)
        self.__romA1_s = (epsA + epsA.T)/2.0
        self.__romA1_a = (epsA - epsA.T)/2.0
        epsA = np.einsum('ijk,lkl->ij', eps, self.__A_indED)
        self.__romA2_a = (epsA - epsA.T)/2.0
        self.__b_s_romA = None
        self.__b_a_romA = None

    def __get_A_indED(self) -> tp.Optional[np.ndarray]:
        """Get A tensor (induced electric dipole)."""
        return self.__A_indED

    A_indED = property(__get_A_indED, __set_A_indED)
    Atensor = property(__get_A_indED, __set_A_indED)

    def __set_G_indMD(self, tens: TypeTensor2D):
        """Set G tensor (induced magnetic dipole).

        Sets the electric dipole-magnetic dipole (induced magnetic
            dipole).
        """
        if isinstance(tens, np.ndarray):
            self.__G_indMD = tens
        elif isinstance(tens, (tuple, list)):
            if len(tens) != 3:
                raise NotImplementedError('Non 2D tensors NYI')
            else:
                self.__G_indMD = np.array(tens)
        else:
            raise NotImplementedError('Non-sequence tensor NYI.')
        self.__calG_s = (self.__G_indMD + self.__G_indMD.T)/2.0
        self.__calG_a = (self.__G_indMD - self.__G_indMD.T)/2.0
        self.__alp_calG = None
        self.__b_s_calG = None
        self.__b_a_calG = None

    def __get_G_indMD(self) -> tp.Optional[np.ndarray]:
        """Get G tensor (induced magnetic dipole)."""
        return self.__G_indMD

    G_indMD = property(__get_G_indMD, __set_G_indMD)

    def __set_A_indEQ(self, tens: TypeTensor3D):
        """Set A tensor (induced electric quadrupole).

        Sets the electric dipole-electric quadrupole (induced electric
            quadrupole).
        """
        if isinstance(tens, np.ndarray):
            self.__A_indEQ = tens
        elif isinstance(tens, (tuple, list)):
            if len(tens) != 3:
                raise NotImplementedError('Non 3D tensors NYI')
            else:
                self.__A_indEQ = np.array(tens)
        else:
            raise NotImplementedError('Non-sequence tensor NYI.')
        eps = levi_civita_tens()
        epsA = np.einsum('ikl,jkl->ij', eps, self.__A_indEQ)
        self.__calA1_s = (epsA + epsA.T)/2.0
        self.__calA1_a = (epsA - epsA.T)/2.0
        # epsA = np.einsum('ijk,llk->ij', eps, self.__A_indEQ)
        epsA = np.einsum('ijk,lkl->ij', eps, self.__A_indEQ)
        self.__calA2_a = (epsA - epsA.T)/2.0
        self.__b_s_calA = None
        self.__b_a_calA = None

    def __get_A_indEQ(self) -> tp.Optional[np.ndarray]:
        """Get A tensor (induced electric quadrupole)."""
        return self.__A_indEQ

    A_indEQ = property(__get_A_indEQ, __set_A_indEQ)

    def __set_incfrq(self, energy: tp.Optional[float]):
        """Set the incident frequency.

        The incident energy is transformed into a frequency internally.
        """
        if energy is not None:
            self.__Winc = 2*pi*_c_au*1.0e-8*PHYSFACT.bohr2ang*energy
        else:
            self.__Winc = None
        self.__b_s_romA = None
        self.__b_a_romA = None

    def __get_incfrq(self) -> tp.Optional[float]:
        return self.__Winc

    inc_frq = property(__get_incfrq, __set_incfrq)

    def __set_scatfrq(self, energy: tp.Optional[float]):
        """Set the scattered frequency.

        The scattered energy is transformed into a frequency internally.
        """
        if energy is not None:
            self.__Wscat = 2*pi*_c_au*1.0e-8*PHYSFACT.bohr2ang*energy
        else:
            self.__Wscat = None
        self.__b_s_calA = None
        self.__b_a_calA = None

    def __get_scatfrq(self) -> tp.Optional[float]:
        return self.__Wscat

    scat_frq = property(__get_scatfrq, __set_scatfrq)


# ==============
# Module Methods
# ==============

def raman_intensities(rinv: RamanInvariants,
                      setup: str = 'SCP(180)u',
                      do_ROA: bool = False,
                      do_FFR: bool = False,
                      use_CID: bool = False,
                      get_activity: bool = False) -> float:
    """Compute Raman/ROA intensities

    Computes and return the Raman or ROA intensities.

    Parameters
    ----------
    rinv
        Raman invariants object.
    setup
        Raman setup, as a string.
    do_ROA
        Compute the ROA intensity instead of the Raman scattering.
    do_FFR
        Use simplified expressions assuming the far-from-resonance regime.
    use_CID
        Use the dimensionless version of ROA.
    get_activity
        Return the Raman/ROA activity instead of the intensity.

    Raises
    ------
    InvariantNA
        Intensity not available because of missing invariants.
    IndexError
        Unrecognized Raman setup.
    """
    res = re.match(r'(?P<mode>ICP|SCP|DCP_?I{1,2})'
                   + r'(?P<scat>\(?90\)?[xzu*]|\(?0\)?|\(?180\)?)',
                   setup,
                   re.I)
    if not res:
        raise IndexError('Unrecognized Raman setup')
    mode = res.group('mode').upper().replace('_', '')
    scat = res.group('scat').lower().replace(')', '').replace('(', '')
    do_CIS = do_ROA is False or (do_ROA and use_CID)
    do_CID = do_ROA
    if do_FFR:
        if mode == 'ICP':
            if scat == '0':
                # CID: 8K/c ( 90 alp_romG + 2 b_s_romG - 2 b_s_romA )
                # CIS: 4K ( 45 alpha2 + 7 b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 90*rinv.alp_romG + 2*rinv.b_s_romG \
                        - 2*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph
            elif scat == '90x':
                # CID: 4K/c ( 45 alp_romG + 7 b_s_romG + b_s_romA )
                # CIS: 2K ( 45 alpha2 + 7 b_s_alph )
                if do_CID:
                    a_CID = 4.
                    x_CID = 45*rinv.alp_romG + 7*rinv.b_s_romG \
                        + rinv.b_s_romA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph
            elif scat == '90z':
                # CID: 8K/c ( 3 b_s_romG - b_s_romA )
                # CIS: 4K ( 3 b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 3*rinv.b_s_romG - rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 3*rinv.b_s_alph
            elif scat == '90u':
                # CID: 4K/c ( 45 alp_romG + 13 b_s_romG - b_s_romA )
                # CIS: 4K ( 45 alpha2 + 13 b_s_alph )
                if do_CID:
                    a_CID = 4.
                    x_CID = 45*rinv.alp_romG + 13*rinv.b_s_romG \
                        - rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 13*rinv.b_s_alph
            elif scat == '90*':
                # CID: 40K/3c ( 9 alp_romG + 2 b_s_romG )
                # CIS: 20K/3 ( 9 alpha2 + 2 b_s_alph )
                if do_CID:
                    a_CID = 40. / 3.
                    x_CID = 9*rinv.alp_romG + 2*rinv.b_s_romG
                if do_CIS:
                    a_CIS = 20. / 3.
                    x_CIS = 9*rinv.alpha2 + 2*rinv.b_s_alph
            elif scat[:3] == '180':
                # CID: 8K/c ( 12 b_s_romG + 4 b_s_romA )
                # CIS: 4K ( 45 alpha2 + 7 b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 12*rinv.b_s_romG + 4*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        elif mode == 'SCP':
            if scat == '0':
                # CID: 8K/c ( 90 alp_romG + 2 b_s_romG - 2 b_s_romA )
                # CIS: 4K ( 45 alpha2 + 7 b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 90*rinv.alp_romG + 2*rinv.b_s_romG \
                        - 2*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph
            elif scat == '90x':
                # CID: 4K/c ( 45 alp_romG + 7 b_s_romG + b_s_romA )
                # CIS: 2K ( 45 alpha2 + 7 b_s_alph )
                if do_CID:
                    a_CID = 4.
                    x_CID = 45*rinv.alp_romG + 7*rinv.b_s_romG \
                        + rinv.b_s_romA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph
            elif scat == '90z':
                # CID: 8K/c ( 3 b_s_romG - b_s_romA )
                # CIS: 4K ( 3 b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 3*rinv.b_s_romG - rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 3*rinv.b_s_alph
            elif scat == '90u':
                # CID: 4K/c ( 45 alp_romG + 13 b_s_romG - b_s_romA )
                # CIS: 4K ( 45 alpha2 + 13 b_s_alph )
                if do_CID:
                    a_CID = 4.
                    x_CID = 45*rinv.alp_romG + 13*rinv.b_s_romG \
                        - rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 13*rinv.b_s_alph
            elif scat == '90*':
                # CID: 40K/3c ( 9 alp_romG + 2 b_s_romG )
                # CIS: 20K/3 ( 9 alpha2 + 2 b_s_alph )
                if do_CID:
                    a_CID = 40. / 3.
                    x_CID = 9*rinv.alp_romG + 2*rinv.b_s_romG
                if do_CIS:
                    a_CIS = 20. / 3.
                    x_CIS = 9*rinv.alpha2 + 2*rinv.b_s_alph
            elif scat[:3] == '180':
                # CID: 8K/c ( 12 b_s_romG + 4 b_s_romA )
                # CIS: 4K ( 45 alpha2 + 7 b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 12*rinv.b_s_romG + 4*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        elif mode == 'DCPI':
            if scat == '0':
                # CID: 8K/c ( 90 alp_romG + 2 b_s_romG - 2 b_s_romA )
                # CIS: 4K ( 45 alpha2 + b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 90*rinv.alp_romG + 2*rinv.b_s_romG \
                        - 2*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + rinv.b_s_alph
            elif scat[:2] == '90':
                # CID: 2K/c ( 90 alp_romG + 26 b_s_romG - 2 b_s_romA )
                # CIS: 2K ( 45 alpha2 + 13 b_s_alph )
                if do_CID:
                    a_CID = 2.
                    x_CID = 90*rinv.alp_romG + 26*rinv.b_s_romG \
                        - 2*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + 13*rinv.b_s_alph
            elif scat[:3] == '180':
                # CID: 16K/c ( 6 b_s_romG + 2 b_s_romA )
                # CIS: 24K ( b_s_alph )
                if do_CID:
                    a_CID = 16.
                    x_CID = 6*rinv.b_s_romG + 2*rinv.b_s_romA
                if do_CIS:
                    a_CIS = 24.
                    x_CIS = rinv.b_s_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        elif mode == 'DCPII':
            if scat == '0':
                # CID: 16K/c ( 6 b_s_romG + 2 b_s_romA )
                # CIS: 24K ( b_s_alph )
                if do_CID:
                    a_CID = 16.
                    x_CID = 0.
                if do_CIS:
                    a_CIS = 24.
                    x_CIS = rinv.b_s_alph
            elif scat[:2] == '90':
                # CID: 2K/c ( 90 alp_romG + 26 b_s_romG - 2 b_s_romA )
                # CIS: 2K ( 45 alpha2 + 13 b_s_alph )
                if do_CID:
                    a_CID = 2.
                    x_CID = 0.
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + 13*rinv.b_s_alph
            elif scat[:3] == '180':
                # CID: 8K/c ( 90 alp_romG + 2 b_s_romG - 2 b_s_romA )
                # CIS: 4K ( 45 alpha2 + b_s_alph )
                if do_CID:
                    a_CID = 8.
                    x_CID = 0.
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + rinv.b_s_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        else:
            raise IndexError('Unrecognized scattering setup.')
    else:
        if mode == 'ICP':
            if scat == '0':
                if do_CID:
                    a_CID = 8.
                    x_CID = 45*rinv.alp_romG + 7*rinv.b_s_romG \
                        + 5*rinv.b_a_romG + rinv.b_s_romA - rinv.b_a_romA \
                        - 45*rinv.alp_calG + 5*rinv.b_s_calG \
                        - 5*rinv.b_a_calG - 3*rinv.b_s_calA + rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph \
                        + 5*rinv.b_a_alph
            elif scat == '90x':
                if do_CID:
                    a_CID = 4.
                    x_CID = 45*rinv.alp_romG + 7*rinv.b_s_romG \
                        + 5*rinv.b_a_romG + rinv.b_s_romA - rinv.b_a_romA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph \
                        + 5*rinv.b_a_alph
            elif scat == '90z':
                if do_CID:
                    a_CID = 8.
                    x_CID = 3*rinv.b_s_romG + 5*rinv.b_a_romG \
                        - rinv.b_s_romA + rinv.b_a_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 3*rinv.b_s_alph + 5*rinv.b_a_alph
            elif scat == '90u':
                if do_CID:
                    a_CID = 4.
                    x_CID = 45*rinv.alp_romG + 13*rinv.b_s_romG \
                        + 15*rinv.b_a_romG - rinv.b_s_romA + rinv.b_a_romA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 13*rinv.b_s_alph \
                        + 15*rinv.b_a_alph
            elif scat == '90*':
                if do_CID:
                    a_CID = 40. / 3.
                    x_CID = 9*rinv.alp_romG + 2*rinv.b_s_romG \
                        + 2*rinv.b_a_romG
                if do_CIS:
                    a_CIS = 20. / 3.
                    x_CIS = 9*rinv.alpha2 + 2*rinv.b_s_alph \
                        + 2*rinv.b_a_alph
            elif scat[:3] == '180':
                if do_CID:
                    a_CID = 8.
                    x_CID = 45*rinv.alp_romG + 7*rinv.b_s_romG \
                        + 5*rinv.b_a_romG + rinv.b_s_romA - rinv.b_a_romA \
                        + 45*rinv.alp_calG - 5*rinv.b_s_calG \
                        + 5*rinv.b_a_calG + 3*rinv.b_s_calA - rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph \
                        + 5*rinv.b_a_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        elif mode == 'SCP':
            if scat == '0':
                if do_CID:
                    a_CID = 8.
                    x_CID = 45*rinv.alp_romG - 5*rinv.b_s_romG \
                        + 5*rinv.b_a_romG - 3*rinv.b_s_romA - rinv.b_a_romA \
                        - 45*rinv.alp_calG - 7*rinv.b_s_calG \
                        - 5*rinv.b_a_calG + rinv.b_s_calA + rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph \
                        + 5*rinv.b_a_alph
            elif scat == '90x':
                if do_CID:
                    a_CID = 4.
                    x_CID = -45*rinv.alp_calG - 7*rinv.b_s_calG \
                        - 5*rinv.b_s_calG + rinv.b_s_calA + rinv.b_a_calA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph \
                        + 5*rinv.b_a_alph
            elif scat == '90z':
                if do_CID:
                    a_CID = 8.
                    x_CID = -3*rinv.b_s_calG - 5*rinv.b_a_calG \
                        - rinv.b_s_calA - rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 3*rinv.b_s_alph + 5*rinv.b_a_alph
            elif scat == '90u':
                if do_CID:
                    a_CID = 4.
                    x_CID = -45*rinv.alp_calG - 13*rinv.b_s_calG \
                        - 15*rinv.b_a_calG - rinv.b_s_calA - rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 13*rinv.b_s_alph \
                        + 15*rinv.b_a_alph
            elif scat == '90*':
                if do_CID:
                    a_CID = 40. / 3.
                    x_CID = -9*rinv.alp_calG - 2*rinv.b_s_calG \
                        - 2*rinv.b_a_calG
                if do_CIS:
                    a_CIS = 20. / 3.
                    x_CIS = 9*rinv.alpha2 + 2*rinv.b_s_alph \
                        + 2*rinv.b_a_alph
            elif scat[:3] == '180':
                if do_CID:
                    a_CID = 8.
                    x_CID = -45*rinv.alp_romG + 5*rinv.b_s_romG \
                        - 5*rinv.b_a_romG + 3*rinv.b_s_romA + rinv.b_a_romA \
                        - 45*rinv.alp_calG - 7*rinv.b_s_calG \
                        - 5*rinv.b_a_calG + rinv.b_s_calA + rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + 7*rinv.b_s_alph \
                        + 5*rinv.b_a_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        elif mode == 'DCPI':
            if scat == '0':
                if do_CID:
                    a_CID = 2.
                    x_CID = 45*rinv.alp_romG + rinv.b_s_romG \
                        + 5*rinv.b_a_romG - rinv.b_s_romA - rinv.b_a_romA \
                        - 45*rinv.alp_calG - rinv.b_s_calG \
                        - 5*rinv.b_a_calG - rinv.b_s_calA + rinv.b_a_calA
                if do_CIS:
                    a_CIS = 4.
                    x_CIS = 45*rinv.alpha2 + rinv.b_s_alph + 5*rinv.b_a_alph
            elif scat[:2] == '90':
                if do_CID:
                    a_CID = 2.
                    x_CID = 45*rinv.alp_romG + 13*rinv.b_s_romG \
                        + 15*rinv.b_a_romG - rinv.b_s_romA + rinv.b_a_romA \
                        - 45*rinv.alp_calG - 13*rinv.b_s_calG \
                        - 15*rinv.b_a_calG - rinv.b_s_calA - rinv.b_a_calA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + rinv.b_s_alph + 5*rinv.b_a_alph
            elif scat[:3] == '180':
                if do_CID:
                    a_CID = 16.
                    x_CID = 3*rinv.b_s_romG + rinv.b_s_romA \
                        - 3*rinv.b_s_calG + rinv.b_s_calA
                if do_CIS:
                    a_CIS = 24.
                    x_CIS = rinv.b_s_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        elif mode == 'DCPII':
            if scat == '0':
                if do_CID:
                    a_CID = 16.
                    x_CID = 3*rinv.b_s_romG + 3*rinv.b_s_calG \
                        + rinv.b_s_romA - rinv.b_s_calA
                if do_CIS:
                    a_CIS = 24.
                    x_CIS = rinv.b_s_alph
            elif scat[:2] == '90':
                if do_CID:
                    a_CID = 2.
                    x_CID = 45*rinv.alp_romG + 13*rinv.b_s_romG \
                        + 15*rinv.b_a_romG - rinv.b_s_romA + rinv.b_a_romA \
                        + 45*rinv.alp_calG + 13*rinv.b_s_calG \
                        + 15*rinv.b_a_calG + rinv.b_s_calA + rinv.b_a_calA
                if do_CIS:
                    a_CIS = 2.
                    x_CIS = 45*rinv.alpha2 + rinv.b_s_alph + 5*rinv.b_a_alph
            elif scat[:3] == '180':
                if do_CID:
                    a_CID = 2.
                    x_CID = 45*rinv.alp_romG + rinv.b_s_romG \
                        + 5*rinv.b_a_romG - rinv.b_s_romA - rinv.b_a_romA \
                        + 45*rinv.alp_calG + rinv.b_s_calG \
                        + 5*rinv.b_a_calG + rinv.b_s_calA - rinv.b_a_calA
                if do_CIS:
                    a_CIS = 24.
                    x_CIS = rinv.b_s_alph
            else:
                raise IndexError('Unrecognized scattering setup.')
        else:
            raise IndexError('Unrecognized scattering setup.')

    res = None
    if do_ROA:
        if use_CID:
            if abs(x_CIS) > sys.float_info.epsilon:
                res = x_CID/x_CIS if get_activity else \
                    (a_CID*x_CID)/(a_CIS*x_CIS)
            else:
                res = 0.0
        else:
            res = x_CID if get_activity else a_CID*x_CID
    else:
        res = x_CIS if get_activity else a_CIS*x_CIS
    return res
