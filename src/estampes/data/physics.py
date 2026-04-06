"""Module on physical constants and related data.

This module provides basic data related to physical constants and
conversion factor.

Attributes
----------
PHYSFACT : dict
    Physical conversion factors.
PHYSCNST : dict
    Physical constants.
"""

from math import pi, sqrt

from estampes.base import ConstDict


# ================
# Module Constants
# ================

PHYSFACT = ConstDict([
    ('bohr2ang', 0.52917720859),  # Conversion bohr -> Angstrom
    ('amu2kg', 1.660538782e-27),  # Conversion atomic mass unit -> kilogram
    ('e2C', 1.602176487e-19),  # Conversion electron charge -> Coulomb
    ('cal2J', 4.184),  # Conversion calorie -> Joule
    ('Eh2J', 4.35974394e-18),  # Conversion Hartree -> Joule
])

PHYSCNST = ConstDict([
    ('planck', 6.62606896e-34),  # Constant: Planck (J.s)
    ('avogadro', 6.02214179e23),  # Constant: Avogadro number (mol^-1)
    ('slight', 2.99792458e10),  # Constant: Speed of light (cm/s)
    ('boltzmann', 1.3806504e-23),  # Constant: Boltzmann (J/K)
    ('finestruct', 1./137.035999679),  # Constant: Fine structure (no unit)
    ('molvol', 22.413996e-3),  # Molar volume of ideal gas (m^3@273.15K)
    ('emagmom', -928.476377e-26),  # Electron Magnetic Moment (J/Tesla)
    ('prestmass', 1.672621637e-27),  # Proton rest mass (kg)
    ('gfactor', 2.0023193043622),  # Free electron g-factor (no unit)
])


# ==============
# Module Methods
# ==============

def phys_fact(factor: str) -> float:
    """Physical conversion factor.

    Returns some common conversion factors, especially for quantities
    often used in spectroscopy.

    Parameters
    ----------
    factor
        String identifier of the requested factor.

    Returns
    -------
    float
        Conversion factor.

    Note
    ----
    Convention:

    :math:`A`
        Angstroms
    :math:`a_0`
        atomic unit of length (Bohr)
    :math:`\\text{aJ}`
        attojoule
    :math:`E_h`
        atomic unit of energy (Hartree)
    :math:`k_B`
        Boltzmann constant
    :math:`N_a`
        Avogadro constant
    :math:`u`
        unified atomic mass
    :math:`abc`
        alias for :math:`u^{1/2}.a_0.cm^{-1/2}`

    Available factors:

    ========  =========================================================
     Label     Unit
    ========  =========================================================
     au2amu    a.u. (:math:`m_e`) to a.m.u ("proton" mass)
     au2cm1    a.u. (:math:`E_h`) to cm^-1
     au2Deb    a.u. (:math:`e^-.a_0`) to Debye conversion factor
     au2ESU    a.u. (:math:`e^-`) electron to electric static unit
     au2kg     a.u. (:math:`m_e`) to kg
     hbar      :math:`h/(2.pi)`
               in :math:`u.A^2.s^{-1}`
     Fac0AU    :math:`1/(h.c)`
               in :math:`\\text{cm}^{-1}.{E_h}^{-1}`
     Fac1AU    :math:`1/(2.\\pi.h^{1/2}.c^{3/2})`
               in :math:`\\text{cm}^{-1}.abc.{E_h}^{-1}`
     Fac2AU    :math:`1/(4.\\pi^2.c^2)`
               in :math:`\\text{cm}^{-1}.(abc)^2.{E_h}^{-1}`
     Fac3AU    :math:`1/(8.\\pi^3.h^{-1/2}.c^{5/2})`
               in :math:`\\text{cm}^{-1}.(abc)^3.{E_h}^{-1}`
     Fac4AU    :math:`1/(16.\\pi^4.h^{-1}.c^3)`
               in :math:`\\text{cm}^{-1}.(abc)^4.{E_h}^{-1}`
     Fact1     :math:`1/(2.\\pi.h^{1/2}.c^{3/2})`
               in :math:`u^{1/2}.A.\\text{cm}^{-3/2}.\\text{aJ}^{-1}`
     Fact2     :math:`1/(4.\\pi^2.c^2)`
               in :math:`u.A^2.\\text{cm}^{-2}.\\text{aJ}^{-1}`
     Fact3     :math:`1/(8.\\pi^3.h^{-1/2}.c^{5/2})`
               in :math:`u^{3/2}.A^3.\\text{cm}^{-5/2}.\\text{aJ}^{-1}`
     Fact4     :math:`1/(16.\\pi^4.h^{-1}.c^3)`
               in :math:`u^2.A^4.\\text{cm}^{-3}.\\text{aJ}^{-1}`
     FactA     :math:`h/(8.\\pi^2.c^2)`
               in :math:`A^2.u.\\text{cm}^{-1}`
     FactB     :math:`h.c/(2.k_B)`
               in cm.K
     FactC     :math:`k_B/c^2`
               in :math:`u.{a_0}^2.\\text{cm}^{-2}.K^{-1}`
     FactG     :math:`4.\\pi^2.c/h`
               in :math:`\\text{cm}.u^{-1}.A^{-2}`
     HC        :math:`h.c`
               in attoJ.cm == mdyn.Ang.cm.
     MWQ2q     mass-weighted to dimensionless normal coordinates,
               :math:`h^{1/2} / (2.\\pi.c^{1/2})`
               in :math:`a_0.u^{1/2}.\\text{cm}^{-1/2}`
     PiCH12    :math:`\\pi.(c/h)^{1/2}`
               in :math:`\\text{cm}^{1/2}.A^{-1}.u^{-1/2}`
     ToUMA     amu to uma conversion factor
    ========  =========================================================
    """

    m2ang = 1.0e10  # meter-to-Angstroms conversion factor
    # hbar: Reduced Planck constant in amu.Ang^2.s^-1
    hbar = PHYSCNST.planck*m2ang**2/(2.*pi*PHYSFACT.amu2kg)
    # fact_g: 4.pi^2.c/h in in cm.amu^-1.Ang^-2
    fact_g = 2.*pi*PHYSCNST.slight/hbar
    # hc, in attoJ.cm == mdyn.Ang.cm
    hxc = PHYSCNST.planck*PHYSCNST.slight*1.0e18

    au2ESU = PHYSFACT.e2C*PHYSCNST.slight/10.
    au2kg = 1.0e4*PHYSFACT.Eh2J \
        / (PHYSCNST.slight*PHYSCNST.finestruct)**2
    key = factor.strip().lower()

    if key == 'fac0au':
        res = PHYSFACT.Eh2J*1.0e18/hxc
    elif key in ('fac1au', 'fact1'):
        res = 1.0/(sqrt(fact_g)*hxc)
        if key == 'fac1au':
            res = res*PHYSFACT.Eh2J*1.0e18/PHYSFACT.bohr2ang
    elif key in ('fac2au', 'fact2'):
        res = 1.0/(fact_g*hxc)
        if key == 'fac2au':
            res = res*PHYSFACT.Eh2J*1.0e18/PHYSFACT.bohr2ang**2
    elif key in ('fac3au', 'fact3'):
        res = 1.0/(fact_g*sqrt(fact_g)*hxc)
        if key == 'fac3au':
            res = res*PHYSFACT.Eh2J*1.0e18/PHYSFACT.bohr2ang**3
    elif key in ('fac4au', 'fact4'):
        res = 1.0/(fact_g**2*hxc)
        if key == 'fac4au':
            res = res*PHYSFACT.Eh2J*1.0e18/PHYSFACT.bohr2ang**4
    elif key == 'factg':
        res = fact_g
    elif key == 'facta':
        res = 1.0/(2.0*fact_g)
    elif key == 'factb':
        res = PHYSCNST.planck*PHYSCNST.slight/(2.0*PHYSCNST.boltzmann)
    elif key == 'factc':
        res = PHYSCNST.boltzmann*m2ang**2\
                / (PHYSCNST.slight**2*PHYSFACT.amu2kg*PHYSFACT.bohr2ang**2)
    elif key == 'hbar':
        res = hbar
    elif key == 'touma':
        res = PHYSFACT.amu2kg/au2kg
    elif key == 'pich12':
        res = pi*sqrt(PHYSCNST.slight*PHYSFACT.amu2kg/PHYSCNST.planck/m2ang**2)
    elif key == 'au2deb':
        res = 1.0e10*au2ESU*PHYSFACT.bohr2ang
    elif key == 'mwq2q':
        res = 1.0/(sqrt(fact_g)*PHYSFACT.bohr2ang)
    elif key == 'au2cm1':
        res = PHYSFACT.Eh2J/(PHYSCNST.planck*PHYSCNST.slight)
    elif key == 'au2amu':
        res = au2kg/PHYSFACT.amu2kg
    elif key == 'au2esu':
        res = au2ESU
    elif key == 'au2kg':
        res = au2kg
    else:
        res = False
    return res
