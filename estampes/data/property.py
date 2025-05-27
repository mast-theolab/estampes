"""Module on atoms-related data.

This module provides basic data related to atoms.
"""

from math import pi
import typing as tp

from estampes.base import TypeQOpt, TypeQTag
from estampes.data.physics import PHYSCNST, PHYSFACT


# ================
# Module Functions
# ================

def property_data(qtag: TypeQTag,
                  qopt: TypeQOpt = None
                  ) -> tp.NamedTuple:
    """Generate property data.

    Generates a named tuple containing basic data for a given
    property in argument given in input.

    Available information:

    name
        Full name of the property.
    dim
        Dimension (as integer, string or tuple).
    der
        Property is derivable (True, False).
    d1q
        Type of analytic 1st derivative wrt normal modes ('p', 'q').
    unit
        Atomic units.

    Parameters
    ----------
    qlabel
        Property label.
    qsublabel
        Property sub-label/option.

    Returns
    -------
    :obj:`tp.NamedTuple`
        Named tuple with the properties listed above.

    Raises
    ------
    KeyError
        Unrecognized atomic symbol.
    """
    item = str(qtag).lower()
    if item == '1':
        qname = 'energy'
        qdim = 1
        qder = True
        qd1q = 'q'
        qunit = 'Eh'
    elif item == '50':
        qname = 'non-adiabatic couplings'
        qdim = ('nat', 3)
        qder = True
        qd1q = 'p'
        qunit = '1/a0'
    elif item == '101':
        qname = 'electric dipole'
        qdim = 3
        qder = True
        qd1q = 'q'
        qunit = 'e.a0'
    elif item == '102':
        qname = 'magnetic dipole'
        qdim = 3
        qder = True
        qd1q = 'p'
        qunit = 'e.hbar/me'
    elif item == '103':
        qname = 'polarizability tensor'
        qdim = 6
        qder = True
        qd1q = 'q'
        qunit = 'a0^3'
    elif item == '104':
        qname = 'optical rotations'
        qdim = 9
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '105':
        qname = 'dipole-quadrupole polarizability'
        qdim = (3, 6)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '106':
        qname = 'hyperpolarizability'
        qdim = (3, 6)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '107':
        qname = 'quadrupole'
        qdim = 6
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '201':
        qname = 'magnetic susceptibility'
        qdim = (3, 3)
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '202':
        qname = 'rotational g-Tensor'
        qdim = (3, 3)
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '203':
        qname = 'NMR shielding tensors'
        qdim = ('nat', (3, 3))
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '204':
        qname = 'spin-rotation tensors'
        qdim = ('nat', (3, 3))
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '205':
        qname = 'anisotropic hyperfine tensors'
        qdim = ('nat', 6)
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '206':
        qname = 'isotropic (Fermi) terms'
        qdim = 'nat'
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '207':
        qname = 'ESR g-tensor'
        qdim = (3, 3)
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '208':
        qname = 'nuclear quadrupole tensors'
        qdim = ('nat', 6)
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '209':
        qname = 'isotropic spin-spin coupling'
        qdim = 'nattt'
        qder = False
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '301':
        qname = 'polarizability alpha(-w,w)'
        qdim = (3, 3)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '302':
        qname = 'optical rotations'
        qdim = (3, 3)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '303':
        qname = 'polarizability alpha(w,0)'
        qdim = (3, 3)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '304':
        qname = 'dipole-quadrupole polarizability'
        qdim = (3, 6)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '305':
        qname = 'hyperpolarizability beta(-w,w,0)'
        qdim = (3, 6)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    elif item == '306':
        qname = 'hyperpolarizability Beta(w,w,-2w)'
        qdim = (3, 6)
        qder = True
        qd1q = 'q'
        qunit = 'a.u.'
    else:
        raise IndexError('Unrecognized property')

    return QBaseInfo(name=qname, dim=qdim, der=qder, d1q=qd1q, unit=qunit)


def property_units(qtag: str, unit: str = 'SI') -> tp.Tuple[float, str]:
    """Return conversion factors and unit labels.

    Returns the conversion factor from atomic unit to new unit, as well
    as the label.

    Parameters
    ----------
    qtag
        Quantity label.
    unit
        Destination unit.
        SI: international unit system
        cgs: centimeter-gram-second unit

    Returns
    -------
    float
        Conversion unit from atomic unit.
    str
        Unit label.

    Raises
    ------
    ValueError
        Unsupported unit system.
    """
    _unit = unit.lower()
    if qtag == '101':
        if _unit == 'si':  # Conversion to C.m
            # conv = phys_fact('au2deb') * 100./PHYSCNST.slight * 1.0e-21
            conv = PHYSFACT.e2C * PHYSFACT.bohr2ang * 1.0e-10
            new_unit = 'C.m'
        elif _unit == 'cgs':  # Conversion to statC.cm
            conv = PHYSFACT.e2C*PHYSCNST.slight/10. * PHYSFACT.bohr2ang * 1.e-8
            new_unit = 'statC.cm'
    elif qtag == '102':
        # Note that we compute the nuclear magnetic dipole (not electronic)
        if _unit == 'si':  # Conversion to A.m^2 == J/T
            conv = PHYSCNST.planck/(2*pi) * PHYSFACT.e2C / PHYSFACT.amu2kg
            new_unit = 'A.m^2'
        elif _unit == 'cgs':  # Conversion to statA.cm^2
            conv = 1.0e4 * PHYSCNST.planck/(2*pi) \
                * (PHYSFACT.e2C*PHYSCNST.slight/10.) / PHYSFACT.amu2kg
            new_unit = 'statA.cm^2'
    return conv, new_unit


# ==============
# Module Classes
# ==============

class QBaseInfo(tp.NamedTuple):
    name: str  # Printable name
    dim: tp.Union[int, str, tp.Tuple[tp.Any]]  # dimension(s)
    der: bool  # Flag if property derivable
    d1q: str  # Type of analytical 1st derivative (p or q)
    unit: str  # Default unit
