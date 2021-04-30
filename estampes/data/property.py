"""Module on atoms-related data.

This module provides basic data related to atoms.

Attributes
----------

Methods
-------
property_data
    Gets property data.
"""

import typing as tp

from estampes.base import TypeQOpt, TypeQTag


# ================
# Module Functions
# ================

def property_data(qlabel: TypeQTag,
                  qsublabel: TypeQOpt = None
                  ) -> tp.NamedTuple:
    """Generates property data.

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
    item = str(qlabel).lower()
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


# ==============
# Module Classes
# ==============

class QBaseInfo(tp.NamedTuple):
    name: str  # Printable name
    dim: tp.Union[int, str, tp.Tuple[tp.Any]]  # dimension(s)
    der: bool  # Flag if property derivable
    d1q: str  # Type of analytical 1st derivative (p or q)
    unit: str  # Default unit
