"""Provide basic types classes.

A basic module providing types specifications and new types for ESTAMPES.

Attributes
----------
TypeRlCx : float, complex
    Type variable designated either float or complex.
Type1Vib : list, np.ndarray
    Static type for 1 vibrational mode, expected in form (NAt, 3).
TypeAtCrd : list, np.ndarray
    Static type for atomic coordinates.
TypeAtCrdM
    Static type for atomic coordinates (multiple molecules).
TypeAtData : dict
    Static type for atom data.
TypeAtLab : list
    Static type for atomic labels.
TypeAtLabM
    Static type for atomic labels (multiple molecules).
TypeAtMas : list, np.ndarray
    Static type for atomic masses.
TypeBonds : list
    Static type for bond list, as (atom1, atom2).
TypeBondsM
    Static type for bonds information (multiple molecules).
TypeColor : float, str, list
    Static type for colors.
TypeDCrd : str, optional
    Static type for derivative coordinate.
TypeDBlocGLog : str, list, optional
    Static type for data blocks extracted from Gaussian log file.
TypeDFChk : dict
    Static type for data from Gaussian fchk file.
TypeDOrd : int, optional
    Static type for derivative order (0: property).
TypeQData : dict
    Static type for data returned by parsers.
TypeQInfo : dict
    Static type for dictionary of quantity full labels.
TypeQLab : tuple
    Static type for quantity label.
TypeQLvl : str, optional
    Static type for level of theory used to compute quantity.
TypeQOpt : str, int, optional
    Static type for quantity option.
TypeQTag : str, int
    Static type for quantity tag.
TypeRSta : str, int, tuple, optional
    Static type for reference state/transition.
TypeVibs : list, np.ndarray
    Static type for vibrational modes, expected in form (Nib, NAt3).
TypeVibsM : list, np.ndarray
    Static type for vibrational modes (multiple molecules).
"""
import typing as tp
from collections.abc import Sequence


# ==============
# Module Classes
# ==============

class ConstDict(dict):
    """Derived type from dict offering attribute style access.

    A type derived from `dict`, which offers the possibility to access
    keys as attributes.

    References
    ----------
    https://goodcode.io/attributes/python-dict-object/
    """

    def __getattr__(self, name):
        """Get attribute 'name'."""
        if name in self:
            return self[name]
        else:
            raise AttributeError(f'No such attribute: {name}')

    def __setattr__(self, name, value):
        """Set attribute 'name' to 'value'."""
        self[name] = value

    def __delattr__(self, name):
        """Delete attribute 'name'."""
        if name in self:
            del self[name]
        else:
            raise AttributeError(f'No such attribute: {name}')


# =================
# Module Attributes
# =================

TypeStrInt = str | int
TypeRlCx = float | complex

# Label-related types
TypeQTag = TypeStrInt
TypeQOpt = TypeStrInt | None
TypeDOrd = int
TypeDCrd = str | None
TypeRSta = str | int | tuple[TypeStrInt, TypeStrInt] | None
TypeQLvl = str | None
TypeQLab = tuple[TypeQTag, TypeQOpt, TypeDOrd, TypeDCrd, TypeRSta,
                 TypeQLvl]

TypeDFChk = dict[str, list[str] | list[int] | list[float]]
TypeDBlocGLog = Sequence[Sequence[str] | str]
TypeColor = Sequence[int] | Sequence[float] | float | str

# Atoms-related data
TypeAtData = dict[TypeStrInt, dict[str, list[tp.Any]]]
TypeBonds = list[tuple[int, int]]
TypeBondsM = TypeBonds | Sequence[TypeBonds]

# pylint: disable=W0611
# flake8: noqa: F401

try:
    import numpy.typing
    from estampes.base.types_npt import (
        TypeAtCrd, TypeAtLab, TypeAtMas, Type1Vib, TypeVibs,
        TypeAtLabM, TypeAtMasM, TypeAtCrdM, TypeVibsM, 
        TypeVib1M)
except ImportError:
    from estampes.base.types_base import (
        TypeAtCrd, TypeAtLab, TypeAtMas, Type1Vib, TypeVibs,
        TypeAtLabM, TypeAtMasM, TypeAtCrdM, TypeVibsM, TypeVib1M)
