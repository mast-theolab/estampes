"""Module providing basic classes and methods for ESTAMPES

Attributes
----------
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
TypeDFChk : dict
    Static type for data from Gaussian fchk file.
TypeDGLog : str, list, optional
    Static type for data from Gaussian log file.
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

Classes
-------
ArgumentError
    Generates an error for inconsistency/errors in arguments.
ParsingError
    Basic container for parsing-specific errors
ParseDataError
    Generates an error if quantity/property not available/supported
ParseKeyError
    Generates an error if keyword not found
QuantityError
    Generates an error if quantity is not supported
"""

# flake8: noqa: F401

from estampes.base.types import Type1Vib, TypeAtCrd, TypeAtCrdM, TypeAtData, \
    TypeAtLab, TypeAtLabM, TypeAtMas, TypeBonds, TypeBondsM, TypeColor, \
    TypeDCrd, TypeDFChk, TypeDGLog, TypeDOrd, TypeQLab, TypeQLvl, TypeQOpt, \
    TypeQTag, TypeRSta, TypeVibs, ConstDict

from estampes.base.errors import ArgumentError, DataError, ParseDataError, \
    ParseKeyError, ParsingError, QuantityError

from estampes.base.qlabel import QLabel, TypeQInfo
from estampes.base.qdata import QData, TypeQData
