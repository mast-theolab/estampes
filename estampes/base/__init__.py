"""Module providing basic classes and methods for ESTAMPES

Attributes
----------
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

from estampes.base.types import TypeAtCrd, TypeAtCrdM, TypeAtData, TypeAtLab, \
    TypeAtLabM, TypeAtMas, TypeBonds, TypeBondsM, TypeColor, TypeQData, \
    TypeDCrd, TypeDFChk, TypeDGLog, TypeDOrd, TypeQInfo, TypeQLab, TypeQLvl, \
    TypeQOpt, TypeQTag, TypeRSta, ConstDict

from estampes.base.errors import ArgumentError, ParseDataError, \
    ParseKeyError, ParsingError, QuantityError
