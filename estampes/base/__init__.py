"""Module providing basic classes and methods for ESTAMPES

A module providing basic submodules to be used by ESTAMPES:
atom
    Atom-related data.
physics
    Basic physical constants and methods.
property
    Molecular property-related data.
molecule
    Provides the `Molecule` class.
state
    Provides the `ElectronicState` and `MolecularProperties` classes.
types
    Provides special static types used in ESTAMPES.

Attributes
----------
TypeQTag : :obj:`typing.TypeVar`
    Static type for quantity tag.
TypeQOpt : :obj:`typing.Optional`
    Static type for quantity option.
TypeDOrd : :obj:`typing.Optional`
    Static type for derivative order (0: property).
TypeDCrd : :obj:`typing.Optional`
    Static type for derivative coordinate.
TypeRSta : :obj:`typing.Optional`
    Static type for reference state/transition.
TypeQLab : :obj:`typing.Tuple`
    Static type for quantity label.
TypeData : :obj:`typing.Dict`
    Static type for data returned by parsers.
TypeAtData : :obj:`typing.Dict`
    Static type for atom data.
TypeQInfo : :obj:`typing.Dict`
    Static type for dictionary of quantity full labels.
TypeDFChk : :obj:`typing.Dict`
    Static type for data from Gaussian fchk file.
TypeDGLog : :obj:`typing.List`
    Static type for data from Gaussian log file.
TypeColor : :obj:`typing.Union`
    Static type for colors.
TypeAtLab : :obj:`typing.Sequence`
    Static type for atomic labels.
TypeBonds : :obj:`typing.List`
    Static type for bond list, as (atom1, atom2).

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

from estampes.base.types import TypeAtData, TypeAtLab, TypeBonds, TypeColor, \
    TypeData, TypeDCrd, TypeDFChk, TypeDGLog, TypeDOrd, TypeQInfo, TypeQLab, \
    TypeQLvl, TypeQOpt, TypeQTag, TypeRSta, ConstDict

from estampes.base.errors import ArgumentError, ParseDataError, \
    ParseKeyError, ParsingError, QuantityError
