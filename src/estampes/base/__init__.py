"""Module providing basic classes and methods for ESTAMPES

Attributes
----------
AtCrdType
    Static type for 1 atomic coordinate.
AtDatType
    Static type for atom data.
AtLabType
    Static type for 1 atomic label.
AtMasType
    Static type for a single atomic mass.
AtsCrdType
    Static type for a set of atomic coordinates.
AtsLabType
    Static type for atomic labels.
AtsMasType
    Static type for a set of atomic mass.
BondType
    Static type for a single bond, as (atom1, atom2).
BondsType
    Static type for a list of bonds.
ColorType
    Static type for colors.
DBlocGLogType
    Static type for data blocks extracted from Gaussian log file.
DBlocFChkType
    Static type for data from Gaussian fchk file.
MAtsCrdType
    Static type for multiple sets of atomic coordinates.
MAtsLabType
    Static type for multiple sets of atomic labels.
MAtsMasType
    Static type for multiple sets of atomic masses.
MBondsType
    Static type for multiple sets of bonds lists.
MVibType
    Static type for multiple sets of a single vibrations.
MVibsType
    Static type for multiple sets of vibrations.
QDataBaseType : dict
    Static type for extracted data.
QInfoType : dict
    Static type for dictionary of full quantity labels.
QLabCrdType
    Static type for derivative coordinate.
QLabDerType
    Static type for derivative order (0: property).
QLabelType
    Static type for quantity label.
QLabLvlType
    Static type for level of theory used to compute quantity.
QLabStaType
    Static type for reference state/transition.
QLabSubType
    Static type for quantity option.
QLabTagType
    Static type for quantity tag.
QParseDataType : dict
    Static type for data returned by parsers.
RealCplxType
    Type variable designated either float or complex.
StrIntType
    Type for string or integer format.
VibType
    Static type for 1 vibrational mode, expected in form (NAt, 3).
VibsType
    Static type for vibrational modes, expected in form (Nib, NAt3).

Classes
-------
ArgumentError
    Generates an error for inconsistency/errors in arguments.
DataError
    Generic error related to data.
InternalError
    Generic error related to some internal inconsistency.
ParsingError
    Basic container for parsing-specific errors.
ParseDataError
    Generates an error if quantity/property not available/supported.
ParseKeyError
    Generates an error if keyword not found.
QuantityError
    Generates an error if quantity is not supported.
"""

# flake8: noqa: F401

from estampes.base.types import (
    AtCrdType, AtDatType, AtLabType, AtMasType,
    AtsCrdType, AtsLabType, AtsMasType,
    BondType, BondsType, ColorType,
    ConstDict, DBlocFChkType, DBlocGLogType,
    MAtsCrdType, MAtsLabType, MAtsMasType,
    MBondsType, MVibType, MVibsType,
    QLabCrdType, QLabDerType, QLabelType, QLabLvlType,
    QLabStaType, QLabSubType, QLabTagType, VibsType, VibType)

from estampes.base.errors import (
    ArgumentError, DataError, InternalError, ParseDataError, ParseKeyError,
    ParsingError, QuantityError)

from estampes.base.qlabel import QLabel, QInfoType
from estampes.base.qdata import QData, QDataBaseType, QParseDataType
