"""Provide basic types classes.

A basic module providing types specifications and new types for ESTAMPES.

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
RealCplxType
    Type variable designated either float or complex.
StrIntType
    Type for string or integer format.
VibType
    Static type for 1 vibrational mode, expected in form (NAt, 3).
VibsType
    Static type for vibrational modes, expected in form (Nib, NAt3).
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



# Atoms-related data
class TypeAtDat(tp.TypedDict):
    """Type for Atomic Data structure."""
    symb: str
    name: str
    num: int
    mass: float
    rcov: tuple[int | None, ...]
    rvdw: float | None
    rvis: float | None
    rgb: tuple[int, ...] | None

# =================
# Module Attributes
# =================

StrIntType = str | int
RealCplxType = float | complex

# Label-related types
QLabTagType = StrIntType
QLabSubType = StrIntType | None
QLabDerType = int
QLabCrdType = str | None
QLabStaType = str | int | tuple[StrIntType, StrIntType] | None
QLabLvlType = str | None
QLabelType = tuple[QLabTagType, QLabSubType, QLabDerType, QLabCrdType,
                   QLabStaType, QLabLvlType]

DBlocFChkType = dict[str, list[str] | list[int] | list[float]]
DBlocGLogType = Sequence[Sequence[str] | str]
ColorType = Sequence[int] | Sequence[float] | float | str

AtLabType = str | int
AtDatType = dict[AtLabType, TypeAtDat]
AtsLabType = Sequence[AtLabType]
MAtsLabType = AtsLabType | Sequence[AtsLabType]
BondType = tuple[int, int]
BondsType = Sequence[BondType]
MBondsType = BondsType | Sequence[BondsType]


# pylint: disable=W0611
# flake8: noqa: F401

try:
    from estampes.base.types_npt import (
        AtCrdType, AtsCrdType, AtMasType, AtsMasType, VibType, VibsType)
except (ImportError, ModuleNotFoundError):
    from estampes.base.types_base import (
        AtCrdType, AtsCrdType, AtMasType, AtsMasType, VibType, VibsType)

MAtsMasType = AtsMasType | Sequence[AtsMasType]
MAtsCrdType = AtsCrdType | Sequence[AtsCrdType]
MVibsType = VibsType | Sequence[VibsType]
MVibType = VibType | Sequence[VibType]
