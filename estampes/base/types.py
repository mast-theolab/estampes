"""Module providing basic types classes

A basic module providing types specifications and new types for ESTAMPES.

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

Classes
-------
ConstDict
    Derived class from `dict` offering attribute form.
"""

import typing as tp


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
        if name in self:
            return self[name]
        else:
            raise AttributeError(f'No such attribute: {name}')

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError(f'No such attribute: {name}')


# =================
# Module Attributes
# =================

_tp_StrInt = tp.TypeVar('_tp_StrInt', str, int)

TypeQTag = _tp_StrInt
TypeQOpt = tp.Optional[_tp_StrInt]
TypeDOrd = tp.Optional[int]
TypeDCrd = tp.Optional[str]
TypeRSta = tp.Optional[tp.Union[str, int, tp.Tuple[_tp_StrInt, _tp_StrInt]]]
TypeQLab = tp.Tuple[TypeQTag, TypeQOpt, TypeDOrd, TypeDCrd, TypeRSta]
TypeData = tp.Dict[str, tp.Dict[str, tp.Any]]
TypeAtData = tp.Dict[str, tp.Dict[str, tp.List[tp.Any]]]
TypeQInfo = tp.Dict[str, tp.List[tp.Any]]
TypeDFChk = tp.Dict[str, tp.List[tp.Union[str, int, float]]]
TypeDGLog = tp.List[tp.Union[tp.List[str], str]]
