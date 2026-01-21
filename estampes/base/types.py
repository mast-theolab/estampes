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
TypeVibsM : list, np.ndarray
    Static type for vibrational modes (multiple molecules).
"""
import typing as tp

try:
    import numpy.typing as npt
    has_np = True
except ImportError:
    has_np = False


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

# _tp_StrInt = tp.TypeVar('_tp_StrInt', str, int)
_tp_StrInt = tp.Union[str, int]

TypeRlCx = tp.TypeVar('TypeRlCx', float, complex)

# Label-related types
TypeQTag = _tp_StrInt
TypeQOpt = tp.Optional[_tp_StrInt]
TypeDOrd = tp.Optional[int]
TypeDCrd = tp.Optional[str]
TypeRSta = tp.Optional[tp.Union[str, int, tp.Tuple[_tp_StrInt, _tp_StrInt]]]
TypeQLvl = tp.Optional[str]
TypeQLab = tp.Tuple[TypeQTag, TypeQOpt, TypeDOrd, TypeDCrd, TypeRSta, TypeQLvl]

TypeDFChk = tp.Dict[str, tp.List[tp.Union[str, int, float]]]
TypeDGLog = tp.List[tp.Union[tp.List[str], str]]
TypeColor = tp.Union[tp.Sequence[int], tp.Sequence[float], float, str]

# Atoms-related data
TypeAtData = tp.Dict[_tp_StrInt, tp.Dict[str, tp.List[tp.Any]]]
if has_np:
    TypeAtCrd = npt.ArrayLike
    TypeAtLab = npt.ArrayLike
    TypeAtMas = npt.ArrayLike
    Type1Vib = npt.ArrayLike
    TypeVibs = npt.ArrayLike
else:
    TypeAtCrd = tp.Sequence[tp.Sequence[float]]
    TypeAtLab = tp.Sequence[_tp_StrInt]
    TypeAtMas = tp.Sequence[float]
    Type1Vib = tp.Sequence[tp.Sequence[float]]
    TypeVibs = tp.Sequence[tp.Sequence[float]]
TypeBonds = tp.List[tp.Tuple[int, int]]

# Basic types extended to support multiple molecules
TypeAtLabM = tp.Union[TypeAtLab, tp.Sequence[TypeAtLab]]
TypeAtCrdM = tp.Union[TypeAtCrd, tp.Sequence[TypeAtCrd]]
TypeBondsM = tp.Union[TypeBonds, tp.Sequence[TypeBonds]]
TypeVibsM = tp.Union[TypeVibs, tp.Sequence[TypeVibs]]
