"""Provide types specific to the glog module."""

import typing as tp
from collections.abc import Callable, Sequence

from estampes.base.types import StrIntType

QKwrdItemType = tuple[
    int,  # Link
    str,  # Keyword
    StrIntType,  # Jump/Skip function
    str,  # Matching pattern for data to extract
    Callable[[str], bool],  # end condition
    int  # Number of occurrences
]
QKwrdListType = tuple[
    Sequence[int],
    Sequence[str],
    Sequence[StrIntType],
    Sequence[str],
    Sequence[Callable[[str], bool]],
    Sequence[int]
]
QKwrdType = QKwrdItemType | QKwrdListType

KDataType = tuple[
    str,  # Keyword
    int,  # Link
    StrIntType,  # Information on lines to skip after keyword
    tp.Pattern,  # Data extraction matching pattern (compiled)
    int,  # which occurrences to extract
    Callable[[str], bool]
]
