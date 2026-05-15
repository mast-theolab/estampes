"""Provide types specific to the glog module."""

import typing as tp
from collections.abc import Callable, Sequence

from estampes.base.types import TypeStrInt

TypeQKwrdItem = tuple[
    int,  # Link
    str,  # Keyword
    TypeStrInt,  # Jump/Skip function
    str,  # Matching pattern for data to extract
    Callable[[str], bool],  # end condition
    int  # Number of occurrences
]
TypeQKwrdList = tuple[
    Sequence[int],
    Sequence[str],
    Sequence[TypeStrInt],
    Sequence[str],
    Sequence[Callable[[str], bool]],
    Sequence[int]
]
TypeQKwrd = TypeQKwrdItem | TypeQKwrdList

TypeKData = tuple[
    str,  # Keyword
    int,  # Link
    TypeStrInt,  # Information on lines to skip after keyword
    tp.Pattern,  # Data extraction matching pattern (compiled)
    int,  # which occurrences to extract
    Callable[[str], bool]
]
