"""Provide types specific to the glog module."""

import typing as tp

_tp_StrInt = tp.TypeVar('_tp_StrInt', str, int)

TypeQKwrd = tp.Tuple[
    tp.Union[int, tp.List[int]],  # Link
    tp.Union[str, tp.List[str]],  # Keyword
    tp.Union[_tp_StrInt, tp.List[_tp_StrInt]],  # Jump/Skip function
    tp.Union[str, tp.List[str]],  # Matching pattern for data to extract
    #  Block end condition
    tp.Union[tp.Callable[[str], bool], tp.List[tp.Callable[[str], bool]]],
    tp.Union[int, tp.List[int]]  # Number of occurrences
]

TypeKData = tp.Tuple[
    str,  # Keyword
    int,  # Link
    _tp_StrInt,  # Information on lines to skip after keyword
    tp.Pattern,  # Data extraction matching pattern (compiled)
    int,  # which occurrences to extract
    tp.Callable[[str], bool]
]
