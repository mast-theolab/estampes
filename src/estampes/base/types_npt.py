"""Provides special types when NumPy is available."""

from collections.abc import Sequence

from numpy.typing import NDArray

AtCrdType = Sequence[float] | NDArray
AtsCrdType = Sequence[AtCrdType] | NDArray

AtMasType = float
AtsMasType = Sequence[AtMasType] | NDArray

VibType = Sequence[float] | Sequence[Sequence[float]] | NDArray
VibsType = Sequence[Sequence[float]] | Sequence[Sequence[Sequence[float]]] \
    | NDArray
