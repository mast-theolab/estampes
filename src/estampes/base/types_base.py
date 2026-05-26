"""Provides basic types if NumPy is not available."""

from collections.abc import Sequence

AtCrdType = Sequence[float]
AtsCrdType = Sequence[AtCrdType]

AtMasType = float
AtsMasType = Sequence[AtMasType]

VibType = Sequence[float] | Sequence[Sequence[float]]
VibsType = Sequence[Sequence[float]] | Sequence[Sequence[Sequence[float]]]
