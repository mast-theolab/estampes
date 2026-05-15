"""Provides basic types if NumPy is not available."""

from collections.abc import Sequence

TypeAtCrd = Sequence[Sequence[float]]
TypeAtLab = Sequence[str | int]
TypeAtMas = Sequence[float]
Type1Vib = Sequence[Sequence[float]]
TypeVibs = Sequence[Sequence[float]]

TypeAtLabM = TypeAtLab | Sequence[TypeAtLab]
TypeAtMasM = TypeAtMas | Sequence[TypeAtMas]
TypeAtCrdM = TypeAtCrd | Sequence[TypeAtCrd]
TypeVibsM = TypeVibs | Sequence[TypeVibs]
TypeVib1M = Type1Vib | Sequence[Type1Vib]
