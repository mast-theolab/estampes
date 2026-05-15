"""Provides special types when NumPy is available."""

from collections.abc import Sequence

from numpy.typing import ArrayLike

TypeAtCrd = ArrayLike
TypeAtLab = ArrayLike
TypeAtMas = ArrayLike
Type1Vib = ArrayLike
TypeVibs = ArrayLike
TypeAtLabM = TypeAtLab | Sequence[TypeAtLab]
TypeAtMasM = TypeAtMas | Sequence[TypeAtMas]
TypeAtCrdM = TypeAtCrd | Sequence[TypeAtCrd]
TypeVibsM = TypeVibs | Sequence[TypeVibs]
TypeVib1M = Type1Vib | Sequence[Type1Vib]
