"""Module for the file parsing in ESTAMPES.

This module provides high-levels classes and methods, independent of
the format:

Sub-modules
-----------
DataFile
    Main class to manipulate data.
parse_qlabels
    Parses a list of qlabels as QLabels objects or strings.
reshape_dblock
    Reshapes data block.

See submodules for details.

Notes
-----
* While the parser is format dependent, the availability of data will
  obviously depend on the file and the way it was generated.
"""

from estampes.parser.base import DataFile  # noqa: F401
from estampes.parser.functions import reshape_dblock, parse_qlabels
