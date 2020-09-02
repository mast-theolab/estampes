"""Module for the file parsing in ESTAMPES.

This module provides high-levels classes and methods, independent of
  the format:
DataFile
    Main class to manipulate data.
parse_qlabel
    Parses a full quantity label into sub-items.
build_qlabel
    Complement `parse_qlabel`, doing the opposite.
reshape_dblock
    Reshapes data block.

See submodules for details.

Notes
-----
* While the parser is format dependent, the availability of data will
  obviously depend on the file and the way it was generated.
"""

from estampes.parser.base import parse_qlabel, build_qlabel, reshape_dblock, DataFile  # NOQA
