"""Basic module to parse XYZ files.

This module provides basic classes and methods to parse XYZ files.

Methods
-------
get_data
    Gets data from a XYZ file for each quantity label.

Classes
-------
FileXYZ
    Main class to handle XYZ files.
"""

import os
import typing as tp

from estampes import parser as ep
from estampes.base import ParseDataError, ParseKeyError, TypeData
from estampes.tools.atom import convert_labsymb
from estampes.data.physics import PHYSFACT


# ================
# Module Constants
# ================

__ang2au = 1.0 / PHYSFACT.bohr2ang


# ==============
# Module Classes
# ==============

class FileXYZ(object):
    """Main class to handle XYZ files.

    Attributes
    ----------
    filename : str
        File name.
    full_version : tuple
        full version:
        * XYZ
        * None

    Methods
    -------
    read_data(to_find, geom, raise_error)
        Extracts 1 or more data blocks from the XYZ file.
    """
    def __init__(self, fname: str) -> None:
        self.filename = fname
        self.analyze_file()

    @property
    def filename(self) -> str:
        """Gets or sets the filename associated to the XYZ object."""
        return self.__fname

    @filename.setter
    def filename(self, name: str) -> None:
        if not os.path.exists(name):
            raise FileNotFoundError('XYZ file not found')
        self.__fname = name

    @property
    def full_version(self) -> tp.Tuple[str, tp.Any]:
        """Returns the full version, for the parser interface"""
        return "XYZ", None

    @property
    def natoms(self) -> int:
        """Returns the number of atoms"""
        return self.__natoms

    @property
    def nstruct(self) -> int:
        """Returns the number of structures stored in the XYZ file"""
        return self.__ngeoms

    def analyze_file(self):
        """Analyzes XYZ file to get basic structural information.

        Analyzes a XYZ file and reports the number of geometries stored
          and the number of atoms to easily parse the file.
        """
        with open(self.__fname, 'r') as fobj:
            self.__ngeoms = 1
            line = fobj.readline()
            try:
                self.__natoms = int(line.strip())
            except ValueError:
                msg = '1st line should contain the number of atoms'
                raise ParseDataError(msg)
            line = fobj.readline()
            self.__title = line.strip()
            for _ in range(self.__natoms):
                line = fobj.readline()
                if not line:
                    msg = 'File ended while reading molecular geometry'
                    raise ParseDataError(msg)
            to_find = ['natoms']
            while True:
                try:
                    datlist = parse_xyz(fobj, to_find, False)
                except ParseDataError:
                    fmt = 'Error in geometry num. {}.\n' \
                        + 'Empty lines found between configurations.'
                    raise ParseDataError(fmt.self.__ngeoms+1) from None
                except ValueError as e:
                    fmt = 'Error in geometry num. {}.\n{}'
                    raise ValueError(fmt.format(self.__ngeoms+1, e)) from None
                except EOFError:
                    fmt = 'File ended while parsing geometry {}'
                    raise ParseDataError(fmt.format(self.__ngeoms+1))
                if datlist is None:
                    break
                else:
                    self.__ngeoms += 1
                    if datlist['natoms'] != self.__natoms:
                        fmt = 'Number of atoms in geometry {} does not match.'
                        raise ValueError(fmt.format(self.__ngeoms))

    def read_data(self, *to_find: tp.Sequence[str],
                  geom: tp.Optional[int] = 1,
                  raise_error: bool = True) -> tp.Dict[str, tp.Any]:
        """Extracts data corresponding to the keys to find.

        Parameters
        ----------
        to_find
            List of keys to find.
        geom
            Geometry block of interest, starting from 1
            -1: last
            0: all (only for coordinates)
        raise_error : bool
            Only raises error if `True`, otherwise proceeds silently.
        """
        ls_keys = set(to_find)
        datlist = {}
        if geom > self.__ngeoms:
            raise ValueError('Structure {} not available'.format(geom))
        if geom == -1:
            geom_ = self.__ngeoms
        else:
            geom_ = geom
        if 'natoms' in to_find:
            datlist['natoms'] = self.__natoms
        rest = ls_keys - {'natoms', 'title', 'atoms', 'atcrd'}
        if len(rest) > 0 and raise_error:
            raise ParseKeyError(rest[0])
        if {'title', 'atoms', 'atcrd'} & ls_keys:
            if geom_ > 0 or 'atcrd' not in ls_keys:
                with open(self.__fname, 'r') as fobj:
                    for _ in range((geom_-1)*(self.__natoms+2)):
                        _ = fobj.readline()
                    datlist = parse_xyz(fobj, to_find)
            else:
                with open(self.__fname, 'r') as fobj:
                    datlist = parse_xyz(fobj, to_find)
                    datlist['atcrd'] = [datlist['atcrd']]
                    i = 1
                    while i < self.__ngeoms:
                        i += 1
                        datlist['atcrd'].append(
                            parse_xyz(fobj, ['atcrd'])['atcrd'])
        return datlist

# ================
# Module Functions
# ================


def parse_xyz(fobj: tp.IO,
              what: tp.Sequence[str],
              blank_ok: bool = True
              ) -> tp.Optional[tp.Dict[str, tp.Any]]:
    """Parses a XYZ configuration.

    Parses an xyz configuration, assuming the first line to read is the
      number of atoms.
    On return, the pointer to `fobj` is on the last line of the
      configuration.

    Parameters
    ----------
    fobj
        File object: file already opened.
    what
        Data to extract: natoms, title, atcrd, atoms.
    blank_ok
        If true, blank lines before a configuration are OK.

    Returns
    -------
    dict
        Extracted data

    Raises
    ------
    EOFError
        End-of-file reached while reading the configuration.
    ParseDataError
        Error in format of XYZ (e.g., leading blank lines)
    ValueError
        Unexpected type of data.
    """
    # Initialize
    data = {}
    if 'natoms' in what:
        data['natoms'] = 0
    if 'title' in what:
        data['title'] = ''
    if 'atoms' in what:
        data['atoms'] = []
    if 'atcrd' in what:
        data['atcrd'] = []
    # Line 1: number of atoms
    # There may be empty lines ahead
    nblanks = 0
    while True:
        line = fobj.readline()
        if not line:
            return None
        else:
            if line.strip():
                break
            nblanks += 1
    if nblanks > 0 and not blank_ok:
        raise ParseDataError('Leading blank lines found.')
    try:
        natoms = int(line)
    except ValueError:
        msg = '1st line should be the number of atoms'
        raise ValueError(msg)
    if 'natoms' in data:
        data['natoms'] = natoms
    # Line 2: title
    line = fobj.readline()
    if 'title' in data:
        data['title'] = line.strip()
    # line 3-: geometry
    for _ in range(natoms):
        line = fobj.readline()
        cols = line.split()
        atom = cols[0]
        try:
            coord = [float(item)*__ang2au for item in cols[1:]]
        except ValueError:
            msg = 'Wrong XYZ geometry specification.'
            raise ValueError(msg)
        if 'atoms' in data:
            data['atoms'].append(atom)
        if 'atcrd' in data:
            data['atcrd'].append(coord)
    # Termination
    return data


def get_data(dfobj: FileXYZ,
             *qlabels: str,
             error_noqty: bool = True) -> TypeData:
    """Gets data from a XYZ file for each quantity label.

    Reads one or more full quantity labels from `qlabels` and returns
      the corresponding data.

    Parameters
    ----------
    dfobj
        Data file object.
    *qlabels
        List of full quantity labels to parse.
    error_noqty
        If True, error is raised if the quantity is not found.
    geom
        Geometry of interest (starting from 1).

    Returns
    -------
    dict
        For each `qlabel`, returns the corresponding data block.

    Raises
    ------
    TypeError
        Wrong type of data file object.
    ParseKeyError
        Missing required quantity in data block.
    IndexError
        State definition inconsistent with available data.
    QuantityError
        Unsupported quantity.
    """
    # First, check that the file is a correct instance
    if not isinstance(dfobj, FileXYZ):
        raise TypeError('FileXYZ instance expected')
    # Check if anything to do
    if len(qlabels) == 0:
        return None
    # Build Keyword List
    # ------------------
    # List of keywords
    keydata = {}
    qty_dict = {}
    geom = None
    for qlabel in qlabels:
        # Label parsing
        # ^^^^^^^^^^^^^
        qty_dict[qlabel] = ep.parse_qlabel(qlabel)
        qlab = qty_dict[qlabel][0]
        if qlab in ('atoms', 'atnum', 'atlab'):
            keydata[qlabel] = 'atoms'
        else:
            keydata[qlabel] = qlab
            if qlab in ('atcrd', 2, 1):
                if qty_dict[qlabel][1] == 'last':
                    geom = -1
                elif qty_dict[qlabel][1] == 'first':
                    geom = 1
                else:
                    geom = 0
    if geom is None:
        geom = 1
    # Data Extraction
    # ---------------
    # Use of set to remove redundant keywords
    datablocks = dfobj.read_data(*set(keydata.values()), geom=geom,
                                 raise_error=error_noqty)
    data = {}
    for qlabel in qlabels:
        key = keydata[qlabel]
        data[qlabel] = {}
        qlab = qty_dict[qlabel][0]
        if qlab in ('atnum', 'atlab'):
            data[qlabel]['data'] = convert_labsymb(qlab == 'atlab',
                                                   *datablocks[key])
        else:
            data[qlabel]['data'] = datablocks[key]
    return data
