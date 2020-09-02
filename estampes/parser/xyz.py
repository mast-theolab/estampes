"""Basic module to parse XYZ files

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
from estampes.base import ParseDataError, ParseKeyError, TypeData
import estampes.parser as ep
from estampes.tools.atom import convert_labsymb
from estampes.data.physics import PHYSFACT
# from estampes.parser import parse_qlabel

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
            raise FileNotFoundError('Formatted checkpoint not found')
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
                raise ValueError('1st line should contain the number of atoms')
            line = fobj.readline()
            self.__title = line.strip()
            for _ in range(self.__natoms):
                line = fobj.readline()
                if not line:
                    msg = 'File ended while reading molecular geometry'
                    raise ParseDataError(msg)
            line = fobj.readline()
            while line:
                self.__ngeoms += 1
                try:
                    natoms = int(line.strip())
                except ValueError:
                    fmt = '1st line of geometry {} should be number of atoms'
                    raise ValueError(fmt.format(self.__ngeoms))
                if natoms != self.__natoms:
                    fmt = 'Number of atoms in geometry {} does not match.'
                    raise ValueError(fmt.format(self.__ngeoms))
                for _ in range(2+self.__natoms):
                    line = fobj.readline()
                    if not line:
                        fmt = 'File ended while parsing geometry {}'
                        raise ParseDataError(fmt.format(self.__ngeoms))

    def read_data(self, *to_find: tp.Tuple[str],
                  geom: tp.Optional[int] = 1,
                  raise_error: bool = True) -> tp.Dict[str, tp.Any]:
        """Extracts data corresponding to the keys to find.

        Parameters
        ----------
        to_find : list
            List of keys to find.
        raise_error : bool
            Only raises error if `True`, otherwise proceeds silently.
        """
        datlist = {}
        if geom > self.__ngeoms:
            raise ValueError('Structure {} not available'.format(geom))
        if 'natoms' in to_find:
            datlist['natoms'] = self.__natoms
        rest = set(to_find) - {'natoms', 'title', 'atoms', 'atcrd'}
        if len(rest) > 0:
            raise ParseKeyError(rest[0])
        if {'title', 'atoms', 'atcrd'} & set(to_find):
            with open(self.__fname, 'r') as fobj:
                for _ in range((geom-1)*(self.__natoms+2)):
                    line = fobj.readline()
                line = fobj.readline()
                line = fobj.readline()
                if 'title' in to_find:
                    datlist['title'] = line.strip()
                if {'atoms', 'atcrd'} & set(to_find):
                    atoms = []
                    coord = []
                    for _ in range(self.__natoms):
                        line = fobj.readline()
                        cols = line.split()
                        atoms.append(cols[0])
                        coord.append([float(item)*__ang2au
                                      for item in cols[1:]])
                    if 'atoms' in to_find:
                        datlist['atoms'] = atoms
                    if 'atcrd' in to_find:
                        datlist['atcrd'] = coord
        return datlist

# ================
# Module Functions
# ================


def get_data(dfobj: FileXYZ,
             *qlabels: str,
             geom: int = 1) -> TypeData:
    """Gets data from a XYZ file for each quantity label.

    Reads one or more full quantity labels from `qlabels` and returns
      the corresponding data.

    Parameters
    ----------
    dfobj
        Data file object.
    *qlabels
        List of full quantity labels to parse.
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
    qty_dict = {}
    for qlabel in qlabels:
        # Label parsing
        # ^^^^^^^^^^^^^
        if qlabel in ('atoms', 'atnum', 'atlab'):
            qty_dict[qlabel] = 'atoms'
        else:
            qty_dict[qlabel] = ep.parse_qlabel(qlabel)[0]
    # Data Extraction
    # ---------------
    # Use of set to remove redundant keywords
    datablocks = dfobj.read_data(*list(set(qty_dict.values())), geom=geom,
                                 raise_error=True)
    data = {}
    for qlabel in qlabels:
        if qlabel in ('atnum', 'atlab'):
            data[qlabel] = convert_labsymb(qlabel == 'atlab',
                                           datablocks[qty_dict[qlabel]])
        else:
            data[qlabel] = datablocks[qty_dict[qlabel]]
    return data
