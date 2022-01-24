"""Basic module to parse CSV files.

This module provides basic classes and methods to parse CSV files.
Since the internal module of Python for CSV parsing is pretty dumb,
  an alternative parser is built by hand.

Methods
-------
get_data
    Gets data from a CSV file for each quantity label.

Classes
-------
FileCSV
    Main class to handle CSV files.
"""

import os
import re
import typing as tp

from estampes import parser as ep
from estampes.base import TypeQData, ParseKeyError, QuantityError

# ================
# Module Constants
# ================

txt_I = r'^\s*\b[+-]?\d+\b\s*$'
txt_F = r'^\s*[+-]?(?P<int>\d+)?[\.,]?(?P<dec>(?(int)\d*|\d+))' \
    + r'(?P<exp>[DdEe])?(?(exp)[+-]?\d+|)\s*$'
PATTERN_I = re.compile(txt_I)
PATTERN_F = re.compile(txt_F)


# ==============
# Module Classes
# ==============

class FileCSV(object):
    """Main class to handle XYZ files.

    Attributes
    ----------
    filename : str
        File name.
    full_version : tuple
        full version:
        * CSV
        * None

    Methods
    -------
    read_data(to_find, geom, raise_error)
        Extracts 1 or more data blocks from the CSV file.
    """
    def __init__(self, fname: str) -> None:
        self.filename = fname
        self.__ncols = 0
        self.__fields = None
        self.__pars = {
            'com': '',  # comment character
            'sep': '',  # separator character
            'comma': False,  # comma used as decimal point
            'fortran': False  # fortran exponential format
        }
        self.__headers = {
            'commments': [],  # Header comments
            'header': '',  # Header descriptor (non-comment)
        }
        self.analyze_file()

    @property
    def filename(self) -> str:
        """Gets or sets the filename associated to the CSV object."""
        return self.__fname

    @filename.setter
    def filename(self, name: str) -> None:
        if not os.path.exists(name):
            raise FileNotFoundError('CSV file not found')
        self.__fname = name

    @property
    def full_version(self) -> tp.Tuple[str, tp.Any]:
        """Returns the full version, for the parser interface"""
        return "CSV", None

    def analyze_file(self):
        """Analyzes a CSV file to get basic structural information.

        Analyzes a CSV file and return the necessary dialect.
        """
        CHARS_COMT = ('#', '!')
        # Separators characters.  ',' is just before space to avoid
        #   case of CSV file with comma decimals: we check first
        #   delimiters which are likely to be used in this case
        #   Space must be at then end, as it is likely to give many
        #   false positives.
        CHARS_SEP = (';', '\t', ',', ' ')

        c_com = ''
        c_sep = ''
        comma = False
        fortran = False
        header_comment = []
        header = None
        with open(self.__fname, 'r') as fobj:
            while True:
                line = fobj.readline()
                try:
                    char0 = line.lstrip()[0]
                    if char0 in CHARS_COMT:
                        if not c_com:
                            c_com = char0
                        elif char0 != c_com:
                            break
                        header_comment.append(line)
                    else:
                        break
                except IndexError:
                    # Empty line
                    pass
            for c in CHARS_SEP:
                if c in line:
                    c_sep = c
                    if c_sep == ' ':
                        c_sep = None
                    break
            fstr_I = None
            fstr_F = None
            iline = 1
            self.__ncols = len(line.split(c_sep))
            while True:
                cols = line.rstrip().split(c_sep)
                for item in cols:
                    if PATTERN_I.match(item):
                        fstr_I = item
                    elif PATTERN_F.match(item):
                        fstr_F = item
                if fstr_F is None:
                    if fstr_I is None and iline == 1:
                        header = line
                    line = fobj.readline()
                    iline += 1
                    if not line or iline == 10:
                        break
                else:
                    break
        if fstr_F is not None:
            comma = c_sep != ',' and ',' in line
            for item in line.rstrip().split(c_sep):
                if PATTERN_F.match(item) and 'd' in item.lower():
                    fortran = True
        self.__pars = {
            'com': c_com or '#',  # Sets a comment char anyway
            'sep': c_sep,
            'comma': comma,
            'fortran': fortran
        }
        self.__headers = {
            'comments': header_comment,
            'header': header
        }
        if self.__ncols == 1:
            self.__fields = ['y']
        elif self.__ncols == 2:
            self.__fields = ['x', 'y']
        else:
            fmt = 'y{{:0{}d}}'.format(len((str(self.__ncols-1))))
            self.__fields = ['x'] + [fmt.format(i)
                                     for i in range(1, self.__ncols)]
        if self.__headers['header']:
            cols = self.__headers['header'].split(c_sep)
            if len(cols) != self.__ncols:
                raise IndexError('Header seems inconsistent')
        elif self.__headers['comments']:
            for line in self.__headers['comments']:
                cols = line.replace(c_com, '', 1).strip().split(c_sep)
                if len(cols) != self.__ncols:
                    cols = None
                else:
                    break
        else:
            cols = None

        if cols is None:
            self.__labels = {key: None for key in self.__fields}
        else:
            self.__labels = {self.__fields[i]: cols[i].strip()
                             for i in range(self.__ncols)}

    def read_data(self, *to_find: tp.Sequence[str],
                  raise_error: bool = True) -> tp.Dict[str, tp.Any]:
        """Extracts data corresponding to the keys to find.

        Parameters
        ----------
        to_find : list
            List of keys to find.
        raise_error : bool
            Only raises error if `True`, otherwise proceeds silently.
        Raises
        ------
        ParseKeyError
            Key not found.
        """
        rest = set(to_find) - {'spcpar', 'spec'}
        if len(rest) > 0 and raise_error:
            raise ParseKeyError(rest[0])
        # Debugging flag
        FIX_FORT = True
        FIX_DEC = True
        datlist = {}
        if 'spcpar' in to_find:
            datlist['spcpar'] = self.__labels
            datlist['spcpar'].update({
                'func': None,
                'hwhm': None,
                'unitx': None,
                'unity': None
            })
        if 'spec' in to_find:
            datlist['spec'] = {key: [] for key in self.__fields}
            ccom = self.__pars['com']
            csep = self.__pars['sep']
            offset = len(self.__headers['comments'])
            if self.__headers['comments']:
                offset += 1
            with open(self.__fname, 'r') as fobj:
                fobj.seek(0)
                for _ in range(offset+1):
                    line = fobj.readline()
                while line:
                    line2 = line.strip()
                    if line2 and not line2.startswith(ccom):
                        if self.__pars['comma'] and FIX_DEC:
                            line2 = line2.replace(',', '.')
                        if self.__pars['fortran'] and FIX_FORT:
                            line2 = line2.replace('D', 'e')
                        for key, item in zip(self.__fields, line2.split(csep)):
                            res = item.strip()
                            if PATTERN_I.match(res):
                                func = int
                            elif PATTERN_F.match(res):
                                func = float
                            else:
                                func = str
                            try:
                                datlist['spec'][key].append(func(res))
                            except ValueError:
                                datlist['spec'][key].append(res)
                    line = fobj.readline()

        return datlist


# ================
# Module Functions
# ================

def get_data(dfobj: FileCSV,
             *qlabels: str,
             error_noqty: bool = True) -> TypeQData:
    """Gets data from a CSV file for each quantity label.

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
    if not isinstance(dfobj, FileCSV):
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
        qtag, qopt = ep.parse_qlabel(qlabel)[:2]
        if qtag != 'anyspc':
            raise QuantityError('Only spectroscopic data supported.')
        qty_dict[qlabel] = qopt.lower()
    # Data Extraction
    # ---------------
    # Use of set to remove redundant keywords
    datablocks = dfobj.read_data(*list(set(qty_dict.values())),
                                 raise_error=error_noqty)
    data = {}
    for qlabel in qlabels:
        data[qlabel] = datablocks[qty_dict[qlabel]]
    return data
