"""Low-level operations on Gaussian formatted checkpoint files.

Provides low-level interfaces to manipulate/extract data in Gaussian
formatted checkpoint files.

Notes
-----
* While this module contains some basic error checks, it is intended to
  be low-level, with no significant performance impacts.  As a low-level
  module, it should not be accessible to users but wrapped by developers
  in charge of controlling that the queries are computationally AND
  physically sound.

"""

import os  # Used for file existence check
import re  # Used to find keys in fchk file
from tempfile import TemporaryFile
from shutil import copyfileobj
import typing as tp
from math import ceil, sqrt

from estampes.base import ArgumentError, ParseDataError, ParseKeyError, \
    QuantityError, TypeDFChk, TypeQData, TypeQInfo, QData, QLabel
from estampes.data import property as edpr
from estampes.parser.functions import parse_qlabels, reshape_dblock

# ================
# Module Constants
# ================

TypeKword = tp.Dict[str, tp.Tuple[str, int, int]]
TypeQKwrd = tp.Union[str, tp.List[str]]

NCOLS_FCHK = {  # Maximum number of columns per type in fchk
    'C': 5,  # Number of columns for character data per line
    'R': 5,  # Number of columns for float data per line
    'I': 6  # Number of columns for integer data per line
}
FCONV_FCHK = {  # Conversion function for each type
    'C': str,
    'I': int,
    'R': float
}
DFMT_FCHK = {  # Data format for each type
    'C': '{:12s}',
    'I': '{:12d}',
    'R': '{:16.8E}'
}

# ==============
# Module Classes
# ==============


class FChkIO(object):
    """Main class to handle formatted checkpoint file operations.

    Main class to manage the parsing and formatting of data stored in
    Gaussian formatted checkpoint file.
    """

    def __init__(self, fname: str,
                 load_keys: bool = True) -> None:
        """Initialize Gaussian fchk IO instance.

        Initializes an instance of FChkIO.

        Parameters
        ----------
        fname
            Formatted checkpoint file name.
        load_keys
            Preloads keys from the fchk file to speed up search.
        """
        self.filename = fname
        if load_keys:
            self.__keys = self.__store_keys()  # type: tp.Optional[TypeKword]
        else:
            self.__keys = None
        try:
            key = 'Gaussian Version'
            qtag = 'swver'
            qdict = {qtag: QLabel(quantity=qtag)}
            qkwrd = {qtag: key}
            qdata = self.read_data(key)
            self.__gversion = parse_data(qdict, qkwrd, qdata)[qtag].data
        except ParseKeyError:
            self.__gversion = {'major': None, 'minor': None}

    @property
    def filename(self) -> str:
        """Get or set the filename associated to the FChk object."""
        return self.__fname

    @filename.setter
    def filename(self, name: str) -> None:
        if not os.path.exists(name):
            raise FileNotFoundError('Formatted checkpoint not found')
        self.__fname = name

    @property
    def version(self) -> tp.Dict[str, str]:
        """Return the version of Gaussian used to generate the FChk file.

        Returns the version as a dictionary with two keys:

        major
            Major revision.
        minor
            Minor revision

        Notes
        -----
        Earlier versions of Gaussian did not support this so this may be
        empty.
        """
        return self.__gversion

    def show_keys(self):
        """Return the available keys (only if loaded)."""
        if self.__keys is None:
            return None
        else:
            return sorted(self.__keys.keys())

    @property
    def full_version(self) -> tp.Tuple[str, tp.Any]:
        """Return the full Gaussian version, for the parser interface.

        Returns the full version of Gaussian used to generate the log file.
        Available keys:

        major
            Major revision.
        minor
            Minor revision.
        mach
            Processor architecture for which Gaussian was compiled.
        release
            Release date of this version of Gaussian.
        """
        return "Gaussian", self.__gversion

    def read_data(self,
                  *to_find: tp.Tuple[str],
                  raise_error: bool = True) -> TypeQData:
        """Extract data corresponding to the keys to find.

        Parameters
        ----------
        to_find
            Key or list of keys to find.
        raise_error
            Only raises error if True, otherwise proceeds silently.

        Raises
        ------
        ParseKeyError
            Key not found.
        """
        keylist = []  # List of keywords to search
        datlist = {}  # type: TypeQData # List of data

        # Fast Search
        # -----------
        # Uses the data in __keys to find pointers.
        if self.__keys is not None:
            # Build keyword list
            # ^^^^^^^^^^^^^^^^^^
            for item in to_find:
                if not item.strip() in self.__keys:
                    if raise_error:
                        raise ParseKeyError(item)
                else:
                    keylist.append([item, *self.__keys[item]])
            # Sort the keys by order of appearance
            keylist.sort(key=lambda x: x[3])
            with open(self.filename, 'r', encoding="utf-8") as fobj:
                for item in keylist:
                    key, dtype, ndata, fpos = item
                    fobj.seek(fpos)
                    line = fobj.readline()
                    datlist[key] = self.__read_datablock(fobj, line, dtype,
                                                         ndata)

        # Sequential Search
        # -----------------
        # Looks for keywords sequentially while reading file
        else:
            nkeys = len(to_find)
            with open(self.filename, 'r', encoding="utf-8") as fobj:
                line = fobj.readline()
                while line and nkeys > 0:
                    line = fobj.readline()
                    for key in to_find:
                        if line.startswith(key):
                            datlist[key] = self.__read_datablock(fobj, line)
                            nkeys -= 1
            remaining = list(set(to_find) - set(datlist))
            if len(remaining) > 0 and raise_error:
                raise ParseKeyError(remaining[0])

        return datlist

    def write_data(self,
                   data: tp.Dict[str, tp.Sequence[tp.Any]],
                   new_file: tp.Optional[str] = None,
                   error_key: bool = True,
                   error_size: bool = True) -> None:
        """Write data corresponding to the keys to find.

        Reads a dictionary of keys and overwrites the data present in
        the file.
        If the key is not present or the size is inconsistent with the
        data present in the file, an error is raised, except if
        `error_key` or `error_size` are False, respectively.

        Parameters
        ----------
        data
            Dictionary with the replacement data for each key.
        new_file
            Name of the file where data are printed.
            If none, the internal file is overwritten.
        error_key
            If true, raises error if key not found.
        error_size
            If true, raises error for inconsistent size.

        Raises
        ------
        ParseKeyError
            Key not found.
        IndexError
            Inconsistency in size between old and new data for a key.
        """
        fmt_scal = {
            'I': '{:<40s}   I     {:12d}\n',
            'R': '{:<40s}   R     {:22.15E}\n',
            'C': '{:<40s}   C     {:12s}\n',
        }
        fmt_head = '{:<40s}   {:1s}   N={:12d}\n'

        # Compared available keys with those from new data set
        # ----------------------------------------------------
        keys_ok = {}
        if self.__keys is not None:
            # Uses the data in __keys to find pointers.
            # Check if overlap between data and stored keys
            keys = set(self.__keys) & set(data)
            for key in keys:
                keys_ok[self.__keys[key][-1]] = key
            keys_no = list(set(data) - keys)
        else:
            nkeys = len(data)
            keys_no = data.keys()
            with open(self.filename, 'r', encoding="utf-8") as fobj:
                line = fobj.readline()
                while line and nkeys > 0:
                    fpos = 0
                    if line[0] != ' ':
                        for index, key in enumerate(keys_no):
                            if line.startswith(key):
                                keys_ok[fpos] = keys_no.pop(index)
                                nkeys -= 1
                    fpos += len(line)
                    line = fobj.readline()
        if keys_no and error_key:
            raise ParseKeyError(', '.join(keys_no))

        # Now set where data are to be saved
        if new_file is None:
            fdest = TemporaryFile()
        else:
            fdest = open(new_file, 'w', encoding="utf-8")

        # Now let us copy the content of the internal file in destination
        # For each key retained after the analysis, we replace with the new
        # data
        with open(self.filename, 'r', encoding="utf-8") as fsrc:
            fpos = 0
            for line in fsrc:
                if fpos in keys_ok:
                    key = keys_ok[fpos]
                    dtype, ndat_ref, ncols, nlin_ref = self.__info_block(line)
                    ndat_new = len(data[key])
                    if ndat_ref == 0:
                        if ndat_new > 1 and error_size:
                            raise IndexError(f'Inconsistency with {key}')
                        else:
                            fdest.write(fmt_scal[dtype].format(key, data[key]))
                    else:
                        fdest.write(fmt_head.format(key, ndat_new))
                        for i in range(0, ndat_new, ncols):
                            N = min(ncols, ndat_new-i)
                            fmt = N*DFMT_FCHK[dtype] + '\n'
                            fdest.write(fmt.format(data[key][i:i+N]))
                        for _ in range(nlin_ref):
                            line = next(fsrc)
                            fpos += len(line)
                else:
                    fdest.write(line)
                    fpos += len(line)

            # Copy back file if requested
            if new_file is not None:
                fdest.seek(0)
                fsrc.seek(0)
                copyfileobj(fdest, fsrc)

    def __store_keys(self) -> TypeKword:
        """Store the keys in the fchk to speed up search.

        Loads the keys present in the file and pointers to their
        position to speed up their search.
        Data type and block information are also stored.

        Returns
        -------
        dict
            For each key, returns a tuple with:
            1. data type (I, R, C)
            2. Number of values (0 for scalar)
            3. position in files
        """
        to_search = re.compile(r'''
            (?P<title>[\w\s/\-]+?)\s*  # Key
            \b(?P<type>[IRC])\b\s*  # Data type
            (?P<block>N=)?\s+  # N= only set for non-scalar data
            (?P<value>[\d\-\+\.E]+)  # Block size (N=) or scalar value
            $''', re.VERBOSE)

        keys = {}
        with open(self.filename, 'r', encoding="utf-8") as fobj:
            fpos = 0
            for line in fobj:
                res = to_search.match(line)
                if res:
                    nval = int(res.group(3) and res.group(4) or 0)
                    keys[res.group(1)] = (res.group(2), nval, fpos)
                fpos += len(line)
        return keys

    def __info_block(self, line: tp.Optional[str] = None,
                     datatype: tp.Optional[str] = None,
                     numdata: tp.Optional[int] = None
                     ) -> tp.List[tp.Any]:
        """Extract information on a given block.

        Extracts information on a block, either from the line or data
        in arguments.

        Parameters
        ----------
        line
            Starting line of a block.
        datatype
            Type of data.
        numdata
            Number of data.

        Returns
        -------
        str
            Type of data.
        int
            Number of data.
        int
            Number of columns.
        int
            Number of lines.

        Raises
        ------
        ArgumentError
            Arguments are insufficient to generate the data.
        ParseDataError
            Unsupported data types.
        """
        if datatype is None and line is None:
            raise ArgumentError('line and datatype cannot be both absent')
        # If data type unknown, line has not been parsed
        if datatype is None:
            cols = line.split()
            if 'N=' in line:
                dtype = cols[-3]
                ndata = int(cols[-1])
            else:
                dtype = cols[-2]
                ndata = 0
        else:
            dtype = datatype
            ndata = numdata
        # Sets parameters:
        try:
            ncols = NCOLS_FCHK[dtype]
        except KeyError as err:
            raise ParseDataError(dtype, 'Unsupported data type') from err
        nlines = int(ceil(ndata/ncols))
        return dtype, ndata, ncols, nlines

    def __read_datablock(self, fobj: tp.TextIO,
                         line: str,
                         datatype: tp.Optional[str] = None,
                         numdata: tp.Optional[int] = None
                         ) -> tp.List[tp.Any]:
        """Read a data block in the formatted checkpoint file.

        Reads a data block from a Gaussian formatted checkpoint file.
        The file "cursor" should be at the "title/section line" and the
        content stored in 'line'.

        Parameters
        ----------
        fobj
            Opened file.
        line
            Current line read from file object.
        datatype
            Type of the scalar or data block.
        numdata
            Size of the data block (0 if scalar).

        Raises
        ------
        ParseDataError
            Unsupported data type.

        Notes
        -----
        * The function uses readline() to extract the actual block.
        * The parsing is mostly format-free for simplicity.

        .. [1] http://gaussian.com/interfacing/?tabid=3
        """
        dtype, _, _, nlines = self.__info_block(line, datatype, numdata)
        # Sets parameters:
        try:
            fconv = FCONV_FCHK[dtype]
        except KeyError as err:
            raise ParseDataError(dtype, 'Unsupported data type') from err
        # Data Extraction
        # ---------------
        if nlines == 0:
            # Scalar
            data = [fconv(line.split()[-1])]
        else:
            # Data Block
            # We use a slightly different scheme for C since Gaussian cuts
            # arbitrarily strings in the middle in the format
            if dtype == 'C':
                block = ''
                for _ in range(nlines):
                    block += fobj.readline().rstrip('\n')  # Remove newline
                data = block.split()
            else:
                data = []
                for _ in range(nlines):
                    line = fobj.readline()
                    data.extend([fconv(item) for item in line.split()])
        return data

# ================
# Module Functions
# ================


def qlab_to_kword(qlab: QLabel) -> TypeQKwrd:
    """Return the keyword(s) relevant for a given quantity.

    Returns the keyword corresponding to the block containing the
    quantity of interest and the list all keywords of interest for
    possible conversions.

    Parameters
    ----------
    qlab
        Quantity label.

    Returns
    -------
    list
        List of keywords for the data to extract.
    list
        Information needed for extracting the quantity of interest.
        1. keyword in the formatted checkpoint file
        2. position of the first element in the data block
        3. offsets for "sub-block" storage (data in several blocks)

    Raises
    ------
    NotImplementedError
        Missing features.
    QuantityError
        Unsupported quantity.
    ValueError
        Unsupported case.

    Notes
    -----
    - `n` refers to all available states.
    """
    keywords = []
    keyword = None
    if qlab.label == 'natoms':
        keyword = 'Number of atoms'
    elif qlab.label == 'nvib':
        keyword = 'Number of Normal Modes'
    elif qlab.label == 'atmas':
        keyword = 'Real atomic weights'
    elif qlab.label == 'atnum':
        keyword = 'Atomic numbers'
    elif qlab.label == 'molsym':
        raise NotImplementedError()
    elif qlab.label == 'atcrd' or qlab.label == 2:
        keyword = 'Current cartesian coordinates'
    elif qlab.label == 'hessvec':
        if qlab.level == 'H':
            keyword = 'Vib-Modes'
        else:
            keyword = 'Anharmonic Vib-Modes'
    elif qlab.label == 'hessdat':
        if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
            if qlab.level == 'H':
                keyword = 'Vib-E2'
                keywords.append('Number of Normal Modes')
            else:
                keyword = 'Anharmonic Vib-E2'
                keywords.append('Anharmonic Number of Normal Modes')
    elif qlab.label == 'swopt':
        keyword = 'Route'
    elif qlab.label == 'swver':
        keyword = 'Gaussian Version'
    elif qlab.label == 'fcdat':
        raise NotImplementedError()
    elif qlab.label == 'vptdat':
        raise NotImplementedError()
    elif qlab.label in ('dipstr', 'rotstr'):
        if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
            if qlab.kind == 'H':
                keyword = 'Vib-E2'
                keywords = ['Number of Normal Modes']
            else:
                keyword = 'Anharmonic Vib-E2'
                keywords = ['Anharmonic Number of Normal Modes']
        else:
            keywords = ['ETran scalars']
            keyword = 'ETran state values'
    else:
        if isinstance(qlab.rstate, tuple):
            keyword = 'ETran state values'
            if qlab.label == 1 and qlab.derord == 0:
                keywords = ['ETran scalars', 'SCF Energy']
        else:
            if qlab.label == 1:
                if qlab.derord == 0:
                    if qlab.rstate == 'c':
                        keyword = 'Total Energy'
                        del keywords[:]
                    elif isinstance(qlab.rstate, int):
                        if qlab.rstate == 0:
                            keyword = 'SCF Energy'
                        else:
                            keyword = 'ETran state values'
                        keywords = ['Total Energy', 'ETran scalars']
                elif qlab.derord == 1:
                    if qlab.dercrd is None or qlab.dercrd == 'X':
                        if qlab.rstate == 'c' or isinstance(qlab.rstate, int):
                            keyword = 'Cartesian Gradient'
                elif qlab.derord == 2:
                    if qlab.dercrd is None or qlab.dercrd == 'X':
                        if qlab.rstate == 'c' or isinstance(qlab.rstate, int):
                            keyword = 'Cartesian Force Constants'
            elif qlab.label == 50:
                keywords = ['ETran scalars']
                keyword = 'Nonadiabatic coupling'
            elif qlab.label == 91:
                raise NotImplementedError()
            elif qlab.label == 92:
                keyword = 'RotTr to input orientation'
            elif qlab.label == 93:
                keyword = 'RotTr to input orientation'
            elif qlab.label == 101:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Dipole Moment'
                elif qlab.derord == 1:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Dipole Derivatives'
            elif qlab.label == 102:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        raise ParseDataError('Magnetic dipole not available')
                elif qlab.derord == 1:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'AAT'
            elif qlab.label == 103:
                raise NotImplementedError()
            elif qlab.label == 104:
                raise NotImplementedError()
            elif qlab.label == 105:
                raise NotImplementedError()
            elif qlab.label == 106:
                raise NotImplementedError()
            elif qlab.label == 107:
                raise NotImplementedError()
            elif qlab.label == 201:
                raise NotImplementedError()
            elif qlab.label == 202:
                raise NotImplementedError()
            elif qlab.label == 203:
                raise NotImplementedError()
            elif qlab.label == 204:
                raise NotImplementedError()
            elif qlab.label == 205:
                raise NotImplementedError()
            elif qlab.label == 206:
                raise NotImplementedError()
            elif qlab.label == 207:
                raise NotImplementedError()
            elif qlab.label == 208:
                raise NotImplementedError()
            elif qlab.label == 209:
                raise NotImplementedError()
            elif qlab.label == 300:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Frequencies for FD properties'
                    else:
                        msg = 'Incident frequencies not available'
                        raise ParseDataError(msg)
                else:
                    keywords = ['Number of atoms']
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Frequencies for DFD properties'
                    else:
                        msg = 'Incident frequencies not available'
                        raise ParseDataError(msg)
            elif qlab.label == 301:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Alpha(-w,w)'
                elif qlab.derord == 1:
                    keywords = ['Number of atoms']
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Derivative Alpha(-w,w)'
            elif qlab.label == 302:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'FD Optical Rotation Tensor'
                elif qlab.derord == 1:
                    keywords = ['Number of atoms']
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Derivative FD Optical Rotation Tensor'
            elif qlab.label == 303:
                if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                    raise ParseDataError('Alpha(w,0) not available')
            elif qlab.label == 304:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'D-Q polarizability'
                elif qlab.derord == 1:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        keyword = 'Derivative D-Q polarizability'
            elif qlab.label == 305:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        raise NotImplementedError()
                elif qlab.derord == 1:
                    keywords = ['Number of atoms']
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        raise NotImplementedError()
            elif qlab.label == 306:
                if qlab.derord == 0:
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        raise NotImplementedError()
                elif qlab.derord == 1:
                    keywords = ['Number of atoms']
                    if isinstance(qlab.rstate, int) or qlab.rstate == 'c':
                        raise NotImplementedError()
            else:
                raise QuantityError('Unknown quantity')
    keywords.insert(0, keyword)
    return keyword, keywords


def _parse_electrans_data(qlab: QLabel, dblocks: TypeDFChk,
                          kword: str) -> TypeQData:
    """Sub-function to parse electronic-transition related data.

    Parses and returns data for a given quantity related to an
    electronic transition.

    Parameters
    ----------
    qlab
        Quantity label.
    dblocks
        Data blocks, by keyword.
    kword
        Keyword for quantity of interest.

    Returns
    -------
    dict
        Data for each quantity.

    Raises
    ------
    ParseKeyError
        Missing required quantity in data block.
    IndexError
        State definition inconsistent with available data.
    QuantityError
        Unsupported quantity.
    """
    # ETran Scalar Definition
    # -----------------------
    # Check that ETran scalars are present and parse relevant values
    key = 'ETran scalars'
    if key in dblocks:
        # Structure of ETran scalars
        # 1. Number of electronic states
        # 2. Number of scalar data stored per state
        # 3. 1 of R==L transition matrix, 2 otherwise
        # 4. Number of header words (irrelevant in fchk)
        # 5. State of interest
        # 6. Number of deriv. (3*natoms + 3: electric field derivatives)
        (nstates, ndata, _, _, iroot,
            _) = [item for item in dblocks[key][:6]]
    else:
        raise ParseKeyError('Missing scalars definition')
    # States Information
    # ------------------
    initial, final = qlab.rstate
    if initial != 0:
        if final != 0:
            raise IndexError('Unsupported transition')
        else:
            initial, final = final, initial
    # Quantity-specific Treatment
    # ---------------------------
    if qlab.label == 2:
        key = 'SCF Energy'
        if key not in dblocks:
            raise ParseKeyError('Missing ground-state energy')
        energy0 = dblocks[key]
        if final == 'a':
            data = [dblocks[key][i*ndata]-energy0 for i in range(nstates)]
        else:
            fstate = final == 'c' and iroot or final
            if fstate > nstates:
                raise IndexError('Missing electronic state')
            data = float(dblocks[key][(fstate-1)*ndata]) - energy0
    elif qlab.label in (101, 102, 103):
        lqty = edpr.property_data(qlab.label).dim
        if qlab.label == 101:
            if qlab.kind == 'len':
                offset = 1
            else:
                offset = 4
        elif qlab.label == 102:
            offset = 7
        else:
            offset = 10
        if qlab.derord == 0:
            if final == 'a':
                data = [dblocks[kword][i*ndata+offset:i*ndata+offset+lqty]
                        for i in range(nstates)]
            else:
                fstate = final == 'c' and iroot or final
                if fstate > nstates:
                    raise IndexError('Missing electronic state')
                i0 = (fstate-1)*ndata + offset
                data = dblocks[kword][i0:i0+lqty]
    else:
        raise QuantityError('Unsupported quantity')
    return data


def _parse_freqdep_data(qlab: QLabel, dblocks: TypeDFChk,
                        kword: str) -> TypeQData:
    """Sub-function to parse data on frequency-dependent properties.

    Parses and returns data on a specific property for one or more
    incident frequencies.

    Parameters
    ----------
    qlab
        Quantity label
    dblocks
        Data blocks, by keyword.
    kword
        Keyword for quantity of interest.

    Returns
    -------
    dict
        Data for each quantity.

    Raises
    ------
    ParseKeyError
        Missing required quantity in data block.
    IndexError
        Error with definition of incident frequency.
    QuantityError
        Unsupported quantity.
    """
    # Check Incident Frequency
    # ------------------------
    if qlab.kind is None:
        qopt_ = 0
    elif not isinstance(qlab.kind, int):
        raise IndexError()
    else:
        qopt_ = qlab.kind
    # Quantity-specific Treatment
    # ---------------------------
    # Check size of derivatives is requested
    if qlab.derord == 0:
        nder = 1
    elif qlab.derord == 1:
        key = 'Number of atoms'
        if key not in dblocks:
            raise ParseKeyError('Missing number of atoms')
        natoms = dblocks[key]
        nder = 3*natoms
    else:
        raise IndexError('Unsupported derivative order')

    if qlab.label == 301:
        lqty = 9*nder
    elif qlab.label == 302:
        lqty = 9*nder
    elif qlab.label == 303:
        lqty = 9*nder
    elif qlab.label == 304:
        lqty = 18*nder
    elif qlab.label == 305:
        lqty = 18*nder
    elif qlab.label == 306:
        lqty = 18*nder
    else:
        raise QuantityError('Unsupported quantity')
    lblock = len(dblocks[kword])
    ndata = lblock // lqty  # Assumed block is correctly built
    if qopt_ == 0:
        data = [dblocks[kword][i*lqty:(i+1)*lqty] for i in range(ndata)]
    else:
        if qopt_ > ndata:
            raise IndexError('Incident frequency index out of range')
        data = dblocks[kword][(qopt_-1)*lqty:qopt_*lqty]

    return data


def parse_data(qdict: TypeQInfo,
               qlab2kword: tp.Dict[str, str],
               datablocks: TypeDFChk,
               gver: tp.Optional[tp.Tuple[str, str]] = None,
               raise_error: bool = True) -> TypeQData:
    """Parse data arrays to extract specific quantities.

    Parses data array to extract relevant information for each quantity.

    Parameters
    ----------
    qdict
        Dictionary of quantities.
    qlab2kword
        Connects each qlabel with the section title in fchk file.
    datablocks
        Data blocks, by keyword.
    gver
        Gaussian version (major, minor).
    raise_error
        If True, error is raised if the quantity is not found.

    Returns
    -------
    dict
        Data for each quantity.

    Raises
    ------
    ParseKeyError
        Missing required quantity in data block.
    IndexError
        State definition inconsistent with available data.
    ValueError
        Data inconsistency with respect to shape.
    QuantityError
        Unsupported quantity.
    """
    def empty_cases_ok(qlabel: QLabel) -> bool:
        return qlabel.label == 'nvib'
    dobjs = {}
    for qkey, qlabel in qdict.items():
        msg_noqty = f'Missing quantity "{qlabel}" in file'
        kword = qlab2kword[qkey]
        # Basic Check: main property present
        # -----------
        if kword not in datablocks and not empty_cases_ok(qlabel):
            if raise_error:
                raise ParseKeyError(msg_noqty)
            else:
                dobjs[qkey] = None
                continue
        dobjs[qkey] = QData(qlabel)
        # Basic Properties/Quantities
        # ---------------------------
        if qlabel.label == 'natoms':
            dobjs[qkey].set(data=int(datablocks[kword][0]))
        elif qlabel.label in ('atcrd', 2):
            dobjs[qkey].set(data=reshape_dblock(datablocks[kword], (3, )))
        elif qlabel.label in ('atmas', 'atnum'):
            dobjs[qkey].set(data=datablocks[kword])
        elif qlabel.label == 'swopt':
            dobjs[qkey].set(data=' '.join(datablocks[kword]))
        elif qlabel.label == 'molsym':
            raise NotImplementedError()
        elif qlabel.label == 'swver':
            pattern = re.compile(r'(?:(\w+)-)?(\w{3})Rev-?([\w.+]+)')
            res = re.match(pattern, ''.join(datablocks[kword])).groups()
            data = {'major': res[-2], 'minor': res[-1]}
            if len(res) == 3:
                data['system'] = res[0]
            data['release'] = None
            dobjs[qkey].set(data=data)
        # Vibrational Information
        # -----------------------
        # Technically state should be checked but considered irrelevant.
        elif qlabel.label == 'nvib':
            if kword in datablocks:
                dobjs[qkey].set(data=int(datablocks[kword][0]))
            else:
                # For a robust def of nvib, we need the symmetry and
                #   the number of frozen atoms. For now, difficult to do.
                raise NotImplementedError()
        elif qlabel.label == 'hessvec':
            if kword in datablocks:
                dobjs[qkey].set(dtype='L.M^{-1/2}')
                dobjs[qkey].set(data=datablocks[kword])
        elif qlabel.label == 'hessdat':
            if qlabel.kind == 'H':
                key = 'Number of Normal Modes'
            else:
                key = 'Anharmonic Number of Normal Modes'
            if key not in datablocks:
                raise ParseKeyError('Missing necessary dimension')
            ndat = int(datablocks[key][0])
            if qlabel.kind == 'freq':
                offset = 0
            elif qlabel.kind == 'redmas':
                offset = ndat
            dobjs[qkey].set(data=datablocks[kword][offset:offset+ndat])
        # Vibronic Information
        # --------------------
        elif qlabel.label == 'fcdat':
            raise NotImplementedError()
        # Anharmonic Information
        # ----------------------
        elif qlabel.label == 'vptdat':
            raise NotImplementedError()
        # State(s)-dependent quantities
        # -----------------------------
        else:
            # Transition moments
            # ^^^^^^^^^^^^^^^^^^
            if isinstance(qlabel.rstate, tuple):
                dobjs[qkey].set(
                    data=_parse_electrans_data(qlabel, datablocks, kword))
            # States-specific Quantities
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^
            else:
                key = 'ETran scalars'
                if key in datablocks:
                    (nstates, ndata, _, _, iroot,
                        _) = [item for item in datablocks[key][:6]]
                else:
                    nstates = 1
                    ndata = 16  # TODO: This may not work in the future.
                    iroot = 0
                curr_sta = qlabel.rstate == 'c' or qlabel.rstate == iroot
                # Only energy is currently computed for all states:
                if qlabel.rstate == 'a' and qlabel.kind == 2:
                    data = [float(datablocks[kword][i*ndata])
                            for i in range(nstates)]
                # Data for current electronic states
                elif curr_sta:
                    if qlabel.kind in ('dipstr', 'rotstr'):
                        if qlabel.level == 'H':
                            key = 'Number of Normal Modes'
                        else:
                            key = 'Anharmonic Number of Normal Modes'
                        if key not in datablocks:
                            raise ParseKeyError('Missing necessary dimension')
                        ndat = int(datablocks[key][0])
                        if qlabel.kind == 'dipstr':
                            offset = 7*ndat
                        else:
                            offset = 8*ndat
                        dobjs[qkey].set(
                            data=datablocks[kword][offset:offset+ndat])
                    elif qlabel.label == 1:
                        if qlabel.derord == 0:
                            dobjs[qkey].set(data=datablocks[kword][0])
                        elif qlabel.derord in (1, 2):
                            dobjs[qkey].set(data=datablocks[kword])
                        else:
                            raise NotImplementedError()
                    elif qlabel.label == 92:
                        dobjs[qkey].set(data=datablocks[kword][:9])
                    elif qlabel.label == 93:
                        dobjs[qkey].set(data=datablocks[kword][9:])
                    elif qlabel.label in (50, 91):
                        dobjs[qkey].set(data=datablocks[kword])
                    elif qlabel.label == 101:
                        if qlabel.derord in (0, 1):
                            dobjs[qkey].set(data=datablocks[kword])
                        else:
                            raise NotImplementedError()
                    elif qlabel.label == 102:
                        if qlabel.derord == 1:
                            dobjs[qkey].set(data=datablocks[kword])
                        else:
                            raise NotImplementedError()
                    elif qlabel.label == 300:
                        if qlabel.derord in (0, 1):
                            if qlabel.kind == 0:
                                dobjs[qkey].set(data=datablocks[kword])
                        else:
                            raise NotImplementedError()
                    else:
                        raise NotImplementedError()

    return dobjs


def get_data(dfobj: FChkIO,
             *qlabels: str,
             error_noqty: bool = True,
             **keys4qlab) -> TypeQData:
    """Get data from a FChk file for each quantity label.

    Reads one or more full quantity labels from `qlabels` and returns
    the corresponding data.

    Parameters
    ----------
    dfobj
        Formatted checkpoint file as `FChkIO` object.

    *qlabels
        List of full quantity labels to parse.
    error_noqty
        If True, error is raised if the quantity is not found.
    **keys4qlab
        Aliases for the qlabels to be used in returned data object.

    Returns
    -------
    dict
        Data for each quantity.

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
    if not isinstance(dfobj, FChkIO):
        raise TypeError('FChkIO instance expected')
    # Check if anything to do
    if len(qlabels) == 0 and len(keys4qlab) == 0:
        return None
    # Build Keyword List
    # ------------------
    # Build full list of qlabels
    qty_dict, dupl_keys = parse_qlabels(qlabels, keys4qlab)
    # List of keywords
    full_kwlist = []
    main_kwlist = {}
    for qkey, qlabel in qty_dict.items():
        # Label parsing
        # ^^^^^^^^^^^^^
        keyword, keywords = qlab_to_kword(qlabel)
        if keyword is not None:
            full_kwlist.extend(keywords)
            main_kwlist[qkey] = keyword
    # Check if list in the end is not empty
    if not main_kwlist:
        raise QuantityError('Unsupported quantities')
    # Data Extraction
    # ---------------
    # Use of set to remove redundant keywords
    datablocks = dfobj.read_data(*list(set(full_kwlist)), raise_error=False)
    # Data Parsing
    # ------------
    gver = (dfobj.version['major'], dfobj.version['minor'])
    try:
        data = parse_data(qty_dict, main_kwlist, datablocks, gver,
                          error_noqty)
    except (QuantityError, NotImplementedError) as err:
        raise QuantityError('Unsupported quantities') from err
    except (ParseKeyError, IndexError) as err:
        raise IndexError('Missing data in FChk') from err

    # Fill redundant keys
    if dupl_keys:
        for item, key in dupl_keys.items():
            data[item] = data[key]

    return data


def get_hess_data(dfobj: tp.Optional[FChkIO] = None,
                  get_evec: bool = True,
                  get_eval: bool = True,
                  pre_data: tp.Optional[TypeQData] = None,
                  force_calc: bool = False
                  ) -> tp.Tuple[tp.Any]:
    """Get or build Hessian data (eigenvectors and values).

    This function retrieves or builds the eigenvectors and eigenvalues.
    Contrary to ``get_data`` which only looks for available data, this
    functions looks for alternative forms to build necessary data.
    It also returns a Numpy array instead of Python lists.
    Preloaded data can be provided to avoid duplicating extraction
    queries.  They are expected in the same format as given by
    FChkIO methods.
    By default, it uses data read from formatted checkpoint file, as
    they are sufficiently accurate and contains multiple corrections:
    projection of gradients, rotation/translation...

    Parameters
    ----------
    dfobj
        Formatted checkpoint file as `FChkIO` object.
    get_evec
        Return the eigenvectors.
    get_eval
        Return the eigenvalues.
    pre_data
        Database with quantities already loaded from previous queries.
    force_calc
        Force the computation of the eigenvectors and eigenvalues.

    Returns
    -------
    :obj:numpy.ndarray
        Eigenvectors (None if not requested).
    :obj:numpy.ndarray
        Eigenvalues (None if not requested).

    Raises
    ------
    ValueError
        Inconsitent values given in input.
    IOError
        Error if file object not set but needed.
    IndexError
        Quantity not found.

    Notes
    -----
    * Numpy is needed to run this function
    * Data can be given in argument or will be extracted from `dfobj`
    """
    import sys
    import numpy as np
    from estampes.tools.math import square_ltmat
    from estampes.data.physics import phys_fact

    def build_evec(fccart, atmas, get_evec, get_eval, nvib=None):
        """Build evec and eval from force constants matrix."""
        inv_sqmas = 1./np.sqrt(atmas)
        nat3 = fccart.shape[0]
        ffx = np.einsum('i,ij,j->ij', inv_sqmas, fccart, inv_sqmas)
        hessval, hessvec = np.linalg.eigh(ffx)
        vibs = np.full(nat3, True)
        freqs = []
        for i, val in enumerate(hessval):
            freq = sqrt(abs(val)*phys_fact('fac2au'))
            if freq < 10.0:
                vibs[i] = False
            else:
                if val < 0:
                    freqs.append(-freq)
                else:
                    freqs.append(freq)
        if nvib is not None:
            if np.count_nonzero(vibs) != nvib:
                msg = 'Unable to identify vibrations from rot/trans'
                raise QuantityError(msg)
        eigvec = norm_evec(hessvec[:, vibs].T) if get_evec else None
        eigval = freqs if get_eval else None

        return eigvec, eigval

    def convert_evec(hessvec, atmas=None, natoms=None, form='L.M^{-1/2}'):
        """Convert eigenvector based on form and available data."""
        eigvec = None
        if form in ('L.M^-1/2', 'L/M^1/2', 'L.M^{-1/2}', 'L/M^{1/2}'):
            if atmas is not None:
                nat3 = atmas.size
                eigvec = norm_evec(np.einsum(
                    'ij,j->ij',
                    np.reshape(hessvec, (-1, nat3)),
                    np.sqrt(atmas)
                ))
            else:
                raise QuantityError('Missing atomic masses to correct evec')
        else:
            if natoms is not None:
                eigvec = np.reshape(hessvec, (-1, 3*natoms))
            else:
                # Compute nat3 based on: 3*nat*(3nat-ntrro) = size(hessvec)
                # ntrro = 6 for non-linear, 5 otherwise
                # The positive root should be last one (ascending order)
                # assume first most common case: non-linear
                N = hessvec.size
                val = np.polynomial.polynomial.polyroots((-N, -6, 1))[-1]
                if val.is_integer():
                    nat3 = int(val)
                else:
                    nat3 = int(np.polynomial.polynomial.polyroots(
                        (-N, -5, 1))[-1])
                eigvec = np.reshape(hessvec, (-1, nat3))
        return eigvec

    def norm_evec(eigvec):
        """Normalize eigenvector (assumed to have shape (nvib, nat3))."""
        res = np.empty(eigvec.shape)
        norms = np.sum(eigvec**2, axis=1)
        for i in range(eigvec.shape[0]):
            if norms[i] > sys.float_info.epsilon:
                res[i, :] = eigvec[i, :] / sqrt(norms[i])
            else:
                res[i, :] = 0.0
        return res

    if not (get_evec or get_eval):
        raise ValueError('Nothing to do')

    natoms = None
    nvib = None
    atmas = None
    hessvec = None
    hessval = None
    fccart = None
    key_ffx = None
    key_evec = None
    eigvec = None
    eigval = None
    do_calc = force_calc
    if pre_data is not None:
        for key in pre_data:
            qlabel = pre_data[key].qlabel
            if qlabel.label == 'natoms':
                natoms = pre_data['natoms'].data
            elif qlabel.label == 'nvib':
                natoms = pre_data['nvib'].data
            elif qlabel.label == 'atmas':
                # np.repeat used to duplicate atmas for each Cart. coord.
                atmas = np.repeat(np.array(pre_data[key].data), 3)
            elif qlabel.label == 'hessvec':
                key_evec = key
                hessvec = np.array(pre_data[key].data)
            elif qlabel.label == 'hessdat' and qlabel.kind == 'freq':
                hessval = np.array(pre_data[key].data)
            elif (qlabel.label == 1 and qlabel.derord == 2
                  and qlabel.dercrd == 'X'):
                key_ffx = key
                fccart = np.array(pre_data[key].data)

    if force_calc and fccart is not None and atmas is not None:
        if (len(fccart.shape) == 1
                or pre_data[key_ffx]['shape'].lower() == 'lt'):
            fccart = square_ltmat(pre_data[key_ffx].data)
        eigvec, eigval = build_evec(fccart, atmas, get_evec, get_eval, nvib)
    elif not force_calc:
        if get_evec and hessvec is not None:
            # Check if eigenvectors need to be corrected
            #  Default is assumed to be Gaussian fchk
            evec_form = pre_data[key_evec].get('form', 'L.M^{-1/2}')
            eigvec = convert_evec(hessvec, atmas, natoms, evec_form)
        if get_eval and hessval is not None:
            eigval = hessval

    calc_evec = get_evec and eigvec is None
    calc_eval = get_eval and eigval is None
    if calc_evec or calc_eval:
        # We are missing data, now extracting data and recomputing.
        read_data = {}
        if calc_evec or calc_eval:
            read_data['natoms'] = QLabel(quantity='natoms')
            read_data['d2EdX2'] = QLabel(quantity=1, derorder=2, dercoord='X')
            read_data['atmas'] = QLabel(quantity='atmas')
            read_data['nvib'] = QLabel(quantity='nvib')
        if calc_evec:
            read_data['hessvec'] = QLabel(quantity='hessvec', level='H')
        if calc_eval:
            read_data['hessval'] = QLabel(quantity='hessdat',
                                          descriptor='freq', level='H')
        tmp_data = get_data(dfobj, error_noqty=False, **read_data)

        do_calc = force_calc
        if not force_calc:
            if calc_evec:
                if (tmp_data['hessvec'] is not None
                        and tmp_data['atmas'] is not None):
                    hessvec = tmp_data['hessvec'].data
                    atmas = np.repeat(np.array(tmp_data['atmas'].data), 3)
                    natoms = tmp_data['natoms'].data
                    evec_form = tmp_data['hessvec'].get('form', 'L.M^{-1/2}')
                    eigvec = convert_evec(hessvec, atmas, natoms, evec_form)
                else:
                    do_calc = True
            if calc_eval:
                if tmp_data['hessval'] is not None:
                    eigval = np.array(tmp_data['hessval'].data)
                else:
                    do_calc = True

        if do_calc:
            if (tmp_data['d2EdX2'] is not None
                    and tmp_data['atmas'] is not None):
                # Rediagonlizing is more accurate, but require being
                atmas = np.repeat(np.array(tmp_data['atmas'].data), 3)
                ffx = square_ltmat(tmp_data['d2EdX2'].data)
                nvib = None if tmp_data['nvib'] is None \
                    else tmp_data['nvib'].data
                eigvec, eigval = build_evec(ffx, atmas, get_evec, get_eval,
                                            nvib)
            else:
                if calc_evec:
                    msg = 'Unable to retrieve force constants eigenvectors'
                    raise QuantityError(msg)
                if calc_eval:
                    msg = 'Unable to retrieve normal mode wavenumbers'
                    raise QuantityError(msg)

    return eigvec, eigval
