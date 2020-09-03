"""Low-level operations on Gaussian output (log) files.

Provides low-level interfaces to manipulate/extract data in Gaussian
  output files.

Methods
-------
get_data
    Gets data from a GLog file for each quantity label.

Classes
-------
GLogIO
    Main class to handle Gaussian output file operations.
"""

import os  # Used for file existence check
import re  # Used to find keys in log file
import typing as tp
import estampes.parser as ep
from estampes.base import ParseKeyError, QuantityError, TypeData, TypeDCrd, \
    TypeDGLog, TypeDOrd, TypeQInfo, TypeQOpt, TypeQTag, TypeRSta
from estampes.data.physics import PHYSFACT


# ================
# Module Constants
# ================

__ang2au = 1.0 / PHYSFACT.bohr2ang

__tpStrInt = tp.TypeVar('_tp_StrInt', str, int)
# TypeSBloc = tp.Optional[tp.Tuple[str, int]]
# TypeQInfos = tp.Tuple[list, tp.List[str, int, TypeSBloc]]
TypeQData = tp.Dict[str, tp.Optional[tp.Any]]
TypeQKwrd = tp.Tuple[
    tp.Union[int, tp.List[int]],  # Link
    tp.Union[str, tp.List[str]],  # Keyword
    tp.Union[__tpStrInt, tp.List[__tpStrInt]],  # Jump/Skip function
    tp.Union[str, tp.List[str]],  # Matching pattern for data to extract
    #  Block end condition
    tp.Union[tp.Callable[[str], bool], tp.List[tp.Callable[[str], bool]]],
    tp.Union[int, tp.List[int]]  # Number of occurrences
]
TypeKData = tp.Tuple[
    str,  # Keyword
    int,  # Link
    __tpStrInt,  # Information on lines to skip after keyword
    tp.Pattern,  # Data extraction matching pattern (compiled)
    int,  # which occurrences to extract
    tp.Callable[[str], bool]
]


# ==============
# Module Classes
# ==============

class GLogIO(object):
    """Main class to handle Gaussian output file operations.

    Attributes
    ----------
    filename : str
        Gaussian output filename.
    version : str
        Version, software-dependent.
    full_version : tuple
        full version:
        * Gaussian
        * Gaussian major and minor revisions, mach and relesase date

    Methods
    -------
    read_data
        Extracts 1 or more data blocks from Gaussian's log file.
    """
    def __init__(self, fname: str,
                 load_pos: bool = True) -> None:
        self.filename = fname
        self.__linkpos = {}
        self.__route = None
        self.__gversion = None
        self.__rte_opt = None
        self.__links = None
        self.get_head()
        if load_pos:
            self.__store_linkpos()
        # try:
        #     txt = 'Gaussian Version'
        #     self.__gversion = self.get_data(txt)[txt]
        # except ParseKeyError:
        #     self.__gversion = None

    @property
    def filename(self) -> str:
        """Gets or sets the filename associated to the GLog object."""
        return self.__fname

    @filename.setter
    def filename(self, name: str) -> None:
        if not os.path.exists(name):
            raise FileNotFoundError('Formatted checkpoint not found')
        self.__fname = name

    @property
    def version(self) -> tp.Dict[str, str]:
        """Returns the version of Gaussian used to generate the log file.
        """
        return {key: self.__gversion[key] for key in ('major', 'minor')}

    @property
    def full_version(self) -> tp.Tuple[str, tp.Any]:
        """Returns the full version, for the parser interface"""
        return "Gaussian", self.__gversion

    def get_head(self):
        """Returns the header information: Version, Route."""
        keydata = []
        qtydata = {}
        key2blk = {}
        i = 0
        for item in ('route', 'swopt', 'swver'):
            qtydata[item] = (item, None, None, None, None)
            key2blk[item] = (i, i)
            i += 1
            link, key, skips, fmt, end, num = qlab_to_linkdata(item)
            keydata.append((key, link, skips, re.compile(fmt), num, end))
        ndata, data = self.read_data(*keydata)
        data = parse_data(qtydata, key2blk, ndata, data)
        self.__route = data['route']
        self.__links = sorted(set([int(item[0]) for item in self.__route]))
        self.__gversion = data['swver']
        self.__rte_opt = data['swopt']
        # Look at verbosity level of output
        i = self.__rte_opt.index('#') + 1
        if i >= len(self.__rte_opt):
            self.__verb = 0
        else:
            key = self.__rte_opt[i].upper()
            if key == 'P':
                self.__verb = 1
            elif key == 'T':
                self.__verb = -1
            else:
                self.__verb = 0

    def read_data(self,
                  *to_find: TypeKData,
                  raise_error: bool = True) -> TypeDGLog:
        """Extracts data corresponding to the keys to find.

        Parameters
        ----------
        to_find
            List of tuples with the following data:
            keyword: str
                keyword to search.
            link: int
                link where keyword should be found (0 if no specific).
            skip: str/int
                lines to skip from the keyword to reach actual data.
            pattern: obj:`re.Pattern`
                Regular expression pattern object.
            niter: int
                Which occurrences of the quantity to extract.
            endcond: function
                End condition function, which takes a string as argument.
        raise_error
            Only raises error if `True`, otherwise proceeds silently.

        Raises
        ------
        ParseKeyError
            Key not found.

        Notes
        -----
        * The system treats each item in to_find separately.
          Post-processing routines should take care of aliases.
        """
        n_tofind = len(to_find)
        keylist = []  # List of keywords to search
        keydata = []  # Data associated to each keyword
        lnklist = []  # List of links involved
        lnkdata = {}
        datlist = [[] for _ in range(n_tofind)]  # Data to return
        ndatblk = [0 for _ in range(n_tofind)]

        # Generate list of links and check if jump fast search possible
        fast_mode = True
        for i in range(n_tofind):
            link = abs(to_find[i][1])
            new = link not in lnklist
            if new:
                if link == 0:
                    fast_mode = False
                    lnklist.append(link)
                else:
                    if self.__links is None or link in self.__links:
                        lnklist.append(link)
                    if link not in self.__linkpos:
                        fast_mode = False
        lnklist.sort()

        ind = 0
        for link in lnklist:
            imin = ind
            if fast_mode:
                lnkdata[link] = [(imin, 0)]
            for i, item in enumerate(to_find):
                if abs(item[1]) == link:
                    key = item[0]
                    if key not in keylist:
                        keylist.append(key)
                        keydata.append([])
                        j = ind
                        ind += 1
                    else:
                        j = keylist.index(key)
                    keydata[j].append((i, *item[2:]))
            if fast_mode:
                lnkdata[link][1] = ind

        # Sequential Search
        # -----------------
        # Looks for keywords sequentially while reading file
        if not fast_mode:
            block2id = []  # stores real indexes in keylist/keydata
            blockskp = []  # stores the "skip" information
            blockfmt = []  # stores the formats
            blockend = []  # stores the end conditions
            nocc = []  # number of occurrences to extract, used to drop search
            dataid = []  # stores the indexes for the data list
            with open(self.filename, 'r') as fobj:
                # line = fobj.readline()
                for line in fobj:
                    for i, kword in enumerate(keylist):
                        if line.startswith(kword):
                            for j, block in enumerate(keydata[i]):
                                # Save data to correct block in keydata
                                block2id.append([i, j])
                                dataid.append(block[0])
                                blockskp.append(block[1])
                                blockfmt.append(block[2])
                                nocc.append(block[3])
                                blockend.append(block[4])
                                if nocc[-1] > 0:
                                    datlist[dataid[-1]].append([])
                                    ndatblk[dataid[-1]] += 1
                                else:
                                    datlist[dataid[-1]] = []
                                    ndatblk[dataid[-1]] = 1

                    if block2id:
                        lblock = len(block2id)
                        iblock = 0
                        while iblock < lblock:
                            if isinstance(blockskp[iblock], str):
                                if line.startswith(blockskp[iblock]):
                                    blockskp[iblock] = 0
                            if blockskp[iblock] == 0:
                                if blockfmt[iblock].match(line):
                                    res = blockfmt[iblock].match(line)
                                    if nocc[iblock] > 0:
                                        datlist[dataid[iblock]][-1].append(
                                            res.groupdict()['val'])
                                    else:
                                        datlist[dataid[iblock]].append(
                                            res.groupdict()['val'])
                                if blockend[iblock](line):
                                    if nocc[iblock] == 0:
                                        i, j = block2id[iblock]
                                        del keydata[i][j]
                                        if not keydata[i]:
                                            del keydata[i]
                                            del keylist[i]
                                        for k in range(lblock):
                                            a, b = block2id[k]
                                            if a > i:
                                                block2id[k][0] -= 1
                                            elif a == i and b > j:
                                                block2id[k][1] -= 1
                                    del block2id[iblock]
                                    del dataid[iblock]
                                    del blockskp[iblock]
                                    del blockfmt[iblock]
                                    del nocc[iblock]
                                    del blockend[iblock]
                                    lblock -= 1
                                else:
                                    iblock += 1
                            else:
                                if isinstance(blockskp[iblock], int):
                                    if blockskp[iblock] > 0:
                                        blockskp[iblock] -= 1
                                iblock += 1
                    if not keylist:
                        break

        # Fast Search
        # -----------
        else:
            raise NotImplementedError('Fast search not yet ready')

        return ndatblk, datlist

    def __store_linkpos(self):
        """Stores the link header positions in the file if available.

        Loads the keys present in the file and pointers to their
          position to speed up their search.
        Data type and block information are also stored.
        """
        link_heads = {
            # 1: 'Entering Gaussian System,',
            1: 'Entering Link 1,',
            601: 'Population analysis using the SCF Density.',
            716: 'Full mass-weighted force constant matrix:',
            717: 'Second-order Perturbative Anharmonic Analysis',
            718: 'Generation of the Franck-Condon spectrum'
        }
        # to_search = re.compile(r'''\
        #     (?P<title>[\w\s]+?)\s*  # Key
        #     \b(?P<type>[IRC])\b\s*  # Data type
        #     (?P<block>N=)?\s+  # N= only set for non-scalar data
        #     (?P<value>[\d\-\+\.E]+)  # Block size (N=) or scalar value
        #     \n''', re.VERBOSE)

        # keys = {}
        # with open(self.filename, 'r') as fobj:
        #     fpos = 0
        #     for line in fobj:
        #         res = to_search.match(line)
        #         if res:
        #             keys[res.group(1)] = (
        #                 res.group(2),
        #                 int(res.group(3) and res.group(4) or 0),
        #                 fpos)
        #         fpos += len(line)
        # return keys

# ================
# Module Functions
# ================


def qlab_to_linkdata(qtag: TypeQTag,
                     qopt: TypeQOpt = None,
                     dord: TypeDOrd = None,
                     dcrd: TypeDCrd = None,
                     rsta: TypeRSta = None,
                     gver: tp.Optional[str] = None) -> TypeQKwrd:
    """Returns the keyword(s) relevant for a given quantity.

    Returns the a tuple, containing:
    1. the link(s), which may provide the quantities.
        The returned value can be:
        * `=0`: No specific link provides the data
        * `>0`: Link which produces the data
    2. specific keyword to search, *from the start of the line*
    3. jump/skip information:
        * `int`: Number of lines to pass.
        * `str`: Sub-string to search (ex: '----').
    4. Regular expression for the data specification.
        The key `val` refers to the data of interest.
    5. Ending condition, as a function taking a string in arg.
    6. Number of occurrences to search, as integer:
        * `0`: only the first occurrence, then stop.
        * `-1`: takes the last, discarding all previous ones.
        * `1`: takes all occurrences.
    If multiple links can provide the data or alternative keywords are
      available, the function returns tuples of data for each quantity.


    Parameters
    ----------
    qtag
        Quantity identifier or label.
    qopt
        Quantity-specific options.
    dord
        Derivative order.
    dcrd
        Reference coordinates for the derivatives.
    rsta
        Reference state or transition:
        * scalar: reference state
        * tuple: transition
    gver
        Gaussian version.

    Returns
    -------
    int or list of int
        Link(s), which may provide the quantities.
        The returned value can be:
        * `=0`: No specific link provides the data.
        * `>0`: Link which produces the data.
        * `<0`: Data may not be always printed.
    str or list of str
        Specific keyword to search.
    str or int or list of str or list of int
        Jump information:
        * `int`: number of lines to pass
        * `str`: sub-string to search (ex: '----')
    str or list of str
        Regular expression for the data specification.
    function or list of function
        End condition, as a function taking a `str` as argument.
    int or list of int
        Number of occurrences to extract. Possible values:
        * `0`: only takes the first occurrence.
        * `-1`: only takes the last occurrence.
        * `1`: takes all possible occurrences.

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
    lnk = None
    key = None
    sub = None
    fmt = None
    num = None
    if qtag == 'route':  # Log specific extraction, not provided in base
        lnk = 1
        key = ' Cite this work as:'
        sub = ' -'
        def end(s): return s.startswith(' Charge =')
        fmt = r'^ (?P<val>\d+\/(?:,?\d+=-?\d+)*\/(?:,?\d+)+(?:\(-\d\))?);\s*$'
        num = 0
    elif qtag == 'natoms':
        lnk = 101
        key = ' NAtoms='
        sub = 0
        def end(s): return True
        fmt = r'^ NAtoms=\s+(?P<val>\d+) NActive=\s+.*$'
        num = 0
    elif qtag == 'nvib':
        raise NotImplementedError()
    elif qtag == 'atmas':
        lnk = (-101, 716)
        key = (' AtmWgt= ', '- Thermochemistry -')
        sub = (0, 3)
        end = (lambda s: s.startswith(' Leave Link'),
               lambda s: s.startswith(' Molecular Mass:'))
        fmt = (r'^ AtmWgt=\s+(?P<val>(\s+\d+\.\d+)+)\s*$',
               r'^ Atom\s+\d+ has atomic number\s+\d+ and '
               + r'mass\s+(?P<val>\d+\.\d+)\s*$')
        num = 0
    elif qtag == 'atnum':
        lnk = (0, 0)
        key = ('                         Standard orientation:',
               '                          Input orientation:')
        sub = (5, 5)
        end = (lambda s: s.startswith(' ------'),
               lambda s: s.startswith(' ------'))
        txt = r'^\s+\d+\s+(?P<val>\d+)\s+\d+(?:\s+-?\d+\.\d+){3}\s*$'
        fmt = (txt, txt)
        num = (0, 0)
    elif qtag == 'molsym':
        lnk = 101
        key = ' Framework group '
        sub = 0
        def end(s): return True
        fmt = r'^ Framework group \s+(?P<val>[^ ]+)\s*$'
        num = 0
    elif qtag == 'atcrd' or qtag == 2:
        lnk = (0, 0)
        key = ('                         Standard orientation:',
               '                          Input orientation:')
        sub = (5, 5)
        end = (lambda s: s.startswith(' ------'),
               lambda s: s.startswith(' ------'))
        txt = r'^(?:\s+\d+){3}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
        fmt = (txt, txt)
        if qopt == 'first':
            num = (0, 0)
        elif qopt == 'last':
            num = (-1, -1)
        else:
            num = (1, 1)
    elif qtag in ('hessvec', 'hessval'):
        lnk = 716
        key = ' and normal coordinates:'
        sub = 1
        def end(s): return s.startswith(' - Thermochemistry')
        fmt = NotImplemented
        raise NotImplementedError()
    elif qtag == 'swopt':
        lnk = 1
        key = ' Cite this work as:'
        sub = ' #'
        def end(s): return s.startswith(' -')
        fmt = r'^ (?P<val>.*)$'
        num = 0
    elif qtag == 'swver':
        lnk = 1
        key = ' ****'
        sub = 1
        def end(s): return True
        fmt = r'^ (?P<val>Gaussian (?:\d\d|DV):\s.*)$'
        num = 0
    elif qtag == 'fcdat':
        if qopt == 'JMat':
            lnk = (-718, -718)
            key = (' Duschinsky matrix', ' Final Duschinsky matrix')
            sub = (2, 2)
            end = (lambda s: not s.strip(),
                   lambda s: not s.strip())
            txt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            fmt = (txt, txt)
            num = (-1, -1)
        elif qopt == 'JMatF':
            lnk = -718
            key = ' Full Duschinsky matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qopt == 'KVec':
            lnk = -718
            key = ' Shift Vector'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$'
            num = 0
        elif qopt == 'SRAMat':
            lnk = -718
            key = ' A Matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qopt == 'SRBVec':
            lnk = -718
            key = ' B Vector'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$'
            num = 0
        elif qopt == 'SRCMat':
            lnk = -718
            key = ' C Matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qopt == 'SRDVec':
            lnk = -718
            key = ' D Vector'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$'
            num = 0
        elif qopt == 'SREMat':
            lnk = -718
            key = ' E Matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qopt == 'Spec':
            lnk = -718
            key = '                       Final Spectrum'
            sub = ' -----------'
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:-?\d+.\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3})+)\s*$'
            num = 0
        elif qopt == 'SpcLeg':
            lnk = 718
            key = '                       Final Spectrum'
            sub = ' ------'
            def end(s): return s.startswith(' -----------')
            fmt = r'^\s+(?P<val>.*: .*)\s*$'
            num = 0
        elif qopt == 'BShape':
            lnk = 718
            key = '                       Final Spectrum'
            sub = 3
            def end(s): return not s.strip()
            # sub = 2
            # def end(s): return s.startswith(' Legend')
            fmt = r'^\s+(?P<val>.*\w.*)\s*$'  # The \w is to exclude empty lines
            num = 0
        elif qopt == 'Conv':
            # 2 sets, one for the main one, one to extract the FCF (HT)
            lnk = (-718, -718)
            key = ('              Calculations of Band Intensities',
                   '              Calculations of Band Intensities')
            sub = (2, 2)
            end = (lambda s: s.startswith('     ===='),
                   lambda s: s.startswith('     ===='))
            fmt = (r'^\s+Spectrum progression:\s+(?P<val>-?\d+\.\d+)%\s*$',
                   r'^\s+.+Franck-Condon Factors:\s+(?P<val>-?\d+\.\d+)%)\s*$')
            num = (0, 0)
        elif qopt == 'Assign':
            lnk = (-718, -718)
            key = ('                 Information on Transitions',
                   '                 Information on Transitions')
            sub = (2, 2)
            end = (lambda s: s.startswith('     ===='),
                   lambda s: s.startswith('     ===='))
            fmt = (r'^\s+Energy =\s+(?P<val>-?\d+\.\d+ cm.-1: .*)\s*$',
                   r'^\s+-. Intensity =\s+(?P<val>.*)\s*$')
            num = (0, 0)
        elif qopt == 'GeomIS':
            lnk = 718
            key = '              New orientation in initial state'
            sub = 5
            def end(s): return s.startswith(' ------')
            fmt = r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
        elif qopt == 'GeomFS':
            lnk = 718
            key = '              New orientation in final state'
            sub = 5
            def end(s): return s.startswith(' ------')
            fmt = r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
        elif qopt == 'GeomMS':
            lnk = 718
            key = '              New orientation in intermediate state'
            sub = 5
            def end(s): return s.startswith(' ------')
            fmt = r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
        elif qopt == 'ExGeom':
            lnk = -718
            key = ' Extrapolated geometry'
            sub = 4
            def end(s): return s.startswith(' ------')
            fmt = r'^\s+\w+(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
    elif qtag == 'vptdat':
        raise NotImplementedError()
    elif qtag in ('dipstr', 'rotstr'):
        raise NotImplementedError()
    else:
        raise NotImplementedError()
        # if type(rsta) is tuple:
        #     raise NotImplementedError()
        #     if qtag == 1 and dord == 0:
        #         keywords = ['ETran scalars', 'SCF Energy']
        # else:
        #     raise NotImplementedError()
        #     if qtag == 1:
        #         if dord == 0:
        #             if rsta == 'c':
        #                 del keywords[:]
        #                 raise NotImplementedError()
        #             elif type(rsta) is int:
        #                 if rsta == 0:
        #                     keyword = 'SCF Energy'
        #                 else:
        #                     keyword = 'ETran state values'
        #                 keywords.append('Total Energy', 'ETran scalars')
        #         elif dord == 1:
        #             if dcrd is None or dcrd == 'X':
        #                 if rsta == 'c' or type(rsta) is int:
        #                     keyword = 'Cartesian Gradient'
        #         elif dord == 2:
        #             if dcrd is None or dcrd == 'X':
        #                 if rsta == 'c' or type(rsta) is int:
        #                     keyword = 'Cartesian Force Constants'
        #     elif qtag == 50:
        #         raise NotImplementedError()
        #     elif qtag == 91:
        #         raise NotImplementedError()
        #     elif qtag == 92:
        #         raise NotImplementedError()
        #     elif qtag == 93:
        #         raise NotImplementedError()
        #     elif qtag == 101:
        #         raise NotImplementedError()
        #     elif qtag == 102:
        #         raise NotImplementedError()
        #     elif qtag == 103:
        #         raise NotImplementedError()
        #     elif qtag == 104:
        #         raise NotImplementedError()
        #     elif qtag == 105:
        #         raise NotImplementedError()
        #     elif qtag == 106:
        #         raise NotImplementedError()
        #     elif qtag == 107:
        #         raise NotImplementedError()
        #     elif qtag == 201:
        #         raise NotImplementedError()
        #     elif qtag == 202:
        #         raise NotImplementedError()
        #     elif qtag == 203:
        #         raise NotImplementedError()
        #     elif qtag == 204:
        #         raise NotImplementedError()
        #     elif qtag == 205:
        #         raise NotImplementedError()
        #     elif qtag == 206:
        #         raise NotImplementedError()
        #     elif qtag == 207:
        #         raise NotImplementedError()
        #     elif qtag == 208:
        #         raise NotImplementedError()
        #     elif qtag == 209:
        #         raise NotImplementedError()
        #     elif qtag == 300:
        #         raise NotImplementedError()
        #     elif qtag == 301:
        #         raise NotImplementedError()
        #     elif qtag == 302:
        #         raise NotImplementedError()
        #     elif qtag == 303:
        #         raise NotImplementedError()
        #     elif qtag == 304:
        #         raise NotImplementedError()
        #     elif qtag == 305:
        #         raise NotImplementedError()
        #     elif qtag == 306:
        #         raise NotImplementedError()
        #     else:
        #         raise QuantityError('Unknown quantity')

    return lnk, key, sub, fmt, end, num


def parse_data(qdict: TypeQInfo,
               key2blocks: tp.Dict[str, tp.Tuple[int, int]],
               ndatablock: tp.Sequence[int],
               datablocks: TypeDGLog,
               gver: tp.Optional[tp.Tuple[str, str]] = None,
               raise_error: bool = True) -> TypeQData:
    """Parses data arrays to extract specific quantity.

    Parses data array to extract relevant information for each quantity.

    Parameters
    ----------
    qdict
        Dictionary of quantities.
    key2blocks
        Range tuples associating each keyword to the data blocks.
    ndatablock
        Number of occurrences of the data blocks extracted.
    datablocks
        Data blocks.  Each data block may contain several blocks if
          multiple occurrences have been extracted.
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
    QuantityError
        Unsupported quantity.
    IndexError
        Unexpected data keyword for dictionary-type data.
    """
    def empty_cases_ok(qtag, qopt):
        return qtag == 'nvib' or \
            (qtag == 'fcdat' and qopt in ('JMat', 'JMatF'))
    data = {}
    for qlabel in qdict:
        qtag, qopt, dord, dcrd, rsta = qdict[qlabel]
        first, last = key2blocks[qlabel]
        # Basic Check: property available
        # -----------
        # Check if some data extracted
        num = 0
        if first == last:
            iref = first
            num = len(datablocks[iref])
        else:
            iref = -1
            for i in range(first, last+1):
                if num == 0 and ndatablock[i] > 0 and datablocks[i]:
                    iref = i  # Ignore first, null indexes
                    num = len(datablocks[i])
        if num == 0 and not empty_cases_ok(qtag, qopt):
            if raise_error:
                raise ParseKeyError('Missing quantity in file')
            else:
                data[qlabel] = None
                continue
        # Basic Properties/Quantities
        # ---------------------------
        if qtag == 'route':
            fmt = '{:d}{:02d}'
            data[qlabel] = []
            for num, line in enumerate(datablocks[iref]):
                ov, iops, links = line.split('/')
                if '(' in links:
                    links, jump = links[:-1].split('(')
                    toline = num + 1 + int(jump)
                else:
                    toline = 0
                iops = tuple(iops.split(','))
                for link in links.split(','):
                    fulllink = fmt.format(int(ov), int(link))
                    data[qlabel].append((fulllink, num+1, iops, toline))
        elif qtag == 'natoms':
            data[qlabel] = int(datablocks[iref][0])
        elif qtag == 'atcrd' or qtag == 2:
            data[qlabel] = []
            # By default, we choose the standard orientation if present
            if datablocks[last]:
                i = last
            else:
                i = first
            for line in datablocks[i]:
                data[qlabel].append([float(item)*__ang2au
                                     for item in line.split()])
        elif qtag == 'atmas':
            data[qlabel] = []
            for line in datablocks[iref]:
                data[qlabel].extend(float(item) for item in line.split())
        elif qtag in ('atnum',):
            data[qlabel] = []
            for line in datablocks[iref]:
                data[qlabel].extend(int(item) for item in line.split())
        elif qtag == 'swopt':
            # The last line are the final dashes
            data[qlabel] = ' '.join(datablocks[iref][:-1])
        elif qtag == 'molsym':
            raise NotImplementedError()
        elif qtag == 'swver':
            txt = r'\s*Gaussian (\w+):\s+(\w+)-(\w{3})Rev([\w.+]+) ' \
                + r'(\d+-\w{3}-\d{4})\s*'
            pattern = re.compile(txt)
            res = re.match(pattern, ''.join(datablocks[iref])).groups()
            data[qlabel] = {'major': res[2], 'minor': res[3],
                            'system': res[1], 'release': res[4]}
        # Vibrational Information
        # -----------------------
        # Technically state should be checked but considered irrelevant.
        elif qtag == 'nvib':
            if iref >= 0:
                data[qlabel] = int(datablocks[iref])
            else:
                # For a robust def of nvib, we need the symmetry and
                #   the number of frozen atoms. For now, difficult to do.
                raise NotImplementedError()
        elif qtag in ('hessvec', 'hessval'):
            if iref >= 0:
                data[qlabel] = datablocks[iref]
        # Vibronic Information
        # --------------------
        elif qtag == 'fcdat':
            if qopt == 'JMat':
                if datablocks[last]:
                    i = last
                else:
                    i = first
                if 'diagonal' in datablocks[i]:
                    data[qlabel] = []
                else:
                    data[qlabel] = []
                    N = 0
                    for line in datablocks[i]:
                        cols = line.split()
                        irow = int(cols[0]) - 1
                        if irow == N:
                            data[qlabel].append([])
                            N += 1
                        data[qlabel][irow].extend(
                            [float(item.replace('D', 'e'))
                             for item in cols[1:]])
            elif qopt == 'JMatF':
                if 'diagonal' in datablocks[iref]:
                    data[qlabel] = []
                else:
                    data[qlabel] = []
                    N = 0
                    for line in datablocks[iref]:
                        cols = line.split()
                        irow = int(cols[0]) - 1
                        if irow == N:
                            data[qlabel].append([])
                            N += 1
                        data[qlabel][irow].extend(
                            [float(item.replace('D', 'e'))
                             for item in cols[1:]])
            elif qopt == 'KVec':
                data[qlabel] = []
                for line in datablocks[iref]:
                    data[qlabel].append(float(line.replace('D', 'e')))
            elif qopt == 'SRAMat':
                data[qlabel] = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data[qlabel].append([])
                        N += 1
                    data[qlabel][irow].extend(
                        [float(item.replace('D', 'e'))
                            for item in cols[1:]])
            elif qopt == 'SRBVec':
                data[qlabel] = []
                for line in datablocks[iref]:
                    data[qlabel].append(float(line.replace('D', 'e')))
            elif qopt == 'SRCMat':
                data[qlabel] = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data[qlabel].append([])
                        N += 1
                    data[qlabel][irow].extend(
                        [float(item.replace('D', 'e'))
                            for item in cols[1:]])
            elif qopt == 'SRDVec':
                data[qlabel] = []
                for line in datablocks[iref]:
                    data[qlabel].append(float(line.replace('D', 'e')))
            elif qopt == 'SREMat':
                data[qlabel] = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data[qlabel].append([])
                        N += 1
                    data[qlabel][irow].extend(
                        [float(item.replace('D', 'e'))
                            for item in cols[1:]])
            elif qopt == 'Spec':
                data[qlabel] = {'x': []}
                ny = len(datablocks[iref][0].split()) - 1
                if ny == 1:
                    yfmt = 'y'
                    data[qlabel]['y'] = []
                else:
                    yfmt = 'y{{idy:0{}d}}'.format(ny)
                    for i in range(ny):
                        data[qlabel][yfmt.format(idy=i+1)] = []
                for line in datablocks[iref]:
                    cols = [float(item.replace('D', 'e'))
                            for item in line.split()]
                    data[qlabel]['x'].append(cols[0])
                    for i, item in enumerate(cols[1:]):
                        data[qlabel][yfmt.format(idy=i+1)].append(item)
            elif qopt == 'SpcLeg':
                data[qlabel] = {}
                if len(datablocks[iref]) <= 3:
                    # X, Y, Int, don't use numbering
                    yfmt = 'y'
                else:
                    yfmt = 'y{{idy:0{}d}}'.format(
                        len(str(len(datablocks[iref])-2)))
                for line in datablocks[iref]:
                    key, title = line.split(':')
                    if key.strip() == '1st col.':
                        _key = 'x'
                    elif key.strip() == 'Intensity':
                        _key = 'I'
                    else:
                        try:
                            idy = int(key.strip()[0]) - 1
                        except ValueError:
                            raise IndexError('Unrecognized key in spc leg.: ' +
                                             key)
                        _key = yfmt.format(idy=idy)
                    data[qlabel][_key] = title.strip()
            elif qopt == 'BShape':
                if len(datablocks[iref]) == 1:
                    if 'No' in datablocks[iref][0]:
                        data[qlabel] = {'func': 'stick', 'HWHM': None}
                    else:
                        raise IndexError('Unrecognized broadening function.')
                elif len(datablocks[iref]) == 2:
                    func = datablocks[iref][0].split()[-3].lower()
                    hwhm = float(datablocks[iref][1].split()[-2])
                    data[qlabel] = {'func': func, 'HWHM': hwhm}
                else:
                    raise IndexError('Unrecognized broadening definition.')
            elif qopt == 'Conv':
                raise NotImplementedError()
            elif qopt == 'Assign':
                raise NotImplementedError()
            elif qopt == 'GeomIS':
                data[qlabel] = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data[qlabel].append([float(item)*__ang2au
                                         for item in line.split()])
            elif qopt == 'GeomFS':
                data[qlabel] = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data[qlabel].append([float(item)*__ang2au
                                         for item in line.split()])
            elif qopt == 'ExGeom':
                data[qlabel] = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data[qlabel].append([float(item)*__ang2au
                                         for item in line.split()])
        # Anharmonic Information
        # ----------------------
        elif qtag == 'vptdat':
            raise NotImplementedError()
        # State(s)-dependent quantities
        # -----------------------------
        else:
            raise NotImplementedError()
            # # Transition moments
            # # ^^^^^^^^^^^^^^^^^^
            # if type(rsta) is tuple:
            #     data[qlabel] = _parse_electrans_data(qtag, datablocks, kword,
            #                                          qopt, dord, dcrd, rsta)
            # # States-specific Quantities
            # # ^^^^^^^^^^^^^^^^^^^^^^^^^^
            # else:
            #     key = 'ETran scalars'
            #     if key in datablocks:
            #         (nstates, ndata, _, _, iroot,
            #             _) = [item for item in datablocks[key][:6]]
            #     curr_sta = rsta == 'c' or rsta == iroot
            #     # Only energy is currently computed for all states:
            #     if rsta == 'a' and qtag == 2:
            #         data = [float(datablocks[kword][i*ndata])
            #                 for i in range(nstates)]
            #     # Data for current electronic states
            #     elif curr_sta:
            #         if qtag in ('dipstr', 'rotstr'):
            #             if qopt == 'H':
            #                 key = 'Number of Normal Modes'
            #             else:
            #                 key = 'Anharmonic Number of Normal Modes'
            #             if key not in datablocks:
            #                 raise ParseKeyError('Missing required dimension')
            #             ndat = int(datablocks[key])
            #             if qtag == 'dipstr':
            #                 offset = 7*ndat
            #             else:
            #                 offset = 8*ndat
            #             data[qlabel] = datablocks[kword][offset:offset+ndat]
            #         elif qtag == 1:
            #             data[qlabel] = datablocks[kword]
            #         elif qtag == 92:
            #             data[qlabel] = datablocks[kword][:9]
            #         elif qtag == 93:
            #             data[qlabel] = datablocks[kword][9:]
            #         elif qtag in (50, 91):
            #             raise NotImplementedError()
            #         elif qtag == 101:
            #             if dord in (0, 1):
            #                 data[qlabel] = datablocks[kword]
            #             else:
            #                 raise NotImplementedError()
            #         elif qtag == 102:
            #             if dord == 1:
            #                 data[qlabel] = datablocks[kword]
            #             else:
            #                 raise NotImplementedError()
            #         elif qtag == 300:
            #             if dord in (0, 1):
            #                 if qopt == 0:
            #                     data[qlabel] = datablocks[kword]
            #             else:
            #                 raise NotImplementedError()
            #         elif qtag == 300:
            #             if dord in (0, 1):
            #                 data[qlabel] = datablocks[kword]
            #             else:
            #                 raise NotImplementedError()
            #         else:
            #             raise NotImplementedError()

    return data


def get_data(dfobj: GLogIO,
             *qlabels: str,
             error_noqty: bool = True) -> TypeData:
    """Gets data from a GLog file for each quantity label.

    Reads one or more full quantity labels from `qlab` and returns the
      corresponding data.

    Parameters
    ----------
    dfobj
        Gaussian output file as `GLogIO` object.

    *qlabels
        List of full quantity labels to parse.
    error_noqty
        If True, error is raised if the quantity is not found.

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
    if not isinstance(dfobj, GLogIO):
        raise TypeError('GLogIO instance expected')
    # Check if anything to do
    if len(qlabels) == 0:
        return None
    # Build Keyword List
    # ------------------
    # List of keywords
    keydata = []
    key2blocks = {}
    idata = 0
    qty_dict = {}
    for qlabel in qlabels:
        # Label parsing
        # ^^^^^^^^^^^^^
        qty_dict[qlabel] = ep.parse_qlabel(qlabel)
        links, keys, skips, fmts, ends, nums = \
            qlab_to_linkdata(*qty_dict[qlabel])
        if isinstance(links, tuple):
            nblocks = len(links)
            key2blocks[qlabel] = (idata, idata+nblocks-1)
            for i in range(nblocks):
                keydata.append((keys[i], links[i], skips[i],
                                re.compile(fmts[i]), nums[i], ends[i]))
            idata += nblocks
        else:
            key2blocks[qlabel] = (idata, idata)
            keydata.append((keys, links, skips, re.compile(fmts), nums, ends))
            idata += 1
    # Data Extraction
    # ---------------
    # Use of set to remove redundant keywords
    ndata, datablocks = dfobj.read_data(*keydata, raise_error=False)
    # Data Parsing
    # ------------
    gver = (dfobj.version['major'], dfobj.version['minor'])
    try:
        data = parse_data(qty_dict, key2blocks, ndata, datablocks, gver,
                          error_noqty)
    except (QuantityError, NotImplementedError):
        raise QuantityError('Unsupported quantities')
    except (ParseKeyError, IndexError):
        raise IndexError('Missing data in Gaussian log')

    return data