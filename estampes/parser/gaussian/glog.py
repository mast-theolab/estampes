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

from estampes import parser as ep
from estampes.base import ParseKeyError, QuantityError, TypeDCrd, TypeDGLog, \
    TypeDOrd, TypeQData, TypeQInfo, TypeQLvl, TypeQOpt, TypeQTag, TypeRSta
from estampes.data.physics import PHYSFACT, phys_fact


# ================
# Module Constants
# ================

__ang2au = 1.0 / PHYSFACT.bohr2ang

_tp_StrInt = tp.TypeVar('_tp_StrInt', str, int)
# TypeSBloc = tp.Optional[tp.Tuple[str, int]]
# TypeQInfos = tp.Tuple[list, tp.List[str, int, TypeSBloc]]
TypeQKwrd = tp.Tuple[
    tp.Union[int, tp.List[int]],  # Link
    tp.Union[str, tp.List[str]],  # Keyword
    tp.Union[_tp_StrInt, tp.List[_tp_StrInt]],  # Jump/Skip function
    tp.Union[str, tp.List[str]],  # Matching pattern for data to extract
    #  Block end condition
    tp.Union[tp.Callable[[str], bool], tp.List[tp.Callable[[str], bool]]],
    tp.Union[int, tp.List[int]]  # Number of occurrences
]
TypeKData = tp.Tuple[
    str,  # Keyword
    int,  # Link
    _tp_StrInt,  # Information on lines to skip after keyword
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
            qtydata[item] = ep.parse_qlabel(ep.build_qlabel(item))
            key2blk[item] = (i, i)
            i += 1
            link, key, skips, fmt, end, num = qlab_to_linkdata(item)
            keydata.append((key, link, skips, re.compile(fmt), num, end))
        ndata, data = self.read_data(*keydata)
        data = parse_data(qtydata, key2blk, ndata, data)
        self.__route = data['route']['data']
        self.__links = sorted(set([int(item[0]) for item in self.__route]))
        self.__gversion = data['swver']
        self.__rte_opt = data['swopt']['data']
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
        def del_block(iblock: int,
                      block2id: tp.Sequence[tp.Sequence[int]],
                      nocc: tp.Sequence[int],
                      dataid: tp.Sequence[int],
                      blockskp: tp.Sequence[tp.Union[str, int]],
                      blockfmt: tp.Sequence[tp.Pattern],
                      blockend: tp.Callable[[str], bool]) -> int:
            """Deletes a block in the lookup tables.

            Returns
            -------
            int
                Status, as integer
                0: block removed
                1: keylist item removed"""
            istat = 0
            if nocc[iblock] == 0:
                i, j = block2id[iblock]
                del keydata[i][j]
                if not keydata[i]:
                    del keydata[i]
                    del keylist[i]
                    istat = 1
                for k in range(len(block2id)):
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
            return istat

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
                for line in fobj:
                    i = -1
                    for kword in keylist:
                        skip = False
                        i += 1
                        if line.startswith(kword):
                            iblock = 0
                            while iblock < len(block2id):
                                if block2id[iblock][0] == i:
                                    res = del_block(iblock, block2id, nocc,
                                                    dataid, blockskp, blockfmt,
                                                    blockend)
                                    if res == 1:  # keylist empty
                                        i -= 1
                                        skip = True
                                else:
                                    iblock += 1
                            if not skip:
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
                                    res = del_block(iblock, block2id, nocc,
                                                    dataid, blockskp, blockfmt,
                                                    blockend)
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
                     qlvl: TypeQLvl = None,
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
    qlvl
        Level of theory use to generate the quantity.
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
        key = (' NAtoms= ', ' - Thermochemistry -')
        sub = (0, 3)
        end = (lambda s: s.startswith(' Leave Link'),
               lambda s: s.startswith(' Molecular Mass:'))
        fmt = (r'^ AtmWgt=\s+(?P<val>(\s+\d+\.\d+)+)\s*$',
               r'^ Atom\s+\d+ has atomic number\s+\d+ and '
               + r'mass\s+(?P<val>\d+\.\d+)\s*$')
        num = (0, 0)
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
        if qopt == 'SimInf':
            lnk = -718
            key = '               Information on the Simulation'
            sub = 2
            def end(s): return s.startswith('     ====')
            fmt = r'^\s+(?P<val>.*\w.*)\s*$'
            num = 0
        elif qopt == 'JMat':
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
            lnk = (-718, -718)
            key = ('                       Final Spectrum',
                   ' Legend:')
            sub = (' -----------', ' -----------')
            end = (lambda s: not s.strip(),
                   lambda s: 'Legend:' in s or not s.strip())
            fmt = (r'^\s+(?P<val>(?:-?\d+.\d+)'
                   + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3})+)\s*$',
                   r'^\s+(?P<val>(?:-?\d+.\d+)'
                   + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3})+)\s*$')
            num = (0, 1)
        elif qopt == 'SpcPar':
            lnk = (718, 718)
            key = ('                       Final Spectrum',
                   ' Legend:')
            sub = (3, 0)
            end = (lambda s: s.startswith(' -----------'),
                   lambda s: s.startswith(' -----------'))
            fmt = (r'^\s+(?P<val>.*\w.*)\s*$',  # \w to exclude empty lines
                   r'^\s+(?P<val>.*\w.*)\s*$')
            num = (0, 1)
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
            # The second block is for my working
            lnk = (718, 718)
            key = ('              New orientation in final state',
                   '               New orientation in final state')
            sub = (5, 5)
            end = (lambda s: s.startswith(' ------'),
                   lambda s: s.startswith(' ------'))
            fmt = (r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$',
                   r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$')
            num = (0, 0)
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
        if qopt == 'CICoef':
            # The second group is to correct the numbering if passive modes
            #     present.
            lnk = (717, 717)
            key = (' Definition of New States w.r.t. Deperturbed States',
                   ' Reduced-Dimensionality on Variational States')
            sub = (4, 3)
            end = (lambda s: not s.strip(),
                   lambda s: not s.strip())
            fmt = (r'^\s*(?P<val>(\d+\s+:|)\s+[+-]?\d\.\d+\s+x\s+'
                   + r'\|[0-9();]+>)\s*$',
                   r'^\s*(?P<val>\d+\s+\|\s+\d+)\s*$')
            num = (0, 0)
        else:
            raise NotImplementedError()
    elif qtag == 'vtrans':
        if qlvl == 'H':
            lnk = (-716, 716, -717)
            key = (' and normal coordinates:',
                   ' and normal coordinates:',
                   ' NOTE: Transition energies are given with')
            sub = (1, 1, 8)
            end = (lambda s: s.startswith(' - Thermochemistry'),
                   lambda s: s.startswith(' - Thermochemistry'),
                   lambda s: s.startswith('     ====='))
            fmt = (r'^\s{16}(?P<val>(?:\s+\d+){1,5})\s*$',
                   r'^\s{16}(?P<val>(?:\s+\d+){1,3})\s*$',
                   r'^\s+\w?\s+(?P<val>\s*\d+\(\d+\))\s+\w+\s+'
                   + r'(?:\s+-?\d+\.\d+|\*+){4}.*\s*$')
            num = (0, -1, 0)
        elif qlvl == 'A':
            lnk = 717
            key = ' NOTE: Transition energies are given with'
            sub = 8
            def end(s): return s.startswith('     =====')
            fmt = r'^\s+\w?\s+(?P<val>(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'\
                  + r'(?:\w+)?\s+(?:\s*-?\d+\.\d+|\*+\s+){4,5}.*\s*$'
            num = 0
        else:
            raise NotImplementedError()
    elif qtag == 'vlevel':
        if qlvl == 'H':
            lnk = (-716, 716, -717)
            key = (' and normal coordinates:',
                   ' and normal coordinates:',
                   ' NOTE: Transition energies are given with')
            sub = (1, 1, 8)
            end = (lambda s: s.startswith(' - Thermochemistry'),
                   lambda s: s.startswith(' - Thermochemistry'),
                   lambda s: s.startswith('     ====='))
            fmt = (r'^\s+Frequencies --- \s*(?P<val>\d.*)\s*$',
                   r'^\s+Frequencies -- \s*(?P<val>\d.*)\s*$',
                   r'^\s+\w?\s+(?:\s*\d+\(\d+\))\s+\w+\s+'
                   + r'(?P<val>-?\d+\.\d+|\*+)(?:\s+-?\d+\.\d+|\*+){4}.*\s*$')
            num = (0, -1, 0)
        elif qlvl == 'A':
            lnk = 717
            key = ' NOTE: Transition energies are given with'
            sub = 8
            def end(s): return s.startswith('     =====')
            fmt = r'^\s+\w?\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+(?:\w+)?\s+' \
                  + r'(?:-?\d+\.\d+|\*+)?\s+(?P<val>-?\d+\.\d+|\*+)' \
                  + r'(?:\s*-?\d+\.\d+|\*+){3}.*\s*$'
            num = 0
        else:
            raise NotImplementedError()
    else:
        lnk0 = -914
        key0 = ' Excitation energies and oscillator strengths:'
        sub0 = 2
        def end0(s): return s.startswith(' End of Minotr F.D. properties file')
        fmt0 = r'^\s+(?P<val>Excited State\s+\d+: .*|' +\
            r'This state for optimization.*)\s*$'
        num0 = 0
        if type(rsta) is tuple:
            Si, Sf = rsta
            if qtag == 1:
                if dord == 0:
                    if Si == 0:
                        lnk = 0
                        key = ' Excitation energies and oscillator strengths:'
                        sub = 1
                        def end(s): return s.startswith(' ****')
                        if isinstance(Sf, int):
                            txt = str(Sf)
                        elif Sf == 'a':
                            txt = r'\d+'
                        else:
                            raise ValueError('Unsupported final state')
                        fmt = r'^ Excited State\s+' + txt \
                            + r':\s+\S+\s+(?P<val>\d+\.\d+)\b\s+eV.*$'
                        num = -1
                    else:
                        raise NotImplementedError()
            elif qtag == 'dipstr':
                if Si == 0:
                    lnk = 0
                    key = ' Ground to excited state transition electric ' \
                        + 'dipole moments (Au):'
                    if isinstance(Sf, int):
                        sub = Sf + 1
                        def end(s): return True
                    elif Sf == 'a':
                        sub = 1
                        def end(s): return s[1] != ' '
                    else:
                        raise ValueError('Unsupported final state')
                    fmt = r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+' \
                        + r'(?P<val>-?\d+\.\d+)\s+-?\d+\.\d+\s*$'
                    num = -1
                else:
                    raise NotImplementedError()
            elif qtag == 'rotstr':
                if Si == 0:
                    lnk = 0
                    key = ' Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm'\
                        + '/Gauss)'
                    if isinstance(Sf, int):
                        sub = Sf+1
                        def end(s): return True
                    elif Sf == 'a':
                        sub = 1
                        def end(s): return s[1] != ' '
                    else:
                        raise ValueError('Unsupported final state')
                    num = -1
                    fmt = r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+' \
                        + r'(?P<val>-?\d+\.\d+)\s+-?\d+\.\d+\s*$'
                else:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()
            #     if qtag == 1 and dord == 0:
            #         keywords = ['ETran scalars', 'SCF Energy']
        else:
            if rsta == 'c':
                lnk1 = [lnk0]
                key1 = [key0]
                sub1 = [sub0]
                end1 = [end0]
                fmt1 = [fmt0]
                num1 = [num0]
            elif isinstance(rsta, int):
                lnk1 = []
                key1 = []
                sub1 = []
                end1 = []
                fmt1 = []
                num1 = []
            else:
                raise NotImplementedError()
            # Add quantity specific subgroups
            if qtag in ('ramact', 'roaact') and qopt != 'static':
                # Information on setup and incident frequency
                # -- Harmonic level
                if qopt != 'static' and qlvl == 'H':
                    lnk1.append(716)
                    key1.append(' and normal coordinates:')
                    sub1.append(1)
                    end1.append(lambda s: not s.startswith(' Incident light:'))
                    num1.append(0)
                    fmt1.append(
                        r'^\s+Incident light \(\S+\):\s+(?P<val>-?\d.*)\s*$')
                # -- Anharmonic level
                lnk1.append(-717)
                if qopt == 'dynamic':
                    key1.append(
                        '          Anharmonic Raman Spectroscopy (Dynamic)')
                    end1.append(
                        lambda s: s.startswith('     =====')
                        or s.startswith(' GradGrad'))
                    sub1.append(2)
                else:
                    fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                    key1.append(fmt.format(qopt))
                    end1.append(
                        lambda s: s.startswith('     =====')
                        or s.startswith(' GradGrad')
                        or s.startswith(' ## INCIDENT WAVENUMBER:'))
                    sub1.append(0)
                fmt1.append(r'^ ##+ (?:INCIDENT WAVENUMBER|MEASUREMENT '
                            + r'INFORMATION):\s+(?P<val>\S+).*\s+##+\s*$')
                num1.append(0)
            if qtag == 1:
                if dord == 0:
                    if rsta == 'c':
                        lnk1.extend((0, 0))
                        key1.extend((' SCF Done:', ' E2'))
                        sub1.extend((0, 0))
                        end1.extend((lambda s: True, lambda s: True))
                        fmt1.extend(
                            (r'^ SCF Done:\s+' +
                             r'(?P<val>E\(.+?\)\s+=\s+-?\d+\.\d+' +
                             r'(?:D[+-]\d+)?)\s*',
                             r'^\s+E2.*=.*' +
                             r'(?P<val>E\(?.+?\)?\s+=\s+-?\d+\.\d+' +
                             r'(?:D[+-]\d+)?)\s*'))
                        if qopt == 'first':
                            num1.extend((0, 0))
                        elif qopt == 'last':
                            num1.extend((-1, -1))
                        else:
                            num1.extend((1, 1))
                    elif type(rsta) is int:
                        raise NotImplementedError()
                elif dord == 1:
                    raise NotImplementedError()
                elif dord == 2:
                    raise NotImplementedError()
            elif qtag == 101:
                if dord == 0:
                    if rsta == 'c':
                        lnk1.append(601)
                        key1.append(' Dipole moment (field-independent basis,'
                                    + ' Debye):')
                        sub1.append(1)
                        end1.append(lambda s: True)
                        fmt1.append(r'^\s+(?P<val>(?:\s*[XYZ]=\s+-?\d+\.\d+)'
                                    + r'{3}).*$')
                        num1.append(0)
            elif qtag == 'dipstr':
                if qlvl == 'H':
                    lnk1.extend([-716, 716, -717])
                    key1.extend([' and normal coordinates:',
                                 ' and normal coordinates:',
                                 '        Dipole strengths (DS)'])
                    sub1.extend([1, 1, 4])
                    end1.extend([
                        lambda s: s.startswith(' - Thermochemistry'),
                        lambda s: s.startswith(' - Thermochemistry'),
                        lambda s: s.startswith(' -----')])
                    fmt1.extend([r'^\s+Dipole strength  --- \s*'
                                 + r'(?P<val>\d.*)\s*$',
                                 r'^\s+Dip. str.   -- \s*(?P<val>\d.*)\s*$',
                                 r'^\s+\d+\(\d+\)\s+'
                                 + r'(?:-?\d+\.\d+\s+|\*+\s+){2}'
                                 + r'(?P<val>-?\d+\.\d+|\*+)\s+'
                                 + r'(?:-?\d+\.\d+|\*+)\s*$'])
                    num1.extend([0, -1, 0])
                elif qlvl == 'A':
                    lnk1.append(717)
                    key1.append('        Dipole strengths (DS)')
                    sub1.append(3)
                    end1.append(lambda s: s.startswith('     =====')
                                or s.startswith(' GradGrad'))
                    fmt1.append(r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'
                                + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$')
                    num1.append(0)
                else:
                    raise NotImplementedError()
            elif qtag == 'rotstr':
                if qlvl == 'H':
                    lnk1.extend([-716, 716, -717])
                    key1.extend([' and normal coordinates:',
                                 ' and normal coordinates:',
                                 '        Rotational strengths (RS)'])
                    sub1.extend([1, 1, 4])
                    end1.extend([
                        lambda s: s.startswith(' Harmonic frequencies'),
                        lambda s: s.startswith(' - Thermochemistry'),
                        lambda s: s.startswith(' -----')])
                    fmt1.extend([r'^\s+Rot. strength --- \s*'
                                 + r'(?P<val>-?\d.*)\s*$',
                                 r'^\s+Rot. str.   -- \s*(?P<val>-?\d.*)\s*$',
                                 r'^\s+\d+\(\d+\)\s+'
                                 + r'(?:-?\d+\.\d+\s+|\*+\s+){2}'
                                 + r'(?P<val>-?\d+\.\d+|\*+)\s+'
                                 + r'(?:-?\d+\.\d+|\*+)\s*$'])
                    num1.extend([0, -1, 0])
                elif qlvl == 'A':
                    lnk1.append(717)
                    key1.append('        Rotational strengths (RS)')
                    sub1.append(3)
                    end1.append(lambda s: s.startswith('     ====='))
                    fmt1.append(r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'
                                + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$')
                    num1.append(0)
                else:
                    raise NotImplementedError()
            elif qtag == 'ramact':
                if qlvl == 'H':
                    lnk1.extend([-716, 716])
                    key1.extend([' and normal coordinates:',
                                 ' and normal coordinates:'])
                    end1.extend([
                        lambda s: s.startswith(' Harmonic frequencies'),
                        lambda s: s.startswith(' - Thermochemistry')])
                    sub1.extend([1, 1])
                    num1.extend([0, -1])
                    if qopt == 'static':
                        fmt1.extend(
                            [r'^\s+Raman Activities ---\s+(?P<val>-?\d.*)\s*$',
                             r'^\s+Raman Activ -- \s*(?P<val>-?\d.*)\s*$'])
                    else:
                        lnk1.append(716)
                        key1.append(' and normal coordinates:')
                        end1.append(
                            lambda s: s.startswith(' - Thermochemistry'))
                        sub1.append(1)
                        num1.append(-1)
                        fmt1.extend(
                            [r'^\s+Raman Activ Fr=\s?\d --- \s*'
                                 + r'(?P<val>-?\d.*)\s*$',
                             r'^\s+RamAct Fr=\s?\d+--\s+(?P<val>\d.*)\s*$',
                             r'^\s+Raman\d Fr=\s?\d+--\s+(?P<val>\d.*)\s*$'])
                    lnk1.append(-717)
                    if qopt in ('static', 'dynamic'):
                        fmt = 'Anharmonic Raman Spectroscopy ({})'
                        txt = fmt.format(qopt.capitalize())
                        key1.append('{:>49s}'.format(txt))
                        end1.append(
                            lambda s: s.startswith('     =====')
                            or s.startswith(' GradGrad'))
                        sub1.append(2)
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qopt))
                        end1.append(
                            lambda s: s.startswith('     =====')
                            or s.startswith(' GradGrad')
                            or s.startswith(' ## INCIDENT WAVENUMBER:'))
                        sub1.append(1)
                    fmt1.append(r'^\s+\d+\(1\)\s+'
                                + r'(?:-?\d+\.\d+\s+|\*+\s+){2}'
                                + r'(?P<val>-?\d+\.\d+|\*+)\s+'
                                + r'(?:-?\d+\.\d+|\*+)\s*$')
                    num1.append(0)
                elif qlvl == 'A':
                    lnk1.append(717)
                    if qopt in ('static', 'dynamic'):
                        fmt = 'Anharmonic Raman Spectroscopy ({})'
                        txt = fmt.format(qopt.capitalize())
                        key1.append('{:>49s}'.format(txt))
                        end1.append(
                            lambda s: s.startswith('     =====')
                            or s.startswith(' GradGrad'))
                        sub1.append(2)
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qopt))
                        end1.append(
                            lambda s: s.startswith('     =====')
                            or s.startswith(' GradGrad')
                            or s.startswith(' ## INCIDENT WAVENUMBER:'))
                        sub1.append(1)
                    fmt1.append(r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'
                                + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$')
                    num1.append(0)
                else:
                    raise NotImplementedError()
            elif qtag == 'roaact':
                if qlvl == 'H':
                    lnk1.extend([716, -717])
                    key1.append(' and normal coordinates:')
                    end1.append(lambda s: s.startswith(' - Thermochemistry'))
                    if qopt == 'dynamic':
                        key1.append(
                            '             Anharmonic Raman Optical Activity')
                        end1.append(
                            lambda s:
                                s.strip() == 'Dimensionless circular '
                                + 'intensity difference (CID)')
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qopt))
                        end1.append(
                            lambda s:
                                s.strip() == 'Dimensionless circular '
                                + 'intensity difference (CID)'
                                or s.startswith(' ## INCIDENT WAVENUMBER:'))
                    sub1.extend([1, 1])
                    fmt1.extend([r'^\s+ROA\d\s+ Fr= ?\d+-- \s*'
                                 + r'(?P<val>-?\d.*)\s*$',
                                 r'^\s+\d+\(1\)\s+'
                                 + r'(?:-?\d+\.\d+\s+|\*+\s+){2}'
                                 + r'(?P<val>-?\d+\.\d+|\*+)\s+'
                                 + r'(?:-?\d+\.\d+|\*+)\s*$'])
                    num1.extend([-1, 0])
                elif qlvl == 'A':
                    lnk1.append(717)
                    if qopt == 'dynamic':
                        key1.append(
                            '             Anharmonic Raman Optical Activity')
                        end1.append(
                            lambda s:
                                s.strip() == 'Dimensionless circular '
                                + 'intensity difference (CID)')
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qopt))
                        end1.append(
                            lambda s:
                                s.strip() == 'Dimensionless circular '
                                + 'intensity difference (CID)'
                                or s.startswith(' ## INCIDENT WAVENUMBER:'))
                    sub1.append(1)
                    fmt1.append(r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'
                                + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$')
                    num1.append(0)
                else:
                    raise NotImplementedError()
        #     elif qtag == 50:
        #         raise NotImplementedError()
        #     elif qtag == 91:
        #         raise NotImplementedError()
        #     elif qtag == 92:
        #         raise NotImplementedError()
        #     elif qtag == 93:
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
            lnk = tuple(lnk1)
            key = tuple(key1)
            sub = tuple(sub1)
            end = tuple(end1)
            fmt = tuple(fmt1)
            num = tuple(num1)

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
        qtag, qopt, dord, dcrd, rsta, qlvl = qdict[qlabel]
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
        data[qlabel] = {}
        # Basic Properties/Quantities
        # ---------------------------
        if qtag == 'route':
            fmt = '{:d}{:02d}'
            _dat = []
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
                    _dat.append((fulllink, num+1, iops, toline))
            data[qlabel]['data'] = _dat
        elif qtag == 'natoms':
            data[qlabel]['data'] = int(datablocks[iref][0])
        elif qtag == 'atcrd' or qtag == 2:
            _dat = []
            # By default, we choose the standard orientation if present
            if datablocks[last]:
                i = last
            else:
                i = first
            for line in datablocks[i]:
                _dat.append([float(item)*__ang2au for item in line.split()])
            data[qlabel]['data'] = _dat
        elif qtag == 'atmas':
            _dat = []
            for line in datablocks[iref]:
                _dat.extend(float(item) for item in line.split())
            data[qlabel]['data'] = _dat
        elif qtag in ('atnum',):
            data[qlabel]['data'] = []
            for line in datablocks[iref]:
                data[qlabel]['data'].extend(int(item) for item in line.split())
        elif qtag == 'swopt':
            # The last line are the final dashes
            data[qlabel]['data'] = ' '.join(datablocks[iref][:-1])
        elif qtag == 'molsym':
            raise NotImplementedError()
        elif qtag == 'swver':
            txt = r'\s*Gaussian (\w+):\s+(\w+)-(\w{3})Rev([\w.+]+) {1,2}' \
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
                data[qlabel]['data'] = int(datablocks[iref])
            else:
                # For a robust def of nvib, we need the symmetry and
                #   the number of frozen atoms. For now, difficult to do.
                raise NotImplementedError()
        elif qtag in ('hessvec', 'hessval'):
            if iref >= 0:
                data[qlabel]['data'] = datablocks[iref]
        # Vibronic Information
        # --------------------
        elif qtag == 'fcdat':
            if qopt == 'SimInf':
                def getval(s):
                    return s.split(':', maxsplit=1)[-1].strip().title()
                for line in datablocks[iref]:
                    if 'Temperature' in line:
                        if 'not' in line:
                            val = None
                        else:
                            val = float(line.split()[-1][:-1])
                        data[qlabel]['temp'] = val
                    elif 'framework' in line:
                        data[qlabel]['frame'] = getval(line)
                    elif 'Spectroscopy' in line:
                        data[qlabel]['spec'] = getval(line)
                    elif 'Model' in line:
                        data[qlabel]['model'] = getval(line)
                    elif 'electronic transition moment' in line:
                        data[qlabel]['tmom'] = \
                            line.split(':', maxsplit=1)[-1].strip()
            elif qopt == 'JMat':
                if datablocks[last]:
                    i = last
                else:
                    i = first
                if 'diagonal' in datablocks[i]:
                    data[qlabel]['data'] = []
                else:
                    data[qlabel]['data'] = []
                    N = 0
                    for line in datablocks[i]:
                        cols = line.split()
                        irow = int(cols[0]) - 1
                        if irow == N:
                            data[qlabel]['data'].append([])
                            N += 1
                        data[qlabel]['data'][irow].extend(
                            [float(item.replace('D', 'e'))
                             for item in cols[1:]])
            elif qopt == 'JMatF':
                if 'diagonal' in datablocks[iref]:
                    data[qlabel]['data'] = []
                else:
                    data[qlabel]['data'] = []
                    N = 0
                    for line in datablocks[iref]:
                        cols = line.split()
                        irow = int(cols[0]) - 1
                        if irow == N:
                            data[qlabel]['data'].append([])
                            N += 1
                        data[qlabel]['data'][irow].extend(
                            [float(item.replace('D', 'e'))
                             for item in cols[1:]])
            elif qopt == 'KVec':
                data[qlabel]['data'] = []
                for line in datablocks[iref]:
                    data[qlabel]['data'].append(float(line.replace('D', 'e')))
            elif qopt == 'SRAMat':
                data[qlabel]['data'] = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data[qlabel]['data'].append([])
                        N += 1
                    data[qlabel]['data'][irow].extend(
                        [float(item.replace('D', 'e'))
                            for item in cols[1:]])
            elif qopt == 'SRBVec':
                data[qlabel]['data'] = []
                for line in datablocks[iref]:
                    data[qlabel]['data'].append(float(line.replace('D', 'e')))
            elif qopt == 'SRCMat':
                data[qlabel]['data'] = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data[qlabel]['data'].append([])
                        N += 1
                    data[qlabel]['data'][irow].extend(
                        [float(item.replace('D', 'e'))
                            for item in cols[1:]])
            elif qopt == 'SRDVec':
                data[qlabel]['data'] = []
                for line in datablocks[iref]:
                    data[qlabel]['data'].append(float(line.replace('D', 'e')))
            elif qopt == 'SREMat':
                data[qlabel]['data'] = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data[qlabel]['data'].append([])
                        N += 1
                    data[qlabel]['data'][irow].extend(
                        [float(item.replace('D', 'e'))
                            for item in cols[1:]])
            elif qopt == 'Spec':
                # Look for last blocks first, which should contain all blocks:
                discard = []
                if last != first:
                    nblocks = 0
                    nyaxes = 0
                    for i, bloc in enumerate(datablocks[last]):
                        lbloc = len(bloc)
                        if lbloc > 1:
                            # 3 for: Legend title, X axis, Intensity unit
                            nblocks += 1
                            nyaxes += len(bloc[0].split()) - 1
                        else:
                            # Incorrect bloc
                            discard.append(i)
                else:
                    nblocks = 1
                    nyaxes = 0
                if nblocks == 0:
                    msg = 'Inconsistency in spectral data. ' +\
                        + 'This should not happen!'
                    raise IndexError(msg)
                if discard:
                    discard.reverse()
                    for i in discard:
                        del(datablocks[last][i])
                # First block contains the full initial legend block
                # Necessary for the different parameters
                iref = first
                data[qlabel] = {'x': []}
                if nyaxes == 0:
                    nyaxes = len(datablocks[iref][0].split()) - 1
                if nyaxes == 1:
                    yfmt = 'y'
                    data[qlabel]['y'] = []
                else:
                    yfmt = 'y{{idy:0{}d}}'.format(len(str(nyaxes)))
                    for i in range(nyaxes):
                        data[qlabel][yfmt.format(idy=i+1)] = []
                # If multiple blocks, the reading with first parsing method
                #   is wrong since it combines all blocks together
                # We fix it by only considering the last one which should be
                #   always right.
                iref = last
                for bloc in range(nblocks):
                    yax = {}
                    if nblocks > 1:
                        data[qlabel]['x'].append([])
                        xax = data[qlabel]['x'][-1]
                        for i in range(nyaxes):
                            y = yfmt.format(idy=i+1)
                            data[qlabel][y].append([])
                            yax[y] = data[qlabel][y][-1]
                    else:
                        xax = data[qlabel]['x']
                        for i in range(nyaxes):
                            y = yfmt.format(idy=i+1)
                            yax[y] = data[qlabel][y]
                    for line in datablocks[iref][bloc]:
                        cols = [float(item.replace('D', 'e'))
                                for item in line.split()]
                        xax.append(cols[0])
                        for i, item in enumerate(cols[1:]):
                            yax[yfmt.format(idy=i+1)].append(item)
            elif qopt == 'SpcPar':
                # Look for last blocks first, which should contain all blocks:
                discard = []
                if last != first:
                    nblocks = 0
                    nyaxes = 0
                    for i, bloc in enumerate(datablocks[last]):
                        lbloc = len(bloc)
                        if lbloc > 3:
                            # 3 for: Legend title, X axis, Intensity unit
                            nblocks += 1
                            nyaxes += lbloc - 3
                        else:
                            # Incorrect bloc
                            discard.append(i)
                else:
                    nblocks = 1
                    nyaxes = 0
                if nblocks == 0:
                    nblocks = 1
                if discard:
                    discard.reverse()
                    for i in discard:
                        del(datablocks[last][i])
                # First block contains the full initial legend block
                # Necessary for the different parameters
                iref = first
                i = 0
                while datablocks[iref][i].strip() != 'Legend:':
                    i += 1
                    if i >= len(datablocks[iref]):
                        raise IndexError('Legend not found')
                if 'No' in datablocks[iref][i-1]:
                    data[qlabel]['func'] = 'stick'
                    data[qlabel]['hwhm'] = None
                elif 'broadening' in datablocks[iref][i-2]:
                    func = datablocks[iref][i-2].split()[-3].lower()
                    hwhm = float(datablocks[iref][i-1].split()[-2])
                    data[qlabel]['func'] = func
                    data[qlabel]['hwhm'] = hwhm
                else:
                    raise IndexError('Unrecognized broadening function.')
                offset = i + 1
                # offset: num. lines for legend block + legend section title
                if nyaxes == 0:
                    nyaxes = len(datablocks[iref]) - offset - 2
                if nyaxes <= 1:
                    # X, Y, Int, don't use numbering
                    yfmt = 'y'
                else:
                    yfmt = 'y{{idy:0{}d}}'.format(len(str(nyaxes)))
                for line in datablocks[iref][offset:]:
                    res = line.split(':', 1)
                    if len(res) == 2:
                        key, title = res
                    elif 'intensity' in res[0]:
                        title = res[0]
                        key = 'Intensity'
                    else:
                        raise NotImplementedError('Unrecognized legend')
                    if key.strip() == '1st col.':
                        _key = 'x'
                        txt = title.rsplit('in ', maxsplit=1)[1].\
                            replace(')', '').strip()
                        data[qlabel]['unitx'] = txt
                    elif key.strip() == 'Intensity':
                        _key = 'I'
                        txt = title.rsplit('in ', maxsplit=1)[1].\
                            replace(')', '').strip()
                        if data[qlabel]['func'] == 'stick':
                            _desc = 'II:'
                            # Fix a stupidity in the unit in some versions
                            #   for the stick spectrum
                            if txt == 'dm^3.mol^-1.cm^-1':
                                txt = 'dm^3.mol^-1.cm^-2'
                        else:
                            _desc = 'I:'
                        data[qlabel]['unity'] = _desc + txt
                    else:
                        try:
                            idy = int(key.strip()[0]) - 1
                        except ValueError:
                            raise IndexError('Unrecognized key in spc leg.: ' +
                                             key)
                        _key = yfmt.format(idy=idy)
                    if nblocks == 1:
                        data[qlabel][_key] = title.strip()
                    else:
                        if _key[0] == 'I':
                            data[qlabel][_key] = title.strip()
                        else:
                            data[qlabel][_key] = [title.strip()]
                # More than 1 block, read the new ones
                for bloc in range(1, nblocks):
                    for line in datablocks[last][bloc][1:]:
                        res = line.split(':', 1)
                        if len(res) == 2:
                            key, title = res
                        elif 'intensity' in res[0]:
                            title = res[0]
                            key = 'Intensity'
                        else:
                            raise NotImplementedError('Unrecognized legend')
                        if key.strip() == '1st col.':
                            _key = 'x'
                        elif 'col.' in key:
                            try:
                                idy = int(key.strip()[0]) - 1
                            except ValueError:
                                msg = 'Unrecognized key in spc leg.: ' + key
                                raise IndexError(msg)
                            _key = yfmt.format(idy=idy)
                        data[qlabel][_key].append(title.strip())
            elif qopt == 'Conv':
                raise NotImplementedError()
            elif qopt == 'Assign':
                data[qlabel]['T'] = []
                data[qlabel]['E'] = []
                data[qlabel]['I'] = []
                if 'DipStr':
                    qty = 'DS'
                elif 'RotStr':
                    qty = 'RS'
                else:
                    msg = 'Unrecognized spectroscopy-specific quantity'
                    raise IndexError(msg)
                data[qlabel]['other'] = qty
                data[qlabel][qty] = []
                for l1, l2 in zip(datablocks[first], datablocks[last]):
                    txt_E, txt_T = l1.split(':')
                    data[qlabel]['E'].append(float(txt_E.split()[0]))
                    trans = []
                    for state in txt_T.split('->'):
                        trans.append([])
                        for item in state.strip(' |>').split(';'):
                            if '^' in item:
                                i, n = [int(val) for val in item.split('^')]
                            else:
                                i = int(item)
                                n = 0
                            trans[-1].append((i, n))
                    data[qlabel]['T'].append((tuple(trans[0]),
                                              tuple(trans[1])))
                    line = l2.strip(')').split()
                    data[qlabel]['I'].append(float(line[0]))
                    data[qlabel][qty].append(float(line[-1]))
            elif qopt == 'GeomIS':
                data[qlabel]['data'] = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data[qlabel]['data'].append([float(item)*__ang2au
                                                 for item in line.split()])
            elif qopt == 'GeomFS':
                data[qlabel]['data'] = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data[qlabel]['data'].append([float(item)*__ang2au
                                                 for item in line.split()])
            elif qopt == 'ExGeom':
                data[qlabel]['data'] = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data[qlabel]['data'].append([float(item)*__ang2au
                                                 for item in line.split()])
        # Vibrational transitions
        # -----------------------
        elif qtag == 'vlevel':
            if qlvl == 'H':
                for i in range(last, first-1, -1):
                    if datablocks[i]:
                        iref = i
                        break
                else:
                    raise ParseKeyError('Missing quantity in file')
                data[qlabel]['unit'] = 'cm-1'
                i = 0
                for line in datablocks[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[qlabel][i] = float(col)
                        except ValueError:
                            data[qlabel][i] = float('inf')
            elif qlvl == 'A':
                if datablocks[last]:
                    iref = last
                else:
                    iref = first
                data[qlabel]['unit'] = 'cm-1'
                i = 0
                for line in datablocks[iref]:
                    i += 1
                    try:
                        data[qlabel][i] = float(line)
                    except ValueError:
                        data[qlabel][i] = float('inf')
            else:
                raise NotImplementedError()
        elif qtag == 'vtrans':
            if qlvl == 'H':
                for i in range(last, first-1, -1):
                    if datablocks[i]:
                        iref = i
                        break
                else:
                    raise ParseKeyError('Missing quantity in file')
                i = 0
                for line in datablocks[iref]:
                    for col in line.strip().split():
                        i += 1
                        res = col.split('(')
                        data[qlabel][i] = [((0, 0), ), ((int(res[0]), 1), )]
            elif qlvl == 'A':
                i = 0
                for line in datablocks[iref]:
                    i += 1
                    val = []
                    for col in line.strip().split():
                        res = col.split('(')
                        if len(res) == 1:
                            val.append((int(res[0]), 0))
                        else:
                            val.append((int(res[0]),
                                        int(res[1].replace(')', ''))))
                        data[qlabel][i] = [((0, 0), ), tuple(val)]
            else:
                raise NotImplementedError()
        # Anharmonic Information
        # ----------------------
        elif qtag == 'vptdat':
            if qopt == 'CICoef':
                # First, check if reduced-dimensionality used.
                if datablocks[last]:
                    eqv = {}
                    for line in datablocks[last]:
                        j, i = [int(item) for item in line.split('|')]
                        eqv[i] = j

                    def idx(i):
                        return eqv.get(i, 0)
                else:

                    def idx(i):
                        return i
                # Now parses the main assignment:
                data[qlabel]['data'] = {}
                i = 0
                for line in datablocks[first]:
                    if ':' in line:
                        key, dat = line.split(':')
                        i = idx(int(key))
                        if i > 0:
                            data[qlabel][i] = []
                        coef, state = dat.split('x')
                    else:
                        coef, state = line.split('x')
                    if i > 0:
                        modes = state.split(';')
                        dat = []
                        for mode in modes:
                            res = re.split('[|>()]+', mode.strip())
                            if res[0]:
                                i0 = 0
                                i1 = 1
                            else:
                                i0 = 1
                                i1 = 2
                            dat.append((int(res[i0]), int(res[i1])))
                        data[qlabel][i].append((float(coef), tuple(dat)))
            else:
                raise NotImplementedError()
        # State(s)-dependent quantities
        # -----------------------------
        else:
            # Transition moments
            # ^^^^^^^^^^^^^^^^^^
            if type(rsta) is tuple:
                Si, Sf = rsta
                if qtag == 1:
                    if dord == 0:
                        if Si == 0:
                            if isinstance(Sf, int):
                                data[qlabel]['data'] = float(datablocks[iref])
                                data[qlabel]['unit'] = 'eV'
                            else:
                                data[qlabel]['unit'] = 'eV'
                                data[qlabel]['data'] = \
                                    [float(item) for item in datablocks[iref]]
                        else:
                            pass
                    else:
                        pass
                elif qtag == 'dipstr':
                    if Si == 0:
                        if isinstance(Sf, int):
                            data[qlabel]['data'] = float(datablocks[iref])
                            data[qlabel]['unit'] = 'DS:au'
                        else:
                            data[qlabel]['unit'] = 'DS:au'
                            data[qlabel]['data'] = \
                                [float(item) for item in datablocks[iref]]
                    else:
                        pass
                elif qtag == 'rotstr':
                    if Si == 0:
                        if isinstance(Sf, int):
                            data[qlabel]['data'] \
                                = float(datablocks[iref])*1.0e-40
                            data[qlabel]['unit'] = 'RS:esu^2.cm^2'
                        else:
                            data[qlabel]['unit'] = 'RS:esu^2.cm^2'
                            data[qlabel]['data'] = \
                                [float(item) * 1.0e-40
                                 for item in datablocks[iref]]
                    else:
                        pass
                else:
                    pass
            # States-specific Quantities
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^
            else:
                if rsta != 'c':
                    # we should check if the state is the right one.
                    raise NotImplementedError()
                if qtag == 1:
                    if dord == 0:
                        if rsta == 'c':
                            def conv(s):
                                return float(s.split('=')[1].replace('D', 'e'))

                            # 1st element is actually transition information
                            # Ignored as we look for the pure energies
                            fmt = re.compile(r'E\(?(.*?)\)?\s+')
                            iref = first + 1
                            N = 2 if datablocks[last] else 1
                            nblocks = len(datablocks[iref])
                            for i in range(iref, iref+N):
                                if nblocks > 1:
                                    txt = datablocks[i][0][0]
                                    val = [conv(item[0])
                                           for item in datablocks[i]]
                                else:
                                    txt = datablocks[i][0]
                                    val = conv(txt)
                                res = fmt.search(txt)
                                try:
                                    tag = res.group(1)
                                except AttributeError:
                                    msg = 'Unsupported energy format'
                                    raise ParseKeyError(msg) from None
                                data[qlabel][tag] = val
                            data[qlabel]['data'] = data[qlabel][tag]
                            data[qlabel]['unit'] = 'Eh'
                elif qtag == 101:
                    if dord == 0:
                        if rsta == 'c':
                            val = [float(item)/phys_fact('au2Deb') for item in
                                   datablocks[last][0].split()[1::2]]
                            data[qlabel]['data'] = val
                            data[qlabel]['unit'] = 'e.a0'
                        else:
                            raise NotImplementedError()
                    else:
                        raise NotImplementedError()
                elif qtag == 'dipstr':
                    if qlvl == 'H':
                        for i in range(last, first-1, -1):
                            if datablocks[i]:
                                iref = i
                                break
                        else:
                            raise ParseKeyError('Missing quantity in file')
                        data[qlabel]['unit'] = 'DS:esu^2.cm^2'
                        i = 0
                        for line in datablocks[iref]:
                            for col in line.strip().split():
                                i += 1
                                try:
                                    data[qlabel][i] = float(col)*1.0e-40
                                except ValueError:
                                    data[qlabel][i] = float('inf')
                    elif qlvl == 'A':
                        data[qlabel]['unit'] = 'DS:esu^2.cm^2'
                        i = 0
                        for line in datablocks[iref]:
                            i += 1
                            try:
                                data[qlabel][i] = float(line)*1.0e-40
                            except ValueError:
                                data[qlabel][i] = float('inf')
                    else:
                        raise NotImplementedError()
                elif qtag == 'rotstr':
                    if qlvl == 'H':
                        for i in range(last, first-1, -1):
                            if datablocks[i]:
                                iref = i
                                break
                        else:
                            raise ParseKeyError('Missing quantity in file')
                        data[qlabel]['unit'] = 'RS:esu^2.cm^2'
                        i = 0
                        for line in datablocks[iref]:
                            for col in line.strip().split():
                                i += 1
                                try:
                                    data[qlabel][i] = float(col)*1.0e-44
                                except ValueError:
                                    data[qlabel][i] = float('inf')
                    elif qlvl == 'A':
                        data[qlabel]['unit'] = 'RS:esu^2.cm^2'
                        i = 0
                        for line in datablocks[iref]:
                            i += 1
                            try:
                                data[qlabel][i] = float(line)*1.0e-44
                            except ValueError:
                                data[qlabel][i] = float('inf')
                    else:
                        raise NotImplementedError()
                elif qtag == 'ramact':
                    data[qlabel] = _parse_logdat_ramact(qopt, qlvl, datablocks,
                                                        first, last)
                elif qtag == 'roaact':
                    data[qlabel] = _parse_logdat_ramact(qopt, qlvl, datablocks,
                                                        first, last, ROA=True)
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
             error_noqty: bool = True) -> TypeQData:
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


def _parse_logdat_ramact(qopt: str,
                         qlvl: str,
                         datablocks: tp.Sequence[tp.Sequence[int]],
                         first: int,
                         last: int,
                         ROA: bool = False) -> tp.Dict[str, tp.Any]:
    """Sub-function to parse GLog data for Raman activity

    Sub-functions dedicated to parsing data stored in datablocks.

    Parameters
    ----------
    """
    data = {}
    for i in range(last, first-1, -1):
        if datablocks[i]:
            iref = i
            break
    else:
        raise ParseKeyError('Missing quantity in file')
    if qlvl == 'H':
        if iref == last:
            data['unit'] = 'ROA:Ang^6' if ROA else 'RA:Ang^6'
            an_blk = True
        else:
            data['unit'] = 'ROA:amu.Ang^4' if ROA else 'RA:amu.Ang^4'
            an_blk = False
        # We create the database based on the quantity
        if qopt == 'static':
            i = 0
            for line in datablocks[iref]:
                for col in line.strip().split():
                    i += 1
                    try:
                        data[i] = float(col)
                    except ValueError:
                        data[i] = float('inf')
        elif an_blk:
            incfrq = None
            setups = []
            for item in datablocks[first+2]:
                try:
                    _ = float(item)
                    incfrq = item
                    data[incfrq] = {}
                except ValueError:
                    if incfrq is None:
                        msg = 'Wrong block structure for Raman/ROA ' \
                            + 'spectroscopy (incident freq/setup)'
                        raise ParseKeyError(msg)
                    data[incfrq][item] = {}
                    setups.append((incfrq, item))
            nlines = len(datablocks[iref])/len(setups)
            block = 0
            iline = 0
            for line in datablocks[iref]:
                iline += 1
                if iline % nlines == 1:
                    d = data[setups[block][0]][setups[block][1]]
                    block += 1
                    i = 0
                for col in line.strip().split():
                    i += 1
                    try:
                        d[i] = float(col)
                    except ValueError:
                        d[i] = float('inf')
        else:
            incfreqs = [item for line in datablocks[first+1]
                        for item in line.split()]
            if qopt == 'dynamic':
                for freq in incfreqs:
                    data[freq] = {
                        'SCP(180)u': {}, 'SCP(90)z': {}, 'DCPI(180)': {}}
                iline = 0
                ifreq = 0
                ioff = 0
                i = 0
                for line in datablocks[iref]:
                    iline += 1
                    block = iline % 3
                    if block == 1:
                        dfreq = data[incfreqs[ifreq % len(incfreqs)]]
                        ifreq += 1
                        ioff += i+1
                    d = dfreq[('SCP(180)u', 'SCP(90)z', 'DCPI(180)')[block-1]]
                    for i, col in enumerate(line.strip().split()):
                        try:
                            d[ioff+i] = float(col)
                        except ValueError:
                            d[ioff+i] = float('inf')
            else:
                if qopt in incfreqs:
                    ref_freq = incfreqs.index(qopt)
                    data[incfreqs[ref_freq]] = {
                        'SCP(180)u': {}, 'SCP(90)z': {}, 'DCPI(180)': {}}
                    iline = 0
                    ifreq = 0
                    iref = 0
                    i = 0
                    for line in datablocks[iref]:
                        iline += 1
                        block = iline % 3
                        if block == 1:
                            if ifreq % len(incfreqs) != ref_freq:
                                ifreq += 1
                                continue
                            dfreq = data[incfreqs[ref_freq]]
                            ifreq += 1
                            iref += i+1
                        d = dfreq[('SCP(180)u', 'SCP(90)z',
                                   'DCPI(180)')[block-1]]
                        for i, col in enumerate(line.strip().split()):
                            try:
                                d[iref+i+1] = float(col)
                            except ValueError:
                                d[iref+i+1] = float('inf')
                else:
                    raise ParseKeyError('Missing incident frequency')
    elif qlvl == 'A':
        incfrq = None
        setups = []
        for item in datablocks[first+1]:
            try:
                _ = float(item)
                incfrq = item
                data[incfrq] = {}
            except ValueError:
                if incfrq is None:
                    msg = 'Wrong block structure for Raman/ROA spectroscopy ' \
                        + '(incident freq/setup)'
                    raise ParseKeyError(msg)
                data[incfrq][item] = {}
                setups.append((incfrq, item))
        nlines = len(datablocks[iref])/len(setups)
        block = 0
        iline = 0
        for line in datablocks[iref]:
            iline += 1
            if iline % nlines == 1:
                d = data[setups[block][0]][setups[block][1]]
                block += 1
                i = 0
            for col in line.strip().split():
                i += 1
                try:
                    d[i] = float(col)
                except ValueError:
                    d[i] = float('inf')
    else:
        raise NotImplementedError()
    return data
