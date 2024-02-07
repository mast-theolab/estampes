"""Low-level operations on Gaussian output (log) files.

Provides low-level interfaces to manipulate/extract data in Gaussian
output files.

"""

import os  # Used for file existence check
import re  # Used to find keys in log file
import typing as tp
from math import sqrt

from estampes.base import ParseKeyError, QuantityError, TypeDGLog, TypeQData, \
    TypeQInfo, QLabel, QData
from estampes.data.physics import PHYSFACT, phys_fact
from estampes.parser.functions import parse_qlabels


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

xyz2id = {'X': 0, 'Y': 1, 'Z': 2}


# ==============
# Module Classes
# ==============

class GLogIO(object):
    """Carry file operations on Gaussian log file.

    Main class to handle Gaussian output file operations.
    """

    def __init__(self, fname: str,
                 load_pos: bool = True) -> None:
        """Initialize Gaussian logfile IO instance.

        Initializes an instance of GLogIO.

        Parameters
        ----------
        fname
            Log file name.
        load_pos
            Preloads links blocks positions if easily recognizable (#P).
        """
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
        """Filename associated to the GLog object."""
        return self.__fname

    @filename.setter
    def filename(self, name: str) -> None:
        if not os.path.exists(name):
            raise FileNotFoundError('Gaussian output file not found')
        self.__fname = name

    @property
    def version(self) -> tp.Dict[str, str]:
        """Version of Gaussian used to generate the log file.

        Returns the version as a dictionary with two keys:

        major
            Major revision.
        minor
            Minor revision.
        """
        return {key: self.__gversion[key] for key in ('major', 'minor')}

    @property
    def full_version(self) -> tp.Tuple[str, tp.Any]:
        """Full version of Gaussian, for the parser interface.

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

    def get_head(self):
        """Return the header information: Version, Route."""
        keydata = []
        qtydata = {}
        key2blk = {}
        i = 0
        for item in ('route', 'swopt', 'swver'):
            qlab = QLabel(quantity=item)
            key2blk[item] = (i, i)
            i += 1
            link, key, skips, fmt, end, num = qlab_to_linkdata(qlab)
            keydata.append((key, link, skips, re.compile(fmt), num, end))
            qtydata[item] = qlab
        ndata, data = self.read_data(*keydata)
        dobjs = parse_data(qtydata, key2blk, ndata, data)
        self.__route = dobjs['route'].data
        self.__links = sorted(set([int(item[0]) for item in self.__route]))
        self.__gversion = dobjs['swver'].data
        self.__rte_opt = dobjs['swopt'].data
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
                  *to_find: TypeKData) -> TypeDGLog:
        """Extract data corresponding to the keys to find.

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
            """Delete a block in the lookup tables.

            Returns
            -------
            int
                Status, as integer
                0: block removed
                1: keylist item removed
            """
            istat = 0
            if nocc[iblock] == 0:
                i, j = block2id[iblock]
                del keydata[i][j]
                if not keydata[i]:
                    del keydata[i]
                    del keylist[i]
                    istat = 1
                for k, blk2id in enumerate(block2id):
                    a, b = blk2id
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
            with open(self.filename, 'r', encoding="utf-8") as fobj:
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
        """Store link header positions in the file if available.

        Loads the keys present in the file and pointers to their
        position to speed up their search.
        Data type and block information are also stored.
        """
        # link_heads = {
        #     # 1: 'Entering Gaussian System,',
        #     1: 'Entering Link 1,',
        #     601: 'Population analysis using the SCF Density.',
        #     716: 'Full mass-weighted force constant matrix:',
        #     717: 'Second-order Perturbative Anharmonic Analysis',
        #     718: 'Generation of the Franck-Condon spectrum'
        # }
        # to_search = re.compile(r'''\
        #     (?P<title>[\w\s]+?)\s*  # Key
        #     \b(?P<type>[IRC])\b\s*  # Data type
        #     (?P<block>N=)?\s+  # N= only set for non-scalar data
        #     (?P<value>[\d\-\+\.E]+)  # Block size (N=) or scalar value
        #     \n''', re.VERBOSE)

        # keys = {}
        # with open(self.filename, 'r', encoding="utf-8") as fobj:
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


def qlab_to_linkdata(qlab: QLabel, gver: tp.Optional[str] = None) -> TypeQKwrd:
    """Return relevant keyword(s) for a given quantity.

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
    qlab
        Quantity label
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
    ParseKeyError
        The quantity is known to not be available.

    Notes
    -----
    - `n` refers to all available states.
    """
    lnk = None
    key = None
    sub = None
    fmt = None
    num = None
    if qlab.label == 'route':  # Log specific extraction, not provided in base
        lnk = 1
        key = ' Cite this work as:'
        sub = ' -'
        def end(s): return s.startswith(' Charge =')
        fmt = r'^ (?P<val>\d+\/(?:,?\d+=-?\d+)*\/(?:,?\d+)+(?:\(-\d\))?);\s*$'
        num = 0
    elif qlab.label == 'natoms':
        lnk = 101
        key = ' NAtoms='
        sub = 0
        def end(_s): return True
        fmt = r'^ NAtoms=\s+(?P<val>\d+)\s+N\w+=\s+.*$'
        num = 0
    elif qlab.label == 'nvib':
        # We load anyway the frequencies to count them.
        # This is a bit of an absurd way to proceed but only way to be
        #   sure that we have the correct number since Gaussian does not
        #   list it explicitly in the output
        lnk = 716
        key = ' and normal coordinates:'
        sub = 1
        def end(s): return s.startswith(' - Thermochemistry')
        fmt = r'^\s+Frequencies -- \s*(?P<val>\d.*)\s*$'
        num = -1
    elif qlab.label == 'atmas':
        lnk = (-101, 716)
        key = (' NAtoms= ', ' - Thermochemistry -')
        sub = (0, 3)
        end = (lambda s: s.startswith(' Leave Link'),
               lambda s: s.startswith(' Molecular Mass:'))
        fmt = (r'^ AtmWgt=\s+(?P<val>(\s+\d+\.\d+)+)\s*$',
               r'^ Atom\s+\d+ has atomic number\s+\d+ and '
               + r'mass\s+(?P<val>\d+\.\d+)\s*$')
        num = (0, 0)
    elif qlab.label == 'atnum':
        lnk = (0, 0)
        key = ('                         Standard orientation:',
               '                          Input orientation:')
        sub = (5, 5)
        end = (lambda s: s.startswith(' ------'),
               lambda s: s.startswith(' ------'))
        txt = r'^\s+\d+\s+(?P<val>\d+)\s+\d+(?:\s+-?\d+\.\d+){3}\s*$'
        fmt = (txt, txt)
        num = (0, 0)
    elif qlab.label == 'molsym':
        lnk = 101
        key = ' Framework group '
        sub = 0
        def end(_s): return True
        fmt = r'^ Framework group \s+(?P<val>[^ ]+)\s*$'
        num = 0
    elif qlab.label == 'atcrd' or qlab.label == 2:
        lnk = (0, 0)
        key = ('                         Standard orientation:',
               '                          Input orientation:')
        sub = (5, 5)
        end = (lambda s: s.startswith(' ------'),
               lambda s: s.startswith(' ------'))
        txt = r'^(?:\s+\d+){3}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
        fmt = (txt, txt)
        if qlab.kind == 'first':
            num = (0, 0)
        elif qlab.kind == 'last':
            num = (-1, -1)
        else:
            num = (1, 1)
    elif qlab.label == 'hessvec':
        lnk = (-716, 716)
        key = (' and normal coordinates:',
               ' and normal coordinates:')
        sub = (1, 1)
        end = (lambda s: s.startswith(' - Thermochemistry'),
               lambda s: s.startswith(' - Thermochemistry'))
        fmt = (r'^\s+(?P<val>(?:\d+\s+){3}(?:\s+-?\d\.\d+){1,5})\s*$',
               r'^\s+(?P<val>(?:\d+\s+){2}(?:\s+-?\d\.\d+){1,9})\s*$')
        num = (0, -1)
    elif qlab.label == 'hessdat':
        lnk = (-716, 716)
        key = (' and normal coordinates:',
               ' and normal coordinates:')
        sub = (1, 1)
        end = (lambda s: s.startswith(' - Thermochemistry'),
               lambda s: s.startswith(' - Thermochemistry'))
        if qlab.kind == 'freq':
            fmt = (r'^\s+Frequencies --- \s*(?P<val>\d.*)\s*$',
                   r'^\s+Frequencies -- \s*(?P<val>\d.*)\s*$')
        elif qlab.kind == 'redmas':
            fmt = (r'^\s+Reduced masses --- \s*(?P<val>\d.*)\s*$',
                   r'^\s+Red. masses -- \s*(?P<val>\d.*)\s*$')
        else:
            raise NotImplementedError('Unknown subopt for HessDat')
        num = (0, -1)
    elif qlab.label == 'swopt':
        lnk = 1
        key = ' Cite this work as:'
        sub = ' #'
        def end(s): return s.startswith(' -')
        fmt = r'^ (?P<val>.*)$'
        num = 0
    elif qlab.label == 'swver':
        lnk = 1
        key = ' ****'
        sub = 1
        def end(_s): return True
        fmt = r'^ (?P<val>Gaussian (?:\d\d|DV):\s.*)$'
        num = 0
    elif qlab.label == 'intens':
        if qlab.kind == 'IR':
            if qlab.level == 'H':
                lnk = (-716, 716, -717)
                key = (' and normal coordinates:',
                       ' and normal coordinates:',
                       '        Integrated intensity (I)')
                sub = (1, 1, 4)
                end = (lambda s: s.startswith(' - Thermochemistry'),
                       lambda s: s.startswith(' - Thermochemistry'),
                       lambda s: s.startswith(' -----'))
                fmt = (r'^\s+IR Intensities --- \s*(?P<val>\d.*)\s*$',
                       r'^\s+IR Inten    -- \s*(?P<val>\d.*)\s*$',
                       r'^\s+\d+\(\d+\)\s+(?:-?\d+\.\d+\s+|\*+\s+){2}'
                       + r'(?P<val>-?\d+\.\d+|\*+)\s+'
                       + r'(?:-?\d+\.\d+|\*+)\s*$')
                num = (0, -1, 0)
            elif qlab.level == 'A':
                lnk = 717
                key = '        Integrated intensity (I)'
                sub = 3
                def end(s): return s.startswith(' Units:')
                fmt = r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+' \
                      + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$'
                num = 0
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()
    elif qlab.label == 'fcdat':
        if qlab.kind == 'SimInf':
            lnk = -718
            key = '               Information on the Simulation'
            sub = 2
            def end(s): return s.startswith('     ====')
            fmt = r'^\s+(?P<val>.*\w.*)\s*$'
            num = 0
        elif qlab.kind == 'JMat':
            lnk = (-718, -718)
            key = (' Duschinsky matrix', ' Final Duschinsky matrix')
            sub = (2, 2)
            end = (lambda s: not s.strip(),
                   lambda s: not s.strip())
            txt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            fmt = (txt, txt)
            num = (-1, -1)
        elif qlab.kind == 'JMatF':
            lnk = -718
            key = ' Full Duschinsky matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qlab.kind == 'KVec':
            lnk = -718
            key = ' Shift Vector'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$'
            num = 0
        elif qlab.kind == 'SRAMat':
            lnk = -718
            key = ' A Matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qlab.kind == 'SRBVec':
            lnk = -718
            key = ' B Vector'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$'
            num = 0
        elif qlab.kind == 'SRCMat':
            lnk = -718
            key = ' C Matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qlab.kind == 'SRDVec':
            lnk = -718
            key = ' D Vector'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$'
            num = 0
        elif qlab.kind == 'SREMat':
            lnk = -718
            key = ' E Matrix'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            num = 0
        elif qlab.kind == 'Spec':
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
        elif qlab.kind == 'SpcPar':
            lnk = (718, 718)
            key = ('                       Final Spectrum',
                   ' Legend:')
            sub = (3, 0)
            end = (lambda s: s.startswith(' -----------'),
                   lambda s: s.startswith(' -----------'))
            fmt = (r'^\s+(?P<val>.*\w.*)\s*$',  # \w to exclude empty lines
                   r'^\s+(?P<val>.*\w.*)\s*$')
            num = (0, 1)
        elif qlab.kind == 'Conv':
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
        elif qlab.kind == 'Assign':
            lnk = (-718, -718)
            key = ('                 Information on Transitions',
                   '                 Information on Transitions')
            sub = (2, 2)
            end = (lambda s: s.startswith('     ===='),
                   lambda s: s.startswith('     ===='))
            fmt = (r'^\s+Energy =\s+(?P<val>-?\d+\.\d+ cm.-1: .*)\s*$',
                   r'^\s+-. Intensity =\s+(?P<val>.*)\s*$')
            num = (0, 0)
        elif qlab.kind == 'RedDim':
            lnk = -718
            key = ' Reduced system'
            sub = 2
            def end(s): return not s.strip()
            fmt = r'^\s+(?P<val>\d+\s+=\s+\d+\s+\d+\s+=\s+\d+)\s*$'
            num = 0
        elif qlab.kind == 'E(0-0)':
            lnk = 718
            key = '                 Information on Transitions'
            sub = 2
            def end(s): return s.startswith('     ====')
            fmt = r'^\s+Energy of the 0-0 transition:\s+' \
                + r'(?P<val>-?\d+\.\d+ cm\S+)\s*$'
            num = 0
        elif qlab.kind == 'GeomIS':
            lnk = 718
            key = '              New orientation in initial state'
            sub = 5
            def end(s): return s.startswith(' ------')
            fmt = r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
        elif qlab.kind == 'GeomFS':
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
        elif qlab.kind == 'GeomMS':
            lnk = 718
            key = '              New orientation in intermediate state'
            sub = 5
            def end(s): return s.startswith(' ------')
            fmt = r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
        elif qlab.kind == 'ExGeom':
            lnk = -718
            key = ' Extrapolated geometry'
            sub = 4
            def end(s): return s.startswith(' ------')
            fmt = r'^\s+\w+(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
            num = 0
    elif qlab.label == 'vptdat':
        if qlab.kind == 'CICoef':
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
    elif qlab.label == 'vtrans':
        if qlab.kind == 'RR':
            lnk = (718, 718)
            key = ('                 Information on Transitions',
                   '                 Information on Transitions')
            sub = (2, 2)
            end = (lambda s: s.startswith('     ====='),
                   lambda s: s.startswith('     ====='))
            fmt = (r'^\s+Energy = \s*-?\d+\.\d+ cm.-1:\s+'
                   + r'(?P<val>\|.+ ->\s+\|.+)\s*$',
                   r'^\s+-> Omega =\s*(?P<val>\d+\.\d+) cm.-1$')
            num = (0, 0)
        elif qlab.kind == 'SOS':
            raise NotImplementedError()
        else:
            if qlab.level == 'H':
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
            elif qlab.level == 'A':
                lnk = 717
                key = ' NOTE: Transition energies are given with'
                sub = 8
                def end(s): return s.startswith('     =====')
                # fmt = r'^\s+\w?\s+(?P<val>(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'\
                #       + r'(?:\w+)?\s+(?:\s*-?\d+\.\d+|\*+\s+){4,5}.*\s*$'
                fmt = r'^\s+\w?\s+(?P<val>(?:\s*\d+\(\d+\)\s+(?:\w+)?)' \
                    + r'{1,3}|\d+)\s+(?:\s*-?\d+\.\d+|\*+\s+){4,5}.*\s*$'
                num = 0
            else:
                raise NotImplementedError()
    elif qlab.label == 'vlevel':
        if qlab.kind == 'RR':
            lnk = (718, 718)
            key = ('                 Information on Transitions',
                   '                 Information on Transitions')
            sub = (2, 2)
            end = (lambda s: s.startswith('     ====='),
                   lambda s: s.startswith('     ====='))
            fmt = (r'^\s+Energy = \s*(?P<val>-?\d+\.\d+) cm.-1:\s+\|.+$',
                   r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1$')
            num = (0, 0)
        elif qlab.kind == 'SOS':
            raise NotImplementedError()
        else:
            if qlab.level == 'H':
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
                       + r'(?P<val>-?\d+\.\d+|\*+)(?:\s+-?\d+\.\d+|\*+){4}'
                       + r'.*\s*$')
                num = (0, -1, 0)
            elif qlab.level == 'A':
                lnk = 717
                key = ' NOTE: Transition energies are given with'
                sub = 8
                def end(s): return s.startswith('     =====')
                fmt = r'^\s+\w?\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+(?:\w+)?' \
                    + r'\s+(?:-?\d+\.\d+|\*+)?\s+(?P<val>-?\d+\.\d+|\*+)' \
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
        if isinstance(qlab.rstate, tuple):
            Si, Sf = qlab.rstate
            if qlab.label == 1:
                if qlab.derord == 0:
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
                else:
                    raise NotImplementedError()
            elif qlab.label == 'dipstr':
                if Si == 0:
                    lnk = 0
                    key = ' Ground to excited state transition electric ' \
                        + 'dipole moments (Au):'
                    if isinstance(Sf, int):
                        sub = Sf + 1
                        def end(_s): return True
                    elif Sf == 'a':
                        sub = 1
                        def end(s): return s[2] != ' '
                    else:
                        raise ValueError('Unsupported final state')
                    fmt = r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+' \
                        + r'(?P<val>-?\d+\.\d+)\s+-?\d+\.\d+\s*$'
                    num = -1
                else:
                    raise NotImplementedError()
            elif qlab.label == 'rotstr':
                if Si == 0:
                    lnk = 0
                    key = ' Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm'\
                        + '/Gauss)'
                    if isinstance(Sf, int):
                        sub = Sf+1
                        def end(_s): return True
                    elif Sf == 'a':
                        sub = 1
                        def end(s): return s[2] != ' '
                    else:
                        raise ValueError('Unsupported final state')
                    if qlab.kind == 'vel':
                        num = 0
                        fmt = r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+' \
                            + r'(?P<val>-?\d+\.\d+)\s+-?\d+\.\d+\s*$'
                    else:
                        num = -1
                        fmt = r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+' \
                            + r'(?P<val>-?\d+\.\d+)\s*$'
                else:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()
            #     if qlab.label == 1 and dord == 0:
            #         keywords = ['ETran scalars', 'SCF Energy']
        else:
            if qlab.rstate == 'c':
                lnk1 = [lnk0]
                key1 = [key0]
                sub1 = [sub0]
                end1 = [end0]
                fmt1 = [fmt0]
                num1 = [num0]
            elif isinstance(qlab.rstate, int):
                lnk1 = []
                key1 = []
                sub1 = []
                end1 = []
                fmt1 = []
                num1 = []
            else:
                raise NotImplementedError()
            # Add quantity specific subgroups
            if qlab.label in ('ramact', 'roaact') and qlab.kind != 'static':
                # Information on setup and incident frequency
                # -- Harmonic level
                if qlab.kind != 'static' and qlab.level == 'H':
                    lnk1.append(716)
                    key1.append(' and normal coordinates:')
                    sub1.append(1)
                    end1.append(lambda s: not s.startswith(' Incident light:'))
                    num1.append(0)
                    fmt1.append(
                        r'^\s+Incident light \(\S+\):\s+(?P<val>-?\d.*)\s*$')
                # -- Anharmonic level
                lnk1.append(-717)
                if qlab.kind == 'dynamic':
                    key1.append(
                        '          Anharmonic Raman Spectroscopy (Dynamic)')
                    end1.append(
                        lambda s: s.startswith('     =====')
                        or s.startswith(' GradGrad'))
                    sub1.append(2)
                else:
                    fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                    key1.append(fmt.format(qlab.kind))
                    end1.append(
                        lambda s: s.startswith('     =====')
                        or s.startswith(' GradGrad')
                        or s.startswith(' ## INCIDENT WAVENUMBER:'))
                    sub1.append(0)
                fmt1.append(r'^ ##+ (?:INCIDENT WAVENUMBER|MEASUREMENT '
                            + r'INFORMATION):\s+(?P<val>\S+).*\s+##+\s*$')
                num1.append(0)
            if qlab.label == 1:
                if qlab.derord == 0:
                    if qlab.rstate == 'c':
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
                        if qlab.kind == 'first':
                            num1.extend((0, 0))
                        elif qlab.kind == 'last':
                            num1.extend((-1, -1))
                        else:
                            num1.extend((1, 1))
                    elif type(qlab.rstate) is int:
                        raise NotImplementedError()
                elif qlab.derord == 1:
                    raise NotImplementedError()
                elif qlab.derord == 2:
                    if qlab.rstate == 'c':
                        lnk1.append(-716)
                        key1.append(
                            ' Force constants in Cartesian coordinates:')
                        sub1.append(1)
                        end1.append(
                            lambda s:
                                not re.match(r'^\s+\d+\s+-?\d', s))
                        fmt1.append(
                            r'^\s+(?P<val>\d+(?:\s+\d+|'
                            + r'\s+-?\d+\.\d+D?[-+]\d+){1,5})\s*$')
                        num1.append(0)
                    elif type(qlab.rstate) is int:
                        raise NotImplementedError()
            elif qlab.label == 101:
                if qlab.derord == 0:
                    if qlab.rstate == 'c':
                        lnk1.append(601)
                        key1.append(' Dipole moment (field-independent basis,'
                                    + ' Debye):')
                        sub1.append(1)
                        end1.append(lambda s: True)
                        fmt1.append(r'^\s+(?P<val>(?:\s*[XYZ]=\s+-?\d+\.\d+)'
                                    + r'{3}).*$')
                        num1.append(0)
            elif qlab.label == 1300:
                if qlab.level == 'VE':
                    if qlab.rstate == 'c':
                        lnk1.append(718)
                        key1.append(
                            '                 Information on Transitions'
                        )
                        sub1.append(2)
                        end1.append(lambda s: s.startswith('     ====='))
                        fmt1.append(
                            r'^\s+-> Omega = \s*(?P<val>\d+\.\d+ cm.-1)$')
                        num1.append(0)
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()
            elif qlab.label == 1301:
                if qlab.level == 'VE':
                    if qlab.derord == 0:
                        if qlab.rstate == 'c':
                            lnk1.extend((718, 718))
                            key1.extend((
                                '     ELECTRIC DIPOLE-ELECTRIC DIPOLE',
                                '                 Information on Transitions'
                            ))
                            sub1.extend((2, 2))
                            end1.extend((
                                lambda s: '--' in s,
                                lambda s: s.startswith('     =====')
                            ))
                            fmt1.extend((
                                r'^\s+(?P<val>[XYZ]-[XYZ]'
                                + r'(?:\s+-?\d+\.\d+E?[-+]\d+){2}\s*$)',
                                r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1$'
                            ))
                            num1.extend((1, 0))
                        else:
                            raise NotImplementedError()
                    else:
                        msg = 'Tensor derivatives not available from vRR.'
                        raise ParseKeyError(msg)
                else:
                    raise NotImplementedError()
            elif qlab.label == 1302:
                if qlab.level == 'VE':
                    if qlab.derord == 0:
                        if qlab.rstate == 'c':
                            lnk1.extend((718, 718))
                            key1.extend((
                                '     ELECTRIC DIPOLE-MAGNETIC DIPOLE',
                                '                 Information on Transitions'
                            ))
                            sub1.extend((2, 2))
                            end1.extend((
                                lambda s: '--' in s,
                                lambda s: s.startswith('     =====')
                            ))
                            fmt1.extend((
                                r'^\s+(?P<val>[XYZ]-[XYZ]'
                                + r'(?:\s+-?\d+\.\d+E?[-+]\d+){2}\s*$)',
                                r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1$'
                            ))
                            num1.extend((1, 0))
                        else:
                            raise NotImplementedError()
                    else:
                        msg = 'Tensor derivatives not available from vRR.'
                        raise ParseKeyError(msg)
                else:
                    raise NotImplementedError()
            elif qlab.label == 1303:
                if qlab.level == 'VE':
                    if qlab.derord == 0:
                        if qlab.rstate == 'c':
                            lnk1.extend((718, 718))
                            key1.extend((
                                '     MAGNETIC DIPOLE-ELECTRIC DIPOLE',
                                '                 Information on Transitions'
                            ))
                            sub1.extend((2, 2))
                            end1.extend((
                                lambda s: '--' in s,
                                lambda s: s.startswith('     =====')
                            ))
                            fmt1.extend((
                                r'^\s+(?P<val>[XYZ]-[XYZ]'
                                + r'(?:\s+-?\d+\.\d+E?[-+]\d+){2}\s*$)',
                                r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1$'
                            ))
                            num1.extend((1, 0))
                        else:
                            raise NotImplementedError()
                    else:
                        msg = 'Tensor derivatives not available from vRR.'
                        raise ParseKeyError(msg)
                else:
                    raise NotImplementedError()
            elif qlab.label == 1304:
                if qlab.level == 'VE':
                    if qlab.derord == 0:
                        if qlab.rstate == 'c':
                            lnk1.extend((718, 718))
                            key1.extend((
                                '   ELECTRIC DIPOLE-ELECTRIC QUADRUPOLE',
                                '                 Information on Transitions'
                            ))
                            sub1.extend((2, 2))
                            end1.extend((
                                lambda s: '--' in s,
                                lambda s: s.startswith('     =====')
                            ))
                            fmt1.extend((
                                r'^\s+(?P<val>[XYZ]-[XYZ]{2}'
                                + r'(?:\s+-?\d+\.\d+E?[-+]\d+){2}\s*$)',
                                r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1$'
                            ))
                            num1.extend((1, 0))
                        else:
                            raise NotImplementedError()
                    else:
                        msg = 'Tensor derivatives not available from vRR.'
                        raise ParseKeyError(msg)
                else:
                    raise NotImplementedError()
            elif qlab.label == 1305:
                if qlab.level == 'VE':
                    if qlab.derord == 0:
                        if qlab.rstate == 'c':
                            lnk1.extend((718, 718))
                            key1.extend((
                                '   ELECTRIC QUADRUPOLE-ELECTRIC DIPOLE',
                                '                 Information on Transitions'
                            ))
                            sub1.extend((2, 2))
                            end1.extend((
                                lambda s: '--' in s,
                                lambda s: s.startswith('     =====')
                            ))
                            fmt1.extend((
                                r'^\s+(?P<val>[XYZ]{2}-[XYZ]'
                                + r'(?:\s+-?\d+\.\d+E?[-+]\d+){2}\s*$)',
                                r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1$'
                            ))
                            num1.extend((1, 0))
                        else:
                            raise NotImplementedError()
                    else:
                        msg = 'Tensor derivatives not available from vRR.'
                        raise ParseKeyError(msg)
                else:
                    raise NotImplementedError()
            elif qlab.label == 'dipstr':
                if qlab.level == 'H':
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
                elif qlab.level == 'A':
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
            elif qlab.label == 'rotstr':
                if qlab.level == 'H':
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
                elif qlab.level == 'A':
                    lnk1.append(717)
                    key1.append('        Rotational strengths (RS)')
                    sub1.append(3)
                    end1.append(lambda s: s.startswith('     ====='))
                    fmt1.append(r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'
                                + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$')
                    num1.append(0)
                else:
                    raise NotImplementedError()
            elif qlab.label == 'ramact':
                if qlab.level == 'H':
                    lnk1.extend([-716, 716])
                    key1.extend([' and normal coordinates:',
                                 ' and normal coordinates:'])
                    end1.extend([
                        lambda s: s.startswith(' Harmonic frequencies'),
                        lambda s: s.startswith(' - Thermochemistry')])
                    sub1.extend([1, 1])
                    num1.extend([0, -1])
                    if qlab.kind == 'static':
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
                    if qlab.kind in ('static', 'dynamic'):
                        fmt = 'Anharmonic Raman Spectroscopy ({})'
                        txt = fmt.format(qlab.kind.capitalize())
                        key1.append('{:>49s}'.format(txt))
                        end1.append(
                            lambda s: s.startswith('     =====')
                            or s.startswith(' GradGrad'))
                        sub1.append(2)
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qlab.kind))
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
                elif qlab.level == 'A':
                    lnk1.append(717)
                    if qlab.kind in ('static', 'dynamic'):
                        fmt = 'Anharmonic Raman Spectroscopy ({})'
                        txt = fmt.format(qlab.kind.capitalize())
                        key1.append('{:>49s}'.format(txt))
                        end1.append(
                            lambda s: s.startswith('     =====')
                            or s.startswith(' GradGrad'))
                        sub1.append(2)
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qlab.kind))
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
            elif qlab.label == 'roaact':
                if qlab.level == 'H':
                    lnk1.extend([716, -717])
                    key1.append(' and normal coordinates:')
                    end1.append(lambda s: s.startswith(' - Thermochemistry'))
                    if qlab.kind == 'dynamic':
                        key1.append(
                            '             Anharmonic Raman Optical Activity')
                        end1.append(
                            lambda s:
                                s.strip() == 'Dimensionless circular '
                                + 'intensity difference (CID)')
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qlab.kind))
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
                elif qlab.level == 'A':
                    lnk1.append(717)
                    if qlab.kind == 'dynamic':
                        key1.append(
                            '             Anharmonic Raman Optical Activity')
                        end1.append(
                            lambda s:
                                s.strip() == 'Dimensionless circular '
                                + 'intensity difference (CID)')
                    else:
                        fmt = ' ## INCIDENT WAVENUMBER:{:>15s} CM^-1 ##'
                        key1.append(fmt.format(qlab.kind))
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
        #     elif qlab.label == 50:
        #         raise NotImplementedError()
        #     elif qlab.label == 91:
        #         raise NotImplementedError()
        #     elif qlab.label == 92:
        #         raise NotImplementedError()
        #     elif qlab.label == 93:
        #         raise NotImplementedError()
        #     elif qlab.label == 102:
        #         raise NotImplementedError()
        #     elif qlab.label == 103:
        #         raise NotImplementedError()
        #     elif qlab.label == 104:
        #         raise NotImplementedError()
        #     elif qlab.label == 105:
        #         raise NotImplementedError()
        #     elif qlab.label == 106:
        #         raise NotImplementedError()
        #     elif qlab.label == 107:
        #         raise NotImplementedError()
        #     elif qlab.label == 201:
        #         raise NotImplementedError()
        #     elif qlab.label == 202:
        #         raise NotImplementedError()
        #     elif qlab.label == 203:
        #         raise NotImplementedError()
        #     elif qlab.label == 204:
        #         raise NotImplementedError()
        #     elif qlab.label == 205:
        #         raise NotImplementedError()
        #     elif qlab.label == 206:
        #         raise NotImplementedError()
        #     elif qlab.label == 207:
        #         raise NotImplementedError()
        #     elif qlab.label == 208:
        #         raise NotImplementedError()
        #     elif qlab.label == 209:
        #         raise NotImplementedError()
        #     elif qlab.label == 300:
        #         raise NotImplementedError()
        #     elif qlab.label == 301:
        #         raise NotImplementedError()
        #     elif qlab.label == 302:
        #         raise NotImplementedError()
        #     elif qlab.label == 303:
        #         raise NotImplementedError()
        #     elif qlab.label == 304:
        #         raise NotImplementedError()
        #     elif qlab.label == 305:
        #         raise NotImplementedError()
        #     elif qlab.label == 306:
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
    """Parse data arrays to extract specific quantity.

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
    dobjs = {}
    all_tags = []
    for qkey, qlabel in qdict.items():
        msg_noqty = f'Missing quantity "{qlabel}" in file'
        all_tags.append(qlabel.label)
        first, last = key2blocks[qkey]
        # Basic Check: property available
        # -----------
        # Check if some data extracted
        num = 0
        if first == last:
            iref = first
            num = len(datablocks[iref])
        else:
            iref = -1
            # Check if transition information may be present:
            if (qlabel.rstate == 'c'
                    and 'Excited State' in ''.join(datablocks[first])):
                realfirst = first + 1
            else:
                realfirst = first
            for i in range(realfirst, last+1):
                if num == 0 and ndatablock[i] > 0 and datablocks[i]:
                    iref = i  # Ignore first, null indexes
                    num = len(datablocks[i])
        if num == 0 and not empty_cases_ok(qlabel.label, qlabel.kind):
            if raise_error:
                raise ParseKeyError(msg_noqty)
            dobjs[qkey] = None
            continue
        dobjs[qkey] = QData(qlabel)
        # Basic Properties/Quantities
        # ---------------------------
        if qlabel.label == 'route':
            fmt = '{:d}{:02d}'
            data = []
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
                    data.append((fulllink, num+1, iops, toline))
            dobjs[qkey].set(data=data)
        elif qlabel.label == 'natoms':
            dobjs[qkey].set(data=int(datablocks[iref][0]))
        elif qlabel.label == 'atcrd' or qlabel.label == 2:
            # By default, we choose the standard orientation if present
            if datablocks[last]:
                i = last
            else:
                i = first
            data = []
            if qlabel.kind in ('all', 'scan'):
                for block in datablocks[i]:
                    _block = []
                    for line in block:
                        _block.append([float(item)*__ang2au
                                       for item in line.split()])
                    data.append(_block)
                num = len(datablocks[i])
            else:
                for line in datablocks[i]:
                    data.append([float(item)*__ang2au
                                 for item in line.split()])
                num = 1
            dobjs[qkey].set(data=data)
            dobjs[qkey].add_field('ngeoms', value=num)
        elif qlabel.label == 'atmas':
            data = []
            for line in datablocks[iref]:
                data.extend(float(item) for item in line.split())
            dobjs[qkey].set(data=data)
        elif qlabel.label in ('atnum',):
            data = []
            for line in datablocks[iref]:
                data.extend(int(item) for item in line.split())
            dobjs[qkey].set(data=data)
        elif qlabel.label == 'swopt':
            # The last line are the final dashes
            dobjs[qkey].set(data=' '.join(datablocks[iref][:-1]))
        elif qlabel.label == 'molsym':
            raise NotImplementedError()
        elif qlabel.label == 'swver':
            txt = r'\s*Gaussian (\w+):\s+(\w+)-(\w{3})Rev([\w.+]+) {1,2}' \
                + r'(\d+-\w{3}-\d{4})\s*'
            pattern = re.compile(txt)
            res = re.match(pattern, ''.join(datablocks[iref])).groups()
            dobjs[qkey].set(data={'major': res[2], 'minor': res[3],
                                  'system': res[1], 'release': res[4]})
        # Vibrational Information
        # -----------------------
        # Technically state should be checked but considered irrelevant.
        elif qlabel.label == 'nvib':
            i = 0
            for line in datablocks[iref]:
                i += len(line.split())

            dobjs[qkey].set(data=i)
        elif qlabel.label == 'hessvec':
            for i in range(first, last+1):
                if datablocks[i]:
                    iref = i
                    break
            else:
                raise ParseKeyError(msg_noqty)
            dobjs[qkey].set(dtype='L.M^{-1/2}')
            data = []
            # We analyze the first line to find if HPModes or normal
            res = datablocks[iref][0].split()[-1]
            if len(res.split('.')[1]) == 2:
                ncols = 3  # maximum number of modes per block
                ioff = -ncols  # Started at -3 to offset increment in 1st block
                for line in datablocks[iref]:
                    cols = line.split()
                    nmodes = (len(cols)-2)//3
                    if cols[0] == '1':
                        ioff += ncols
                        data.extend(
                            [[] for _ in range(nmodes)])
                    for i in range(nmodes):
                        data[ioff+i].extend(
                            [float(item) for item in cols[2+i*3:2+(i+1)*3]])
            else:
                ncols = 5
                ioff = -ncols
                for line in datablocks[iref]:
                    cols = line.split()
                    nmodes = len(cols) - 3
                    if cols[0] == '1' and cols[1] == '1':
                        ioff += ncols
                        data.extend(
                            [[] for _ in range(nmodes)])
                    for i in range(nmodes):
                        data[ioff+i].append(float(cols[3+i]))
            dobjs[qkey].set(data=data)
        elif qlabel.label == 'hessdat':
            for i in range(first, last+1):
                if datablocks[i]:
                    iref = i
                    break
            else:
                raise ParseKeyError(msg_noqty)
            if qlabel.kind == 'freq':
                dobjs[qkey].set(unit='cm-1')
            elif qlabel.kind == 'redmas':
                dobjs[qkey].set(unit='amu')
            else:
                raise NotImplementedError('Unknown subopt for HessDat')
            data = []
            i = 0
            for line in datablocks[iref]:
                data.extend(
                    [float(item) if '*' not in item else float('inf')
                     for item in line.split()])
            dobjs[qkey].set(data=data)
        # General Spectroscopy
        # --------------------
        elif qlabel.label == 'intens':
            if qlabel.level == 'H':
                for i in range(last, first-1, -1):
                    if datablocks[i]:
                        iref = i
                        break
                else:
                    raise ParseKeyError(msg_noqty)
                if qlabel.kind == 'IR':
                    dobjs[qkey].set(unit='II:km.mol-1')
                else:
                    dobjs[qkey].set(unit='II:N/A')
                i = 0
                data = {}
                for line in datablocks[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[i] = float(col)
                        except ValueError:
                            data[i] = float('inf')
            elif qlabel.level == 'A':
                if qlabel.kind == 'IR':
                    dobjs[qkey].set(unit='II:km.mol-1')
                else:
                    dobjs[qkey].set(unit='II:N/A')
                i = 0
                for line in datablocks[iref]:
                    i += 1
                    try:
                        data[i] = float(line)
                    except ValueError:
                        data[i] = float('inf')
            else:
                raise NotImplementedError()
            dobjs[qkey].set(data=data)

        # Vibronic Information
        # --------------------
        elif qlabel.label == 'fcdat':
            if qlabel.kind == 'SimInf':
                def getval(s):
                    return s.split(':', maxsplit=1)[-1].strip().title()
                for line in datablocks[iref]:
                    if 'Temperature' in line:
                        if 'not' in line:
                            val = None
                        else:
                            val = float(line.split()[-1][:-1])
                        dobjs[qkey].add_field('temp', value=val)
                    elif 'framework' in line:
                        dobjs[qkey].add_field('frame', value=getval(line))
                    elif 'Spectroscopy' in line:
                        dobjs[qkey].add_field('spec', value=getval(line))
                    elif 'Model' in line:
                        dobjs[qkey].add_field('model', value=getval(line))
                    elif 'electronic transition moment' in line:
                        dobjs[qkey].add_field(
                            'tmom',
                            value=line.split(':', maxsplit=1)[-1].strip())
            elif qlabel.kind == 'JMat':
                if datablocks[last]:
                    i = last
                else:
                    i = first
                if 'diagonal' in datablocks[i]:
                    dobjs[qkey].set(data=[])
                else:
                    data = []
                    N = 0
                    for line in datablocks[i]:
                        cols = line.split()
                        irow = int(cols[0]) - 1
                        if irow == N:
                            data.append([])
                            N += 1
                        data[irow].extend([float(item.replace('D', 'e'))
                                           for item in cols[1:]])
                    dobjs[qkey].set(data=data)
            elif qlabel.kind == 'JMatF':
                if 'diagonal' in datablocks[iref]:
                    dobjs[qkey].set(data=[])
                else:
                    data = []
                    N = 0
                    for line in datablocks[iref]:
                        cols = line.split()
                        irow = int(cols[0]) - 1
                        if irow == N:
                            data.append([])
                            N += 1
                        data[irow].extend([float(item.replace('D', 'e'))
                                           for item in cols[1:]])
                    dobjs[qkey].set(data=data)
            elif qlabel.kind == 'KVec':
                data = []
                for line in datablocks[iref]:
                    data.append(float(line.replace('D', 'e')))
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'SRAMat':
                data = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data.append([])
                        N += 1
                    data[irow].extend([float(item.replace('D', 'e'))
                                       for item in cols[1:]])
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'SRBVec':
                data = []
                for line in datablocks[iref]:
                    data.append(float(line.replace('D', 'e')))
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'SRCMat':
                data = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data.append([])
                        N += 1
                    data[irow].extend([float(item.replace('D', 'e'))
                                       for item in cols[1:]])
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'SRDVec':
                data = []
                for line in datablocks[iref]:
                    data.append(float(line.replace('D', 'e')))
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'SREMat':
                data = []
                N = 0
                for line in datablocks[iref]:
                    cols = line.split()
                    irow = int(cols[0]) - 1
                    if irow == N:
                        data.append([])
                        N += 1
                    data[irow].extend([float(item.replace('D', 'e'))
                                       for item in cols[1:]])
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'Spec':
                # Look for last blocks first, which should contain all blocks:
                # Note: For the vibronic spectrum, X is necessarily the same
                #       over multiple blocks.  We take this for granted to
                #       simplify the processing and the output.
                discard = []
                if last != first:
                    nblocks = 0
                    nyaxes = 0
                    for i, bloc in enumerate(datablocks[last]):
                        lbloc = len(bloc)
                        if lbloc > 1:
                            # 1 for the x axis
                            nblocks += 1
                            nyaxes += len(bloc[0].split()) - 1
                        else:
                            # Incorrect bloc
                            discard.append(i)
                else:
                    nblocks = 1
                    nyaxes = 0
                if nblocks == 0:
                    msg = 'Inconsistency in spectral data. ' \
                            + 'This should not happen!'
                    raise IndexError(msg)
                if discard:
                    discard.reverse()
                    for i in discard:
                        del datablocks[last][i]
                data = {}
                # First block contains the full initial legend block
                # Necessary for the different parameters
                iref = first
                data['x'] = []
                if nyaxes == 0:
                    nyaxes = len(datablocks[iref][0].split()) - 1
                if nyaxes == 1:
                    yfmt = 'y'
                    data['y'] = []
                else:
                    yfmt = f'y{{idy:0{len(str(nyaxes))}d}}'
                    for i in range(nyaxes):
                        data[yfmt.format(idy=i)] = []
                # If multiple blocks, the reading with first parsing method
                #   is wrong since it combines all blocks together
                # We fix it by only considering the last one which should be
                #   always right.
                iref = last
                yoffset = 0
                for bloc in range(nblocks):
                    yax = {}
                    if nblocks > 1:
                        data['x'].append([])
                        xax = data['x'][-1]
                        for i in range(nyaxes):
                            y = yfmt.format(idy=i+1)
                            data[y].append([])
                            yax[y] = data[y][-1]
                    else:
                        xax = data['x']
                        for i in range(nyaxes):
                            y = yfmt.format(idy=i+1)
                            yax[y] = data[y]
                    for line in datablocks[iref][bloc]:
                        cols = [float(item.replace('D', 'e'))
                                for item in line.split()]
                        if block == 0:
                            xax.append(cols[0])
                        for i, item in enumerate(cols[1:]):
                            yax[yfmt.format(idy=yoffset+i+1)].append(item)
                    yoffset += len(cols) - 1
                for key, val in data.items():
                    dobjs[qkey].add_field(key, value=val)
            elif qlabel.kind == 'SpcPar':
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
                        del datablocks[last][i]
                # First block contains the full initial legend block
                # Necessary for the different parameters
                iref = first
                i = 0
                while datablocks[iref][i].strip() != 'Legend:':
                    i += 1
                    if i >= len(datablocks[iref]):
                        raise IndexError('Legend not found')
                if 'No' in datablocks[iref][i-1]:
                    dobjs[qkey].add_field('func', value='stick')
                    dobjs[qkey].add_field('hwhm', value=None)
                elif 'broadening' in datablocks[iref][i-2]:
                    func = datablocks[iref][i-2].split()[-3].lower()
                    hwhm = float(datablocks[iref][i-1].split()[-2])
                    dobjs[qkey].add_field('func', value=func)
                    dobjs[qkey].add_field('hwhm', value=hwhm)
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
                    yfmt = f'y{{idy:0{len(str(nyaxes))}d}}'
                yoffset = 0
                nyi = 0
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
                        dobjs[qkey].add_field('unitx', value=txt)
                    elif key.strip() == 'Intensity':
                        _key = 'I'
                        txt = title.rsplit('in ', maxsplit=1)[1].\
                            replace(')', '').strip()
                        if dobjs[qkey].get('func') == 'stick':
                            _desc = 'II:'
                            # Fix a stupidity in the unit in some versions
                            #   for the stick spectrum
                            if txt == 'dm^3.mol^-1.cm^-1':
                                txt = 'dm^3.mol^-1.cm^-2'
                        else:
                            _desc = 'I:'
                        dobjs[qkey].add_field('unity', value=_desc + txt)
                    else:
                        try:
                            idy = int(key.strip()[0]) - 1 + yoffset
                        except ValueError as err:
                            raise IndexError(
                                f'Unrecognized key in spc leg.: {key}') \
                                    from err
                        _key = yfmt.format(idy=idy)
                        nyi += 1
                    if nblocks == 1:
                        dobjs[qkey].add_field(_key, value=title.strip())
                    else:
                        if _key[0] == 'I':
                            dobjs[qkey].add_field(_key, value=title.strip())
                        else:
                            dobjs[qkey].add_field(_key, value=title.strip())
                # More than 1 block, read the new ones
                for bloc in range(1, nblocks):
                    nyi = 0
                    for line in datablocks[last][bloc][1:]:
                        res = line.split(':', 1)
                        if len(res) == 2:
                            key, title = res
                            if 'col.' in key and key.strip() != '1st col.':
                                try:
                                    idy = int(key.strip()[0]) - 1 + yoffset
                                except ValueError as err:
                                    msg = f'Unknown key in spc leg.: {key}'
                                    raise IndexError(msg) from err
                                _key = yfmt.format(idy=idy)
                                dobjs[qkey].add_field(_key,
                                                      value=title.strip())
                        elif 'intensity' not in res[0]:
                            raise NotImplementedError('Unrecognized legend')
                    yoffset += nyi
            elif qlabel.kind == 'Conv':
                raise NotImplementedError()
            elif qlabel.kind == 'Assign':
                data = {'T': [], 'E': [], 'I': []}
                if qlabel.label == 'DipStr':
                    qty = 'DS'
                elif qlabel.label == 'RotStr':
                    qty = 'RS'
                else:
                    msg = 'Unrecognized spectroscopy-specific quantity'
                    raise IndexError(msg)
                dobjs[qkey].add_field('other', value=qty)
                data[qty] = []
                for l1, l2 in zip(datablocks[first], datablocks[last]):
                    txt_E, txt_T = l1.split(':')
                    data['E'].append(float(txt_E.split()[0]))
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
                    data['T'].append((tuple(trans[0]), tuple(trans[1])))
                    line = l2.strip(')').split()
                    data['I'].append(float(line[0]))
                    data[qty].append(float(line[-1]))
                for key, val in data.items():
                    dobjs[qkey].add_field(key, value=val)
            elif qlabel.kind == 'GeomIS':
                data = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data.append([float(item)*__ang2au
                                 for item in line.split()])
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'GeomFS':
                data = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data.append([float(item)*__ang2au
                                 for item in line.split()])
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'ExGeom':
                data = []
                # By default, we choose the standard orientation if present
                for line in datablocks[iref]:
                    data.append([float(item)*__ang2au
                                 for item in line.split()])
                dobjs[qkey].set(data=data)
            elif qlabel.kind == 'E(0-0)':
                if len(set(datablocks[iref])) > 1:
                    msg = 'Excessive information on 0-0 energy.'
                    raise ParseKeyError(msg)
                val, unit = datablocks[iref][0].split()
                if re.match(r'cm\^.?-1.?\s*', unit, re.I):
                    dobjs[qkey].set(unit='cm^-1')
                else:
                    dobjs[qkey].set(unit='???')
                dobjs[qkey].set(data=float(val))
            elif qlabel.kind == 'RedDim':
                if len(datablocks[iref][0].split()) % 3 != 0:
                    raise ParseKeyError('Expected reddim structure: num = num')
                nstates = len(datablocks[iref][0].split())//3
                data = {}
                for i in range(nstates):
                    data['state{}'.format(i+1)] = {}
                for line in datablocks[iref]:
                    cols = line.split()
                    for i in range(nstates):
                        data['state{}'.format(i+1)][int(cols[i*3])] = \
                            int(cols[i*3+2])
                for key, val in data.items():
                    dobjs[qkey].add_field(key, value=val)
            else:
                raise NotImplementedError('Unknown option for FCDat')
        # Vibrational transitions
        # -----------------------
        elif qlabel.label == 'vlevel':
            if qlabel.kind == 'RR':
                if len(datablocks[last]) != len(datablocks[last-1]):
                    msg = 'Incident frequencies data and transition energies' \
                        + ' do not match.'
                    raise ParseKeyError(msg)
                dobjs[qkey].set(unit='cm-1')
                counts = {}
                data = {}
                for incfrq, trans in zip(datablocks[last], datablocks[last-1]):
                    if incfrq not in data:
                        data[incfrq] = {}
                        counts[incfrq] = 1
                    else:
                        counts[incfrq] += 1
                    data[incfrq][counts[incfrq]] = float(trans)
                for key, val in data.items():
                    dobjs[qkey].add_field(key, value=val)
            else:
                if qlabel.level == 'H':
                    for i in range(last, first-1, -1):
                        if datablocks[i]:
                            iref = i
                            break
                    else:
                        raise ParseKeyError(msg_noqty)
                    dobjs[qkey].set(unit='cm-1')
                    i = 0
                    data = {}
                    for line in datablocks[iref]:
                        for col in line.strip().split():
                            i += 1
                            try:
                                data[i] = float(col)
                            except ValueError:
                                data[i] = float('inf')
                    dobjs[qkey].set(data=data)
                elif qlabel.level == 'A':
                    if datablocks[last]:
                        iref = last
                    else:
                        iref = first
                    dobjs[qkey].set(unit='cm-1')
                    i = 0
                    data = {}
                    for line in datablocks[iref]:
                        i += 1
                        try:
                            data[i] = float(line)
                        except ValueError:
                            data[i] = float('inf')
                    dobjs[qkey].set(data=data)
                else:
                    raise NotImplementedError()
        elif qlabel.label == 'vtrans':
            if qlabel.kind == 'RR':
                if len(datablocks[last]) != len(datablocks[last-1]):
                    msg = 'Incident frequencies data and transition data ' \
                        + 'do not match.'
                    raise ParseKeyError(msg)
                counts = {}
                data = {}
                for incfrq, trans in zip(datablocks[last], datablocks[last-1]):
                    if incfrq not in data:
                        data[incfrq] = {}
                        counts[incfrq] = 1
                    else:
                        counts[incfrq] += 1
                    svals = []
                    for i, sdat in enumerate(trans.split('->')):
                        sdesc = sdat.strip(' |>')
                        if sdesc == '0':
                            svals.append((0, 0))
                        else:
                            state = []
                            for osc in sdesc.split(','):
                                state.append(tuple([
                                    int(i) for i in osc.split('^')
                                ]))
                            svals.append(tuple(state))
                    data[incfrq][counts[incfrq]] = tuple(svals)
                dobjs[qkey].set(data=data)
            else:
                if qlabel.level == 'H':
                    for i in range(last, first-1, -1):
                        if datablocks[i]:
                            iref = i
                            break
                    else:
                        raise ParseKeyError(msg_noqty)
                    i = 0
                    data = {}
                    for line in datablocks[iref]:
                        cols = line.strip().split()
                        if cols[-1] in ('active', 'inactive', 'passive'):
                            del cols[-1]
                        for col in cols:
                            i += 1
                            res = col.split('(')
                            data[i] = (((0, 0), ), ((int(res[0]), 1), ))
                    dobjs[qkey].set(data=data)
                elif qlabel.level == 'A':
                    i = 0
                    data = {}
                    for line in datablocks[iref]:
                        i += 1
                        val = []
                        cols = line.strip().split()
                        if cols[-1] in ('active', 'inactive', 'passive'):
                            status = cols[-1]
                            del cols[-1]
                        else:
                            status = None
                        for col in cols:
                            res = col.split('(')
                            if len(res) == 1:
                                val.append((int(res[0]), 0))
                            else:
                                val.append((int(res[0]),
                                            int(res[1].replace(')', ''))))
                            data[i] = (((0, 0), ), tuple(val), status)
                    dobjs[qkey].set(data=data)
                else:
                    raise NotImplementedError()
        # Anharmonic Information
        # ----------------------
        elif qlabel.label == 'vptdat':
            if qlabel.kind == 'CICoef':
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
                data = {}
                i = 0
                for line in datablocks[first]:
                    if ':' in line:
                        key, dat = line.split(':')
                        i = idx(int(key))
                        if i > 0:
                            data[i] = []
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
                        data[i].append((float(coef), tuple(dat)))
                dobjs[qkey].set(data=data)
            else:
                raise NotImplementedError()
        # State(s)-dependent quantities
        # -----------------------------
        else:
            # Transition moments
            # ^^^^^^^^^^^^^^^^^^
            if isinstance(qlabel.rstate, tuple):
                Si, Sf = qlabel.rstate
                if qlabel.label == 1:
                    if qlabel.derord == 0:
                        if Si == 0:
                            if isinstance(Sf, int):
                                dobjs[qkey].set(data=float(datablocks[iref]))
                                dobjs[qkey].set(unit='eV')
                            else:
                                dobjs[qkey].set(unit='eV')
                                dobjs[qkey].set(
                                    data=[float(item)
                                          for item in datablocks[iref]])
                        else:
                            pass
                    else:
                        pass
                elif qlabel.label == 'dipstr':
                    if Si == 0:
                        if isinstance(Sf, int):
                            dobjs[qkey].set(data=float(datablocks[iref]))
                            dobjs[qkey].set(unit='DS:au')
                        else:
                            dobjs[qkey].set(unit='DS:au')
                            dobjs[qkey].set(
                                data=[float(item)
                                      for item in datablocks[iref]])
                    else:
                        pass
                elif qlabel.label == 'rotstr':
                    if Si == 0:
                        if isinstance(Sf, int):
                            dobjs[qkey].set(
                                data=float(datablocks[iref])*1.0e-40)
                            dobjs[qkey].set(unit='RS:esu^2.cm^2')
                        else:
                            dobjs[qkey].set(unit='RS:esu^2.cm^2')
                            dobjs[qkey].set(
                                data=[float(item) * 1.0e-40
                                      for item in datablocks[iref]])
                    else:
                        pass
                else:
                    pass
            # States-specific Quantities
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^
            else:
                if qlabel.rstate != 'c':
                    # we should check if the state is the right one.
                    raise NotImplementedError()
                if qlabel.label == 1:
                    if qlabel.derord == 0:
                        if qlabel.rstate == 'c':
                            def conv(s):
                                return float(s.split('=')[1].replace('D', 'e'))

                            fmt = re.compile(r'E\(?(.*?)\)?\s+')
                            N = 2 if datablocks[last] else 1
                            nblocks = len(datablocks[iref])
                            data = {}
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
                                data[tag] = val
                            dobjs[qkey].set(data=data[tag])
                            dobjs[qkey].set(unit='Eh')
                            for key, val in data.items():
                                dobjs[qkey].add_field(key, value=val)
                    elif qlabel.derord == 2:
                        if qlabel.rstate == 'c':
                            maxcols = 5
                            nbloc = 0
                            dobjs[qkey].set(unit='Eh.a0^{-2}')
                            dobjs[qkey].set(shape='lt')
                            data = {}
                            # store in linear form
                            for line in datablocks[iref]:
                                if '.' not in line:
                                    nbloc += 1
                                else:
                                    cols = line.split()
                                    i = int(cols[0])
                                    if nbloc == 1:
                                        data.extend(
                                            [float(item.replace('D', 'e'))
                                             for item in cols[1:]])
                                        if i > 5:
                                            data.extend(
                                                0.0 for _ in range(maxcols, i))
                                    else:
                                        ioff = i*(i-1)//2 + (nbloc-1)*maxcols
                                        ncols = len(cols) - 1
                                        data[ioff:ioff+ncols] \
                                            = [float(item.replace('D', 'e'))
                                               for item in cols[1:]]
                            dobjs[qkey].set(data=data)
                elif qlabel.label == 101:
                    if qlabel.derord == 0:
                        if qlabel.rstate == 'c':
                            val = [float(item)/phys_fact('au2Deb') for item in
                                   datablocks[last][0].split()[1::2]]
                            dobjs[qkey].set(data=val)
                            dobjs[qkey].set(unit='e.a0')
                        else:
                            raise NotImplementedError()
                    else:
                        raise NotImplementedError()
                elif qlabel.label == 1300:
                    unit = datablocks[iref][0].split()[-1]
                    dobjs[qkey].set(unit=unit.lower())
                    data = {'data': [], 'keys': []}
                    for line in datablocks[iref]:
                        key, unit = line.split()
                        if key not in data['keys']:
                            data['keys'].append(key)
                            data['data'].append(float(key))
                    dobjs[qkey].set(data=data['data'])
                    dobjs[qkey].add_field('keys', value=data['keys'])
                elif qlabel.label in (1301, 1302, 1303, 1304, 1305):
                    dobjs[qkey] = _parse_logdat_vtransprop(qlabel, datablocks,
                                                           first, last)
                elif qlabel.label == 'dipstr':
                    if qlabel.level == 'H':
                        for i in range(last, first-1, -1):
                            if datablocks[i]:
                                iref = i
                                break
                        else:
                            raise ParseKeyError(msg_noqty)
                        dobjs[qkey].set(unit='DS:esu^2.cm^2')
                        i = 0
                        data = {}
                        for line in datablocks[iref]:
                            for col in line.strip().split():
                                i += 1
                                try:
                                    data[i] = float(col)*1.0e-40
                                except ValueError:
                                    data[i] = float('inf')
                        dobjs[qkey].set(data=data)
                    elif qlabel.level == 'A':
                        dobjs[qkey].set(unit='DS:esu^2.cm^2')
                        i = 0
                        data = {}
                        for line in datablocks[iref]:
                            i += 1
                            try:
                                data[i] = float(line)*1.0e-40
                            except ValueError:
                                data[i] = float('inf')
                        dobjs[qkey].set(data=data)
                    else:
                        raise NotImplementedError()
                elif qlabel.label == 'rotstr':
                    if qlabel.level == 'H':
                        for i in range(last, first-1, -1):
                            if datablocks[i]:
                                iref = i
                                break
                        else:
                            raise ParseKeyError(msg_noqty)
                        dobjs[qkey].set(unit='RS:esu^2.cm^2')
                        i = 0
                        data = {}
                        for line in datablocks[iref]:
                            for col in line.strip().split():
                                i += 1
                                try:
                                    data[i] = float(col)*1.0e-44
                                except ValueError:
                                    data[i] = float('inf')
                        dobjs[qkey].set(data=data)
                    elif qlabel.level == 'A':
                        dobjs[qkey].set(unit='RS:esu^2.cm^2')
                        i = 0
                        data = {}
                        for line in datablocks[iref]:
                            i += 1
                            try:
                                data[i] = float(line)*1.0e-44
                            except ValueError:
                                data[i] = float('inf')
                        dobjs[qkey].set(data=data)
                    else:
                        raise NotImplementedError()
                elif qlabel.label == 'ramact':
                    dobjs[qkey] = _parse_logdat_ramact(qlabel, datablocks,
                                                       first, last)
                elif qlabel.label == 'roaact':
                    dobjs[qkey] = _parse_logdat_ramact(qlabel, datablocks,
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

    # final checks/corrections
    if not {'vlevel', 'vtrans'} - set(all_tags):
        args = {}
        for qkey, qlabel in qdict.items():
            if qlabel.label == 'vlevel' and qlabel.level == 'A':
                args['vlevel'] = qkey
            elif qlabel.label == 'vtrans' and qlabel.level == 'A':
                args['vtrans'] = qkey
        if 'vtrans' in args and 'vlevel' in args:
            dobjs[args['vtrans']], dobjs[args['vlevel']] = \
                __del_nonactive_modes(dobjs[args['vtrans']],
                                      dobjs[args['vlevel']])

    return dobjs


def get_data(dfobj: GLogIO,
             *qlabels: tp.Union[str, QLabel],
             error_noqty: bool = True,
             **keys4qlab) -> TypeQData:
    """Get data from a GLog file for each quantity label.

    Reads one or more full quantity labels from `qlab` and returns the
    corresponding data.

    Parameters
    ----------
    dfobj
        Gaussian output file as `GLogIO` object.

    *qlabels
        List of quantity labels, as strings to parse or QLabel objects.
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
    ArgumentError
        Unrecognized qlabel structure.
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
    if len(qlabels) == 0 and len(keys4qlab) == 0:
        return None
    # Build Keyword List
    # ------------------
    # Build full list of qlabels
    qty_dict, dupl_keys = parse_qlabels(qlabels, keys4qlab)
    # List of keywords
    keydata = []
    key2blocks = {}
    idata = 0
    for qkey, qlabel in qty_dict.items():
        # Label parsing
        # ^^^^^^^^^^^^^
        links, keys, skips, fmts, ends, nums = qlab_to_linkdata(qlabel)
        if isinstance(links, tuple):
            nblocks = len(links)
            key2blocks[qkey] = (idata, idata+nblocks-1)
            for i in range(nblocks):
                keydata.append((keys[i], links[i], skips[i],
                                re.compile(fmts[i]), nums[i], ends[i]))
            idata += nblocks
        else:
            key2blocks[qkey] = (idata, idata)
            keydata.append((keys, links, skips, re.compile(fmts), nums, ends))
            idata += 1
    # Data Extraction
    # ---------------
    # Use of set to remove redundant keywords
    ndata, datablocks = dfobj.read_data(*keydata)
    # Data Parsing
    # ------------
    gver = (dfobj.version['major'], dfobj.version['minor'])
    try:
        data = parse_data(qty_dict, key2blocks, ndata, datablocks, gver,
                          error_noqty)
    except (QuantityError, NotImplementedError) as err:
        raise QuantityError('Unsupported quantities') from err
    except (ParseKeyError, IndexError) as err:
        raise IndexError(f'Missing data in Gaussian log: {dfobj.filename}') \
            from err
    # Fill redundant keys
    if dupl_keys:
        for item, key in dupl_keys.items():
            data[item] = data[key]

    return data


def _parse_logdat_ramact(qlabel: QLabel,
                         datablocks: TypeDGLog,
                         first: int,
                         last: int,
                         ROA: bool = False) -> QData:
    """Parse GLog data for Raman activity.

    Sub-functions dedicated to parsing data stored in datablocks.

    Parameters
    ----------
    qlabel
        QLabel instance.
    datablock
        Datablocks generated from Gaussian logfile.
    first
        Index of first block of interest in `datablocks`.
    last
        Index of last block of interest in `datablocks`.
    ROA
        If true, data are related to ROA.

    Returns
    -------
    QData
        QData instance
    """
    dobj = QData(qlabel)
    data = {}
    for i in range(last, first-1, -1):
        if datablocks[i]:
            iref = i
            break
    else:
        raise ParseKeyError('Missing data for Raman/ROA activity in file')
    # Data in Gaussian log are multiplied by 10^4 for ROA, so we need to take
    #   this into account
    yfactor = 1.0e-4 if ROA else 1.0
    if qlabel.level == 'H':
        if iref == last:
            dobj.set(unit='ROA:Ang^6' if ROA else 'RA:Ang^6')
            an_blk = True
        else:
            dobj.set(unit='ROA:amu.Ang^4' if ROA else 'RA:amu.Ang^4')
            an_blk = False
        # We create the database based on the quantity
        if qlabel.kind == 'static':
            i = 0
            for line in datablocks[iref]:
                for col in line.strip().split():
                    i += 1
                    try:
                        data[i] = float(col)*yfactor
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
                except ValueError as err:
                    if incfrq is None:
                        msg = 'Wrong block structure for Raman/ROA ' \
                            + 'spectroscopy (incident freq/setup)'
                        raise ParseKeyError(msg) from err
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
                        d[i] = float(col)*yfactor
                    except ValueError:
                        d[i] = float('inf')
        else:
            incfreqs = [item for line in datablocks[first+1]
                        for item in line.split()]
            if qlabel.kind == 'dynamic':
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
                            d[ioff+i] = float(col)*yfactor
                        except ValueError:
                            d[ioff+i] = float('inf')
            else:
                if qlabel.kind in incfreqs:
                    ref_freq = incfreqs.index(qlabel.kind)
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
                                d[iref+i+1] = float(col)*yfactor
                            except ValueError:
                                d[iref+i+1] = float('inf')
                else:
                    raise ParseKeyError('Missing incident frequency')
    elif qlabel.level == 'A':
        dobj.set(unit='ROA:Ang^6' if ROA else 'RA:Ang^6')
        incfrq = None
        setups = []
        for item in datablocks[first+1]:
            try:
                _ = float(item)
                incfrq = item
                data[incfrq] = {}
            except ValueError as err:
                if incfrq is None:
                    msg = 'Wrong block structure for Raman/ROA spectroscopy ' \
                        + '(incident freq/setup)'
                    raise ParseKeyError(msg) from err
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
                    d[i] = float(col)*yfactor
                except ValueError:
                    d[i] = float('inf')
    else:
        raise NotImplementedError()
    dobj.set(data=data)
    return dobj


def _parse_logdat_vtransprop(qlabel: QLabel,
                             datablocks: TypeDGLog,
                             first: int,
                             last: int) -> tp.Dict[str, tp.Any]:
    """Parse GLog data for vib. transition moment.

    Sub-function dedicated to parsing data related to vibrational,
    transition moments of properties, stored in datablocks.

    Parameters
    ----------
    qlabel
        QLabel instance.
    datablock
        Datablocks generated from Gaussian logfile.
    first
        Index of first block of interest in `datablocks`.
    last
        Index of last block of interest in `datablocks`.

    Returns
    -------
    QData
        QData instance
    """
    if qlabel.label in range(1301, 1310):
        if len(datablocks[last]) != len(datablocks[last-1]):
            msg = 'Incident frequencies data and properties moments ' \
                + 'do not match.'
            raise ParseKeyError(msg)
    counts = {}
    dobj = QData(qlabel)
    data = {}
    dosym = False
    for incfrq, dtens in zip(datablocks[last], datablocks[last-1]):
        if incfrq not in data:
            data[incfrq] = {}
            counts[incfrq] = 1
        else:
            counts[incfrq] += 1
        if qlabel.label in (1301, 1302, 1303):
            tensor = [
                [None, None, None],
                [None, None, None],
                [None, None, None]
            ]
        elif qlabel.label in (1304, 1305):
            dosym = len(dtens) == 18
            tensor = [
                [[None, None, None],
                 [None, None, None],
                 [None, None, None]],
                [[None, None, None],
                 [None, None, None],
                 [None, None, None]],
                [[None, None, None],
                 [None, None, None],
                 [None, None, None]]
            ]
        else:
            raise NotImplementedError()
        for line in dtens:
            cols = line.strip().split()
            if len(cols) == 3:
                val = complex(float(cols[-2]), float(cols[-1]))
            else:
                val = float(cols[-1])
            crds = cols[0].split('-')
            ixyz = [xyz2id[xyz] for crd in crds for xyz in crd]
            if len(ixyz) == 3:
                tensor[ixyz[0]][ixyz[1]][ixyz[2]] = val
                if dosym:
                    if len(crds[0]) == 2:
                        tensor[ixyz[1]][ixyz[0]][ixyz[2]] = val
                    else:
                        tensor[ixyz[0]][ixyz[2]][ixyz[1]] = val
            else:
                tensor[ixyz[0]][ixyz[1]] = val
        data[incfrq][counts[incfrq]] = tensor

    dobj.set(data=data)
    return dobj


def __del_nonactive_modes(dobj_vtrans: QData,
                          dobj_vlevel: QData) -> tuple[QData, QData]:
    """Delete non-active modes.

    Selects and extracts only active mode from `vlevel`.

    Parameters
    ----------
    dobj_vtrans
        QData object with information on transition levels.
    dobj_vlevel
        QData object with Vibrational transition energies.

    Returns
    -------
    QData
        `dobj_vtrans` with only active modes kept.
    QData
        `dobj_vlevel` with only active modes kept.
    """
    vtrans = dobj_vtrans.data
    vlevel = dobj_vlevel.data
    nonact_modes = []
    for _, dtrans in vtrans.items():
        if dtrans[-1] is not None:
            if dtrans[-1] in ('inactive', 'passive'):
                to = dtrans[1]
                if len(to) > 1 or to[0][1] != 1:
                    msg = 'Status information should be only on fundamentals'
                    raise ValueError(msg)
                nonact_modes.append(to[0][0])
    if nonact_modes:
        key0 = 1
        newtrans = {key: val for key, val in vtrans.items() if key != 'qlabel'}
        for key, dtrans in newtrans.items():
            to = dtrans[1]
            remove = False
            if len(to) > 1 or to[0][1] != 1:
                for mode, _ in to:
                    if mode in nonact_modes:
                        remove = True
            if remove:
                del vlevel[key]
                del vtrans[key]
            else:
                vtrans[key0] = vtrans.pop(key)
                vlevel[key0] = vlevel.pop(key)
                key0 += 1
    new_vtrans = dobj_vtrans.copy(exclude=('data',))
    new_vlevel = dobj_vlevel.copy(exclude=('data',))
    new_vtrans.set(data=vtrans)
    new_vlevel.set(data=vlevel)

    return new_vtrans, new_vlevel


def get_hess_data(dfobj: tp.Optional[GLogIO] = None,
                  get_evec: bool = True,
                  get_eval: bool = True,
                  pre_data: tp.Optional[TypeQData] = None
                  ) -> tp.Tuple[tp.Any]:
    """Get or build Hessian data (eigenvectors and values).

    This function retrieves or builds the eigenvectors and eigenvalues.
    Contrary to ``get_data`` which only looks for available data, this
    functions looks for alternative forms to build necessary data.
    It also returns a Numpy array instead of Python lists.
    Preloaded data can be provided to avoid duplicating extraction
    queries.  They are expected in the same format as given by
    GLogIO methods.

    Parameters
    ----------
    dfobj
        Gaussian output file as `GLogIO` object.
    get_evec
        Return the eigenvectors.
    get_eval
        Return the eigenvalues.
    pre_data
        Database with quantities already loaded from previous queries.

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
    if pre_data is not None:
        for key in pre_data:
            qlab = pre_data[key].get('qlabel', key)
            if qlab.label == 'natoms':
                natoms = pre_data['natoms'].data
            elif qlab.label == 'nvib':
                natoms = pre_data['nvib'].data
            elif qlab.label == 'atmas':
                # np.repeat used to duplicate atmas for each Cart. coord.
                atmas = np.repeat(np.array(pre_data[key].data), 3)
            elif qlab.label == 'hessvec':
                key_evec = key
                hessvec = np.array(pre_data[key].data)
            elif qlab.label == 'hessval' and qlab.kind == 'freq':
                hessval = np.array(pre_data[key].data)
            elif qlab.label == 1 and qlab.derord == 2 and qlab.dercrd == 'X':
                key_ffx = key
                fccart = np.array(pre_data[key].data)

    if fccart is not None and atmas is not None:
        if (len(fccart.shape) == 1
                or pre_data[key_ffx].shape.lower() == 'lt'):
            fccart = square_ltmat(pre_data[key_ffx].data)
        eigvec, eigval = build_evec(fccart, atmas, get_evec, get_eval, nvib)
    else:
        if get_evec and hessvec is not None:
            # Check if eigenvectors need to be corrected
            #  Default is assumed to be Gaussian fchk
            evec_form = pre_data[key_evec].dtype == 'L.M^{-1/2}'
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
            read_data['nvib'] = QLabel('nvib')
        if calc_evec:
            read_data['hessvec'] = QLabel(quantity='hessvec', level='H')
        if calc_eval:
            read_data['hessval'] = QLabel(quantity='hessdat',
                                          descriptor='freq', level='H')
        tmp_data = get_data(dfobj, error_noqty=False, **read_data)

        if tmp_data['d2EdX2'] is not None and tmp_data['atmas'] is not None:
            # Rediagonlizing is more accurate, but require being
            atmas = np.repeat(np.array(tmp_data['atmas'].data), 3)
            ffx = square_ltmat(tmp_data['d2EdX2'].data)
            nvib = None if tmp_data['nvib'] is None else tmp_data['nvib'].data
            eigvec, eigval = build_evec(ffx, atmas, get_evec, get_eval, nvib)
        else:
            if calc_evec:
                if (tmp_data['hessvec'] is not None
                        and tmp_data['atmas'] is not None):
                    hessvec = tmp_data['hessvec'].data
                    atmas = np.repeat(np.array(tmp_data['atmas'].data), 3)
                    natoms = tmp_data['natoms'].data
                    evec_form = tmp_data['hessvec'].dtype == 'L.M^{-1/2}'
                    eigvec = convert_evec(hessvec, atmas, natoms, evec_form)
                else:
                    msg = 'Unable to retrieve force constants eigenvectors'
                    raise QuantityError(msg)
            if calc_eval:
                if tmp_data['hessval'] is not None:
                    eigval = np.array(tmp_data['hessval'].data)
                else:
                    msg = 'Unable to retrieve normal mode wavenumbers'
                    raise QuantityError(msg)

    return eigvec, eigval
