"""Low-level operations on Gaussian output (log) files.

Provides low-level interfaces to manipulate/extract data in Gaussian
output files.
"""

import os
import re
import typing as tp
from math import sqrt

from estampes.base import QLabel, \
    ParseKeyError, QuantityError, \
    TypeDGLog, TypeQData

from estampes.data.physics import phys_fact
from estampes.parser.functions import parse_qlabels
from estampes.parser.gaussian.glog.seeker import qlab_to_linkdata
from estampes.parser.gaussian.glog.extractor import parse_data
from estampes.parser.gaussian.glog.types import TypeKData


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
            links, keys, skips, fmts, ends, nums = qlab_to_linkdata(qlab)
            if isinstance(links, tuple):
                nblocks = len(links)
                key2blk[item] = (i, i+nblocks-1)
                for link, key, skip, fmt, end, num in zip(links, keys, skips,
                                                          fmts, ends, nums):
                    keydata.append((key, link, skip, re.compile(fmt), num,
                                    end))
                i += nblocks
            else:
                keydata.append((keys, links, skips, re.compile(fmts), nums,
                                ends))
                key2blk[item] = (i, i)
                i += 1
            qtydata[item] = qlab
        ndata, data = self.read_data(*keydata)
        dobjs = parse_data(qtydata, key2blk, ndata, data)
        self.__route = dobjs['route'].data
        # We merge all routes as one, assuming the blocks are coherent (user
        # did not run the same type of job linked)
        self.__links = sorted(set([int(line[0])
                                   for route in self.__route
                                   for line in route]))
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
                else:
                    for k, blk2id in enumerate(block2id):
                        a, b = blk2id
                        if a == i and b > j:
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


def get_hess_data(dfobj: tp.Optional[GLogIO] = None,
                  get_evec: bool = True,
                  get_eval: bool = True,
                  pre_data: tp.Optional[TypeQData] = None
                  ) -> tuple[tp.Any]:
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
    QuantityError
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
