"""Provide the necessary information to extract data.

This module provides mainly an entry function to `datakeys` functions,
which acts as a wrapper.
"""

import re
import typing as tp

from estampes.base import QLabel, \
    ParseKeyError
from estampes.parser.gaussian.glog.search_keys import keys_prp_3xx
from estampes.parser.gaussian.glog.types import TypeQKwrd


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
                if qlab.level == 'H':
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
                    key1.append('          ' +
                                'Anharmonic Raman Spectroscopy (Dynamic)')
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
                    elif isinstance(qlab.rstate, int):
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
                    elif isinstance(qlab.rstate, int):
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
            elif qlab.label in range(300, 400):
                lnk2, key2, sub2, fmt2, end2, num2 = keys_prp_3xx(qlab)
                lnk1.extend(lnk2)
                key1.extend(key2)
                sub1.extend(sub2)
                fmt1.extend(fmt2)
                end1.extend(end2)
                num1.extend(num2)
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
                        key1.append(f'{txt:>49s}')
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
            lnk = tuple(lnk1)
            key = tuple(key1)
            sub = tuple(sub1)
            end = tuple(end1)
            fmt = tuple(fmt1)
            num = tuple(num1)

    return lnk, key, sub, fmt, end, num
