"""Provide the necessary information to extract data.

This module provides mainly an entry function to `datakeys` functions,
which acts as a wrapper.
"""

import re
from collections.abc import Sequence

from estampes.base import QLabel, \
    InternalError, ParseKeyError, QuantityError
from estampes.parser.gaussian.glog.search_keys import keys_arch, keys_prp_3xx
from estampes.parser.gaussian.glog.types import QKwrdType
from estampes.parser.gaussian.glog.logkeys import RR_OMEGA_LINE, \
    RR_OMEGA_UNIT, RR_OMEGA_VAL, KEY_DP, KEY_FP, KEY_UINT


def qlab_to_linkdata(qlab: QLabel,
                     gver: Sequence[str] | None = None) -> QKwrdType:
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
        Gaussian version, as (major, minor).

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
    lnk1 = []
    key1 = []
    sub1 = []
    end1 = []
    fmt1 = []
    num1 = []
    if qlab.label == 'route':  # Log specific extraction, not provided in base
        lnk1.extend((1, 1))
        # key = (' Cite this work as:', ' Link1:  Proceeding to internal job')
        key1.extend((' Cite this work as:', ' Normal termination of Gaussian'))
        sub1.extend((' -', ' ---'))
        end1.extend((
            lambda s: s.startswith(' Charge ='),
            lambda s: s.startswith(' Charge =')))
        fmt1.extend((
            r'^ (?P<val>\d+\/(?:,?\d+=-?\d+)*\/(?:,?\d+)+(?:\(-\d\))?);\s*$',
            r'^ (?P<val>\d+\/(?:,?\d+=-?\d+)*\/(?:,?\d+)+(?:\(-\d\))?);\s*$'))
        num1.extend((0, 1))
    elif qlab.label == 'natoms':
        lnk1.append(101)
        key1.append(' NAtoms=')
        sub1.append(0)
        end1.append(lambda _s: True)
        fmt1.append(r'^ NAtoms=\s+(?P<val>\d+)\s+N\w+=\s+.*$')
        num1.append(0)
    elif qlab.label == 'nvib':
        # We load anyway the frequencies to count them.
        # This is a bit of an absurd way to proceed but only way to be
        #   sure that we have the correct number since Gaussian does not
        #   list it explicitly in the output
        lnk1.append(716)
        key1.append(' and normal coordinates:')
        sub1.append(1)
        end1.append(lambda s: s.startswith(' - Thermochemistry'))
        fmt1.append(r'^\s+Frequencies -- \s*(?P<val>\d.*)\s*$')
        num1.append(-1)
    elif qlab.label == 'atmas':
        lnk1.extend((-101, 716))
        key1.extend((' NAtoms= ', ' - Thermochemistry -'))
        sub1.extend((0, 3))
        end1.extend((lambda s: s.startswith(' Leave Link'),
                     lambda s: s.startswith(' Molecular Mass:')))
        fmt1.extend((r'^ AtmWgt=\s+(?P<val>(\s+\d+\.\d+)+)\s*$',
                     r'^ Atom\s+\d+ has atomic number\s+\d+ and '
                     + r'mass\s+(?P<val>\d+\.\d+)\s*$'))
        num1.extend((0, 0))
    elif qlab.label == 'atnum':
        lnk1.extend((0, 0))
        key1.extend(('                         Standard orientation:',
                     '                          Input orientation:'))
        sub1.extend((5, 5))
        end1.extend((lambda s: s.startswith(' ------'),
                     lambda s: s.startswith(' ------')))
        txt = r'^\s+\d+\s+(?P<val>\d+)\s+\d+(?:\s+-?\d+\.\d+){3}\s*$'
        fmt1.extend((txt, txt))
        num1.extend((0, 0))
    elif qlab.label == 'molsym':
        lnk1.append(101)
        key1.append(' Framework group ')
        sub1.append(0)
        end1.append(lambda _s: True)
        fmt1.append(r'^ Framework group \s+(?P<val>[^ ]+)\s*$')
        num1.append(0)
    elif qlab.label == 'atcrd' or qlab.label == 2:
        lnk1.extend((0, 0))
        key1.extend(('                         Standard orientation:',
                     '                          Input orientation:'))
        sub1.extend((5, 5))
        end1.extend((lambda s: s.startswith(' ------'),
                     lambda s: s.startswith(' ------')))
        txt = r'^(?:\s+\d+){3}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'
        fmt1.extend((txt, txt))
        if qlab.kind == 'first':
            num1.extend((0, 0))
        elif qlab.kind in ('last', 'orient'):
            num1.extend((-1, -1))
        else:
            num1.extend((1, 1))
    elif qlab.label == 'hessvec':
        lnk1.extend((-716, 716))
        key1.extend((' and normal coordinates:',
                     ' and normal coordinates:'))
        sub1.extend((1, 1))
        end1.extend((lambda s: s.startswith(' - Thermochemistry'),
                     lambda s: s.startswith(' - Thermochemistry')))
        fmt1.extend((r'^\s+(?P<val>(?:\d+\s+){3}(?:\s+-?\d\.\d+){1,5})\s*$',
                     r'^\s+(?P<val>(?:\d+\s+){2}(?:\s+-?\d\.\d+){1,9})\s*$'))
        num1.extend((0, -1))
    elif qlab.label == 'hessdat':
        lnk1.extend((-716, 716))
        key1.extend((' and normal coordinates:',
                     ' and normal coordinates:'))
        sub1.extend((1, 1))
        end1.extend((lambda s: s.startswith(' - Thermochemistry'),
                     lambda s: s.startswith(' - Thermochemistry')))
        if qlab.kind == 'freq':
            fmt1.extend((r'^\s+Frequencies --- \s*(?P<val>-?\d.*)\s*$',
                         r'^\s+Frequencies -- \s*(?P<val>-?\d.*)\s*$'))
        elif qlab.kind == 'redmas':
            fmt1.extend((r'^\s+Reduced masses --- \s*(?P<val>\d.*)\s*$',
                         r'^\s+Red. masses -- \s*(?P<val>\d.*)\s*$'))
        else:
            raise NotImplementedError('Unknown subopt for HessDat')
        num1.extend((0, -1))
    elif qlab.label == 'swopt':
        lnk1.append(1)
        key1.append(' Cite this work as:')
        sub1.append(' #')
        end1.append(lambda s: s.startswith(' -'))
        fmt1.append(r'^ (?P<val>.*)$')
        num1.append(0)
    elif qlab.label == 'swver':
        lnk1.append(1)
        key1.append(' ****')
        sub1.append(1)
        end1.append(lambda _s: True)
        fmt1.append(r'^ (?P<val>Gaussian (?:\d\d|DV):\s.*)$')
        num1.append(0)
    elif qlab.label == 'intens':
        if qlab.kind == 'IR':
            if qlab.level == 'H':
                lnk1.extend((-716, 716, -717))
                key1.extend((' and normal coordinates:',
                             ' and normal coordinates:',
                             '        Integrated intensity (I)'))
                sub1.extend((1, 1, 4))
                end1.extend((lambda s: s.startswith(' - Thermochemistry'),
                             lambda s: s.startswith(' - Thermochemistry'),
                             lambda s: s.startswith(' -----')))
                fmt1.extend((r'^\s+IR Intensities --- \s*(?P<val>\d.*)\s*$',
                             r'^\s+IR Inten    -- \s*(?P<val>\d.*)\s*$',
                             r'^\s+\d+\(\d+\)\s+(?:-?\d+\.\d+\s+|\*+\s+){2}'
                             + r'(?P<val>-?\d+\.\d+|\*+)\s+'
                             + r'(?:-?\d+\.\d+|\*+)\s*$'))
                num1.extend((0, -1, 0))
            elif qlab.level == 'A':
                lnk1.append(717)
                key1.append('        Integrated intensity (I)')
                sub1.append(3)
                end1.append(lambda s: s.startswith(' Units:'))
                fmt1.append(r'^\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'
                            + r' .*\s+(?P<val>-?\d+\.\d+|\*+)\s*$')
                num1.append(0)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()
    elif qlab.label == 'fcdat':
        if qlab.kind == 'SimInf':
            lnk1.append(-718)
            key1.append('               Information on the Simulation')
            sub1.append(2)
            end1.append(lambda s: s.startswith('     ===='))
            fmt1.append(r'^\s+(?P<val>.*\w.*)\s*$')
            num1.append(0)
        elif qlab.kind == 'ExcState':
            lnk1.append(-718)
            key1.append('                  Treatment of Input Data')
            sub1.append(2)
            end1.append(lambda s: s.startswith('     ===='))
            fmt1.append(r'^\s+(?P<val>No electronic transition.*$|'
                        + r'NOTE: Using excited electronic state number.*)$')
            num1.append(0)
        elif qlab.kind == 'JMat':
            lnk1.extend((-718, -718))
            key1.extend((' Duschinsky matrix', ' Final Duschinsky matrix'))
            sub1.extend((2, 2))
            end1.extend((lambda s: not s.strip(),
                         lambda s: not s.strip()))
            txt = r'^\s+(?P<val>(?:\d+)' \
                + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$'
            fmt1.extend((txt, txt))
            num1.extend((-1, -1))
        elif qlab.kind == 'JMatF':
            lnk1.append(-718)
            key1.append(' Full Duschinsky matrix')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>(?:\d+)'
                        + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$')
            num1.append(0)
        elif qlab.kind == 'KVec':
            lnk1.append(-718)
            key1.append(' Shift Vector')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$')
            num1.append(0)
        elif qlab.kind == 'SRAMat':
            lnk1.append(-718)
            key1.append(' A Matrix')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>(?:\d+)'
                        + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$')
            num1.append(0)
        elif qlab.kind == 'SRBVec':
            lnk1.append(-718)
            key1.append(' B Vector')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$')
            num1.append(0)
        elif qlab.kind == 'SRCMat':
            lnk1.append(-718)
            key1.append(' C Matrix')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>(?:\d+)'
                        + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$')
            num1.append(0)
        elif qlab.kind == 'SRDVec':
            lnk1.append(-718)
            key1.append(' D Vector')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?:\d+)\s+(?P<val>-?\d\.\d+D?[\+-]\d{2,3})\s*$')
            num1.append(0)
        elif qlab.kind == 'SREMat':
            lnk1.append(-718)
            key1.append(' E Matrix')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>(?:\d+)'
                        + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$')
            num1.append(0)
        elif qlab.kind == 'Spec':
            lnk1.extend((-718, -718))
            key1.extend(('                       Final Spectrum',
                         ' Legend:'))
            sub1.extend((' -----------', ' -----------'))
            end1.extend((lambda s: not s.strip(),
                         lambda s: 'Legend:' in s or not s.strip()))
            fmt1.extend((r'^\s+(?P<val>(?:-?\d+.\d+)'
                         + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3})+)\s*$',
                         r'^\s+(?P<val>(?:-?\d+.\d+)'
                         + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3})+)\s*$'))
            num1.extend((0, 1))
        elif qlab.kind == 'SpcPar':
            lnk1.extend((718, 718))
            key1.extend(('                       Final Spectrum',
                        ' Legend:'))
            sub1.extend((3, 0))
            end1.extend((lambda s: s.startswith(' -----------'),
                         lambda s: s.startswith(' -----------')))
            # In format, \w is used to exclude empty lines
            fmt1.extend((r'^\s+(?P<val>.*\w.*)\s*$',
                         r'^\s+(?P<val>.*\w.*)\s*$'))
            num1.extend((0, 1))
        elif qlab.kind == 'Conv':
            # 2 sets, one for the main one, one to extract the FCF (HT)
            lnk1.extend((-718, -718))
            key1.extend(('              Calculations of Band Intensities',
                         '              Calculations of Band Intensities'))
            sub1.extend((2, 2))
            end1.extend((lambda s: s.startswith('     ===='),
                         lambda s: s.startswith('     ====')))
            fmt1.extend((
                r'^\s+Spectrum progression:\s+(?P<val>-?\d+\.\d+)%\s*$',
                r'^\s+.+Franck-Condon Factors:\s+(?P<val>-?\d+\.\d+)%)\s*$'))
            num1.extend((0, 0))
        elif qlab.kind == 'Assign':
            lnk1.extend((-718, -718))
            key1.extend(('                 Information on Transitions',
                         '                 Information on Transitions'))
            sub1.extend((2, 2))
            end1.extend((lambda s: s.startswith('     ===='),
                         lambda s: s.startswith('     ====')))
            fmt1.extend((r'^\s+Energy =\s+(?P<val>-?\d+\.\d+ cm.-1: .*)\s*$',
                         r'^\s+-. Intensity =\s+(?P<val>.*)\s*$'))
            num1.extend((0, 0))
        elif qlab.kind == 'RedDim':
            lnk1.append(-718)
            key1.append(' Reduced system')
            sub1.append(2)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>\d+\s+=\s+\d+\s+\d+\s+=\s+\d+)\s*$')
            num1.append(0)
        elif qlab.kind == 'E(0-0)':
            lnk1.append(718)
            key1.append('                 Information on Transitions')
            sub1.append(2)
            end1.append(lambda s: s.startswith('     ===='))
            fmt1.append(r'^\s+Energy of the 0-0 transition:\s+'
                        + r'(?P<val>-?\d+\.\d+ cm\S+)\s*$')
            num1.append(0)
        elif qlab.kind == 'GeomIS':
            lnk1.append(718)
            key1.append('              New orientation in initial state')
            sub1.append(5)
            end1.append(lambda s: s.startswith(' ------'))
            fmt1.append(r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$')
            num1.append(0)
        elif qlab.kind == 'GeomFS':
            # The second block is for post G16 versions
            lnk1.extend((718, 718))
            key1.extend(('              New orientation in final state',
                         '               New orientation in final state'))
            sub1.extend((5, 5))
            end1.extend((lambda s: s.startswith(' ------'),
                         lambda s: s.startswith(' ------')))
            fmt1.extend((
                r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$',
                r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'))
            num1.extend((0, 0))
        elif qlab.kind == 'GeomMS':
            # The second block is for post G16 versions
            lnk1.extend((718, 718))
            key1.extend(('              New orientation in intermediate state',
                         '           New orientation in intermediate state'))
            sub1.extend((5, 5))
            end1.extend((lambda s: s.startswith(' ------'),
                         lambda s: s.startswith(' ------')))
            fmt1.extend((
                r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$',
                r'^(?:\s+\d+){2}(?P<val>(?:\s+-?\d+\.\d+){3})\s*$'))
            num1.extend((0, 0))
        elif qlab.kind == 'ExGeom':
            lnk1.append(-718)
            key1.append(' Extrapolated geometry')
            sub1.append(4)
            end1.append(lambda s: s.startswith(' ------'))
            fmt1.append(r'^\s+\w+(?P<val>(?:\s+-?\d+\.\d+){3})\s*$')
            num1.append(0)
    elif qlab.label == 'vptdat':
        if qlab.kind == 'CICoef':
            # The second group is to correct the numbering if passive modes
            #     present.
            lnk1.extend((717, 717))
            key1.extend((' Definition of New States w.r.t. Deperturbed States',
                         ' Reduced-Dimensionality on Variational States'))
            sub1.extend((4, 3))
            end1.extend((lambda s: not s.strip(),
                         lambda s: not s.strip()))
            fmt1.extend((r'^\s*(?P<val>(\d+\s+:|)\s+[+-]?\d\.\d+\s+x\s+'
                         + r'\|[0-9();]+>)\s*$',
                         r'^\s*(?P<val>\d+\s+\|\s+\d+)\s*$'))
            num1.extend((0, 0))
        elif qlab.kind == 'NMOrder':
            lnk1.append(717)
            key1.append('        Vibro-Rotational Analysis Based on Symmetry')
            sub1.append(2)
            end1.append(lambda s: '=====' in s)
            fmt1.append(r'^\s*(?P<val>\([HA]\)\s+\|(?:\s+\d+\|)+)\s*$')
            num1.append(0)
        elif qlab.kind == 'NMFlags':
            lnk1.append(717)
            key1.append(' Definition of the model system: Active modes')
            sub1.append(1)
            end1.append(lambda s: not s.strip())
            fmt1.append(
                r'^\s*(?P<val>'
                + r'(?:The \d+ (?:Active|Inactive|Passive) Modes are:|'
                + r'(?:\s*\d+)+))\s*$')
            num1.append(0)
        elif qlab.kind == 'XMat':
            # The second group is to correct the numbering if passive modes
            #     present.
            lnk1.append(717)
            key1.append(' Total Anharmonic X Matrix (in cm^-1)')
            sub1.append(1)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>\d+(?:\s+\d+|'
                        + r'\s+-?\d+\.\d+D?[-+]\d+){1,5})\s*$')
            num1.append(0)
        elif qlab.kind == 'YMat':
            # The code is very similar to XMat except the matrix is not
            # necessarily symmetric anymore.
            # Hence, we also need to check for passive modes present.
            lnk1.append(717)
            key1.append('                    Anharmonic Y Matrix')
            sub1.append(1)
            end1.append(lambda s: not s.strip())
            fmt1.append(r'^\s+(?P<val>(?:\d+)'
                        + r'(?:\s+-?\d\.\d+D?[\+-]\d{2,3}){1,5})\s*$')
            num1.append(0)
        elif qlab.kind == 'freq':
            lnk1.append(717)
            key1.append(' :      QUADRATIC FORCE CONSTANTS IN NORMAL MODES')
            sub1.append(6)
            end1.append(lambda s: '...........' in s)
            fmt1.append(r'^\s+(?P<val>(?:' + KEY_UINT + r'){2}\s+'
                        + KEY_FP + r')(?:\s+' + KEY_FP + r'){2}\s*$')
            num1.append(0)
        elif qlab.kind == 'cubic':
            lnk1.append(717)
            key1.append(' :        CUBIC FORCE CONSTANTS IN NORMAL MODES')
            sub1.append(6)
            end1.append(lambda s: '...........' in s)
            fmt1.append(r'^\s+(?P<val>(?:' + KEY_UINT + r'){3}\s+'
                        + KEY_FP + r')(?:\s+' + KEY_FP + r'){2}\s*$')
            num1.append(0)
        elif qlab.kind == 'quartic':
            lnk1.append(717)
            key1.append(' :       QUARTIC FORCE CONSTANTS IN NORMAL MODES')
            sub1.append(6)
            end1.append(lambda s: '=====' in s)
            fmt1.append(r'^\s+(?P<val>(?:' + KEY_UINT + r'){4}\s+'
                        + KEY_FP + r')(?:\s+' + KEY_FP + r'){2}\s*$')
            num1.append(0)
        else:
            raise NotImplementedError()
    elif qlab.label == 'vtrans':
        if qlab.kind == 'RR':
            lnk1.extend((718, 718))
            key1.extend(('                 Information on Transitions',
                         '                 Information on Transitions'))
            sub1.extend((2, 2))
            end1.extend((lambda s: s.startswith('     ====='),
                         lambda s: s.startswith('     =====')))
            fmt1.extend((r'^\s+Energy = \s*-?\d+\.\d+ cm.-1:\s+'
                         + r'(?P<val>\|.+ ->\s+\|.+)\s*$',
                         RR_OMEGA_LINE))
            num1.extend((0, 0))
        elif qlab.kind == 'SOS':
            raise NotImplementedError()
        else:
            if qlab.level == 'H':
                lnk1.extend((-716, 716, -717))
                key1.extend((' and normal coordinates:',
                             ' and normal coordinates:',
                             ' NOTE: Transition energies are given with'))
                sub1.extend((1, 1, 8))
                end1.extend((lambda s: s.startswith(' - Thermochemistry'),
                             lambda s: s.startswith(' - Thermochemistry'),
                             lambda s: s.startswith('     =====')))
                fmt1.extend((r'^\s{16}(?P<val>(?:\s+\d+){1,5})\s*$',
                             r'^\s{16}(?P<val>(?:\s+\d+){1,3})\s*$',
                             r'^\s+\w?\s+(?P<val>\s*\d+\(\d+\))\s+\w+\s+'
                             + r'(?:\s+-?\d+\.\d+|\*+){4}.*\s*$'))
                num1.extend((0, -1, 0))
            elif qlab.level == 'A':
                lnk1.append(717)
                key1.append(' NOTE: Transition energies are given with')
                sub1.append(8)
                end1.append(lambda s: s.startswith('     ====='))
                # fmt = r'^\s+\w?\s+(?P<val>(?:\s*\d+\(\d+\)){1,3}|\d+)\s+'\
                #       + r'(?:\w+)?\s+(?:\s*-?\d+\.\d+|\*+\s+){4,5}.*\s*$'
                fmt1.append(
                    r'^\s+\w?\s+(?P<val>(?:\s*\d+\(\d+\)\s+(?:\w+)?)'
                    + r'{1,3}|\d+)\s+(?:\s*-?\d+\.\d+|\*+\s+){4,5}.*\s*$')
                num1.append(0)
            else:
                raise NotImplementedError()
    elif qlab.label == 'vlevel':
        if qlab.kind == 'RR':
            lnk1.extend((718, 718))
            key1.extend(('                 Information on Transitions',
                         '                 Information on Transitions'))
            sub1.extend((2, 2))
            end1.extend((lambda s: s.startswith('     ====='),
                         lambda s: s.startswith('     =====')))
            fmt1.extend((r'^\s+Energy = \s*(?P<val>-?\d+\.\d+) cm.-1:\s+\|.+$',
                         RR_OMEGA_LINE))
            num1.extend((0, 0))
        elif qlab.kind == 'SOS':
            raise NotImplementedError()
        else:
            if qlab.level == 'H':
                lnk1.extend((-716, 716, -717))
                key1.extend((' and normal coordinates:',
                             ' and normal coordinates:',
                             ' NOTE: Transition energies are given with'))
                sub1.extend((1, 1, 8))
                end1.extend((lambda s: s.startswith(' - Thermochemistry'),
                             lambda s: s.startswith(' - Thermochemistry'),
                             lambda s: s.startswith('     =====')))
                fmt1.extend((
                    r'^\s+Frequencies --- \s*(?P<val>-?\d.*)\s*$',
                    r'^\s+Frequencies -- \s*(?P<val>-?\d.*)\s*$',
                    r'^\s+\w?\s+(?:\s*\d+\(\d+\))\s+\w+\s+'
                    + r'(?P<val>-?\d+\.\d+|\*+)(?:\s+-?\d+\.\d+|\*+){4}.*\s*$'
                    ))
                num1.extend((0, -1, 0))
            elif qlab.level == 'A':
                lnk1.append(717)
                key1.append(' NOTE: Transition energies are given with')
                sub1.append(8)
                end1.append(lambda s: s.startswith('     ====='))
                fmt1.append(
                    r'^\s+\w?\s+(?:(?:\s*\d+\(\d+\)){1,3}|\d+)\s+(?:\w+)?'
                    + r'\s+(?:-?\d+\.\d+|\*+)?\s+(?P<val>-?\d+\.\d+|\*+)'
                    + r'(?:\s*-?\d+\.\d+|\*+){3}.*\s*$')
                num1.append(0)
            else:
                raise NotImplementedError()
    else:
        lnk0 = -914
        key0 = ' Excitation energies and oscillator strengths:'
        sub0 = 2

        def end0(s):
            return s.startswith(' End of Minotr F.D. properties file')
        fmt0 = r'^\s+(?P<val>Excited State\s+\d+: .*|' +\
            r'This state for optimization.*)\s*$'
        num0 = 0
        if isinstance(qlab.rstate, tuple):
            state_i, state_f = qlab.rstate
            if qlab.label == 1:
                if qlab.derord == 0:
                    if state_i == 0:
                        lnk1.append(0)
                        key1.append(
                            ' Excitation energies and oscillator strengths:')
                        sub1.append(1)
                        end1.append(lambda s: s.startswith(' ****'))
                        if isinstance(state_f, int):
                            txt = str(state_f)
                        elif state_f == 'a':
                            txt = r'\d+'
                        else:
                            raise ValueError('Unsupported final state')
                        fmt1.append(
                            r'^ Excited State\s+' + txt
                            + r':\s+\S+\s+(?P<val>\d+\.\d+)\b\s+eV.*$')
                        num1.append(-1)
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()
            elif qlab.label in (101, 102, 107):
                if qlab.derord == 0:
                    if qlab.label == 101:
                        if qlab.kind == 'len':
                            key_prp = 'electric dipole'
                        elif qlab.kind == 'vel':
                            key_prp = 'velocity dipole'
                        else:
                            raise NotImplementedError(
                                'Unsupported dipolar formalism')
                        n_xyz = 3
                        n_extra = 2
                    elif qlab.label == 102:
                        key_prp = 'magnetic dipole'
                        n_xyz = 3
                        n_extra = 0
                    elif qlab.label == 107:
                        key_prp = 'velocity quadrupole'
                        n_xyz = 6
                        n_extra = 0
                    else:
                        raise NotImplementedError(
                            'Unsupported property for electronic transition '
                            'moments')
                    # To avoid duplicating the code, we construct a general
                    # system where we build sequentially the string to
                    # search.  The property-dependent part needs to account
                    # for different cases:
                    # for electric dipole: x, y, z, dip. str., osc. str.
                    # for magnetic dipole: x, y, z
                    # for quadrupole: xx, yy, zz, xy, xz, yz
                    # Note: final ) in first block to close the (?:P<val>
                    key_xyz = rf'(?:{KEY_FP}){{{n_xyz}}})'
                    if n_extra > 0:
                        key_xyz += rf'(?:{KEY_FP}){{{n_extra}}}\s*$'
                    else:
                        key_xyz += r'\s*$'
                    lnk1.append(0)
                    num1.append(-1)
                    if state_i == 0:
                        key1.append(' Ground to excited state transition '
                                    + f'{key_prp} moments (Au):')
                        if isinstance(state_f, int):
                            sub1.append(state_f + 1)
                            end1.append(lambda _s: True)
                        elif state_f == 'a':
                            sub1.append(1)
                            end1.append(lambda s: s[2] != ' ')
                        else:
                            raise ValueError('Unsupported final state')
                        fmt1.append(rf'^\s+\d+(?P<val>{key_xyz}')
                    elif isinstance(state_i, int):
                        key1.append(' Excited to excited state transition '
                                    + f'{key_prp} moments (Au):')
                        sub1.append(1)
                        end1.append(lambda s: s[2] != ' ')
                        if isinstance(state_f, int):
                            key_state = rf'^(?P<val>\s+{state_f}\s+{state_i}'
                        elif state_f == 'a':
                            key_state = rf'^(?P<val>\s+\d+\s+{state_i}'
                        else:
                            raise ValueError('Unsupported final state')
                        fmt1.append(key_state + key_xyz)
                    elif state_i == 'a':
                        key1.append(' Excited to excited state transition '
                                    + f'{key_prp} moments (Au):')
                        sub1.append(1)
                        end1.append(lambda s: s[2] != ' ')
                        if isinstance(state_f, int):
                            key_state = rf'^(?P<val>\s+{state_f}\s+\d+'
                        elif state_f == 'a':
                            key_state = r'^(?P<val>\s+\d\s+\d+'
                        else:
                            raise ValueError('Unsupported final state')
                        fmt1.append(key_state + key_xyz)
                    else:
                        raise NotImplementedError('Unsupported initial state')
                elif qlab.derord == 1:
                    # Check if transition specification is compatible
                    if state_i == 'a' or state_i != 0:
                        msg = 'Parsing of electronic transition moments ' \
                            'from excited states not yet supported for ' \
                            'properties'
                        raise NotImplementedError(msg)
                    if isinstance(state_f, int):
                        # We need to check if final state is the right one,
                        # but this can only be done a posteriori with processed
                        # data.
                        lnk1.append(lnk0)
                        key1.append(key0)
                        sub1.append(sub0)
                        end1.append(end0)
                        fmt1.append(fmt0)
                        num1.append(num0)
                    # The block seems only available with #P and is more a
                    # matrix dump, containing all the block of derivatives,
                    # with: energy, edip (len), edip (vel), mdip, equad, in
                    # rows
                    # The first 3 columns are electric-field derivatives.
                    # There is no specific end to the block, so we simply look
                    # if there is a character in column 2, since the number
                    # of properties should be limited to 16 or a relatively
                    # low number anyway.
                    # To find the right property, we simply add a regex on
                    # the row index
                    if qlab.label == 101:
                        if qlab.kind == 'len':
                            fmt_index = '(?: 2 | 3 | 4 )'
                        elif qlab.kind == 'vel':
                            fmt_index = '(?: 5 | 6 | 7 )'
                        else:
                            raise NotImplementedError(
                                'Unsupported dipolar formalism')
                    elif qlab.label == 102:
                        fmt_index = '(?: 8 | 9 | 10 )'
                    elif qlab.label == 107:
                        fmt_index = '(?: 11 | 12 | 13 | 14 | 15 | 16 )'
                    else:
                        raise NotImplementedError(
                            'Unsupported property for electronic transition '
                            'moments')
                    lnk1.append(9999)
                    key1.append(' Electronic Transition Derivatives')
                    sub1.append(1)
                    num1.append(0)
                    end1.append(lambda s: s[1] != ' ')
                    fmt1.append(
                        rf'^\s+{fmt_index}(?P<val>(?:{KEY_DP}){{1,5}})\s*$')
                else:
                    raise NotImplementedError(
                        'Electronic transition moments derivatives NYI')
            elif qlab.label == 'dipstr':
                if state_i == 0:
                    lnk1.append(0)
                    key1.append(
                        ' Ground to excited state transition electric '
                        + 'dipole moments (Au):')
                    if isinstance(state_f, int):
                        sub1.append(state_f + 1)
                        end1.append(lambda _s: True)
                    elif state_f == 'a':
                        sub1.append(1)
                        end1.append(lambda s: s[2] != ' ')
                    else:
                        raise ValueError('Unsupported final state')
                    fmt1.append(r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+'
                                + r'(?P<val>-?\d+\.\d+)\s+-?\d+\.\d+\s*$')
                    num1.append(-1)
                else:
                    raise NotImplementedError()
            elif qlab.label == 'rotstr':
                if state_i == 0:
                    lnk1.append(0)
                    key1.append(
                        ' Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm'
                        + '/Gauss)')
                    if isinstance(state_f, int):
                        sub1.append(state_f+1)
                        end1.append(lambda _s: True)
                    elif state_f == 'a':
                        sub1.append(1)
                        end1.append(lambda s: s[2] != ' ')
                    else:
                        raise ValueError('Unsupported final state')
                    if qlab.kind == 'vel':
                        num1.append(0)
                        fmt1.append(r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+'
                                    + r'(?P<val>-?\d+\.\d+)\s+-?\d+\.\d+\s*$')
                    else:
                        num1.append(-1)
                        fmt1.append(r'^\s+\d+(?:\s+-?\d+\.\d+){3}\s+'
                                    + r'(?P<val>-?\d+\.\d+)\s*$')
                else:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()
        else:
            if qlab.rstate == 'c':
                lnk1.append(lnk0)
                key1.append(key0)
                sub1.append(sub0)
                end1.append(end0)
                fmt1.append(fmt0)
                num1.append(num0)
            elif not isinstance(qlab.rstate, int):
                raise NotImplementedError()
            # Add quantity specific subgroups
            if qlab.label in ('ramact', 'roaact') and qlab.kind != 'static':
                # Information on setup and incident frequency
                if qlab.kind == 'RR':
                    lnk1.append(718)
                    key1.append('                 Information on Transitions')
                    sub1.append(2)
                    end1.append(lambda s: s.startswith('     ====='))
                    num1.append(0)
                    fmt1.append(RR_OMEGA_LINE)
                else:
                    # -- Harmonic level
                    if qlab.level == 'H':
                        lnk1.append(716)
                        key1.append(' and normal coordinates:')
                        sub1.append(1)
                        end1.append(
                            lambda s: not s.startswith(' Incident light:'))
                        num1.append(0)
                        fmt1.append(r'^\s+Incident light \(\S+\):'
                                    + r'\s+(?P<val>-?\d.*)\s*$')
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
                        arckey = keys_arch()
                        lnk1.extend((-716, arckey[0]))
                        key1.extend((
                            ' Force constants in Cartesian coordinates:',
                            arckey[1]))
                        sub1.extend((1, arckey[2]))
                        end1.extend((
                            lambda s:
                                not re.match(r'^\s+\d+(\s+-?\d|\s*$)', s),
                            arckey[3]))
                        fmt1.extend((
                            r'^\s+(?P<val>(?:\s*\d+|'
                            + r'\s*\d+\s+-?\d+\.\d+D?[-+]\d+){1,5})\s*$',
                            arckey[4]))
                        num1.extend((0, arckey[5]))
                    elif isinstance(qlab.rstate, int):
                        raise NotImplementedError()
                elif qlab.derord == 3:
                    if qlab.rstate == 'c':
                        lnk1.append(-717)
                        key1.append(
                            ' :        CUBIC FORCE CONSTANTS IN NORMAL MODES')
                        sub1.append(6)
                        end1.append(lambda s: '...........' in s)
                        if qlab.dercrd == 'QRED':
                            fmt1.append(
                                r'^\s+(?P<val>(?:' + KEY_UINT + r'){3}\s+'
                                + KEY_FP + r')(?:\s+' + KEY_FP + r'){2}\s*$')
                        elif qlab.dercrd == 'Q':
                            fmt1.append(
                                r'^\s+(?P<val>(?:' + KEY_UINT + r'){3}(?:\s+'
                                + KEY_FP + r'){3})\s*$')
                        else:
                            msg = 'Cubic force constants only supported wrt ' \
                                + 'normal coordinates'
                            raise NotImplementedError(msg)
                        num1.append(0)
                    elif isinstance(qlab.rstate, int):
                        raise NotImplementedError()
                elif qlab.derord == 4:
                    if qlab.rstate == 'c':
                        lnk1.append(-717)
                        key1.append(
                            ' :       QUARTIC FORCE CONSTANTS IN NORMAL MODES')
                        sub1.append(6)
                        end1.append(lambda s: '=====' in s)
                        if qlab.dercrd == 'QRED':
                            fmt1.append(
                                r'^\s+(?P<val>(?:' + KEY_UINT + r'){4}\s+'
                                + KEY_FP + r')(?:\s+' + KEY_FP + r'){2}\s*$')
                        elif qlab.dercrd == 'Q':
                            fmt1.append(
                                r'^\s+(?P<val>(?:' + KEY_UINT + r'){4}(?:\s+'
                                + KEY_FP + r'){3})\s*$')
                        else:
                            msg = 'Quartic force constants only supported ' \
                                + 'wrt normal coordinates'
                            raise NotImplementedError(msg)
                        num1.append(0)
                    elif isinstance(qlab.rstate, int):
                        raise NotImplementedError()
                else:
                    raise NotImplementedError(
                        'Unsupported derivative order of energy')
            elif qlab.label == 101:
                if qlab.derord == 0:
                    if qlab.rstate == 'c':
                        arckey = keys_arch()
                        lnk1.extend((601, arckey[0]))
                        key1.extend((' Dipole moment (field-independent basis,'
                                     + ' Debye):',
                                     arckey[1]))
                        sub1.extend((1, arckey[2]))
                        end1.extend((lambda s: True, arckey[3]))
                        fmt1.extend((r'^\s+(?P<val>(?:\s*[XYZ]=\s+-?\d+\.\d+)'
                                    + r'{3}).*$', arckey[4]))
                        num1.extend((0, arckey[5]))
                elif qlab.derord == 1:
                    if qlab.rstate == 'c':
                        arckey = keys_arch()
                        lnk1.append(arckey[0])
                        key1.append(arckey[1])
                        sub1.append(arckey[2])
                        end1.append(arckey[3])
                        fmt1.append(arckey[4])
                        num1.append(arckey[5])
                else:
                    raise NotImplementedError(
                        'Unsupported derivative order of electric dipole')
            elif qlab.label == 102:
                if qlab.derord == 0:
                    if qlab.rstate == 'c':
                        raise QuantityError('Magnetic dipole not available.')
                elif qlab.derord == 1:
                    if qlab.rstate == 'c':
                        arckey = keys_arch()
                        lnk1.append(arckey[0])
                        key1.append(arckey[1])
                        sub1.append(arckey[2])
                        end1.append(arckey[3])
                        fmt1.append(arckey[4])
                        num1.append(arckey[5])
                else:
                    raise NotImplementedError(
                        'Unsupported derivative order of magnetic dipole')
            elif qlab.label in range(300, 400):
                lnk2, key2, sub2, fmt2, end2, num2 = keys_prp_3xx(qlab, gver)
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
                        fmt1.append(RR_OMEGA_UNIT)
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
                                RR_OMEGA_VAL
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
                                RR_OMEGA_VAL
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
                                RR_OMEGA_VAL
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
                                RR_OMEGA_VAL
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
                                RR_OMEGA_VAL
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
                if qlab.kind == 'RR':
                    lnk1.append(718)
                    key1.append('                 Information on Transitions')
                    sub1.append(2)
                    end1.append(lambda s: s.startswith('     ====='))
                    num1.append(0)
                    fmt1.append(
                        r'^\s+-> (?P<val>Omega =\s*\d+\.\d+ cm.-1.* '
                        + r'Sigma =\s*-?\d+\.\d+E?[+-]\d{2,3})\s*$')
                else:
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
                                [r'^\s+Raman Activities ---\s+(?P<val>-?\d.*)'
                                 + r'\s*$',
                                 r'^\s+Raman Activ -- \s*(?P<val>-?\d.*)\s*$'])
                        else:
                            lnk1.append(716)
                            key1.append(' and normal coordinates:')
                            end1.append(
                                lambda s: s.startswith(' - Thermochemistry'))
                            sub1.append(1)
                            num1.append(-1)
                            fmt1.extend(
                                [r'^\s+(?P<val>Raman Activ Fr=\s?\d --- \s*'
                                 + r'-?\d.*)\s*$',
                                 r'^\s+(?P<val>RamAct Fr=\s?\d+--\s+\d.*)\s*$',
                                 r'^\s+(?P<val>Raman\d Fr=\s?\d+--\s+\d.*)'
                                 + r'\s*$'])
                        lnk1.append(-717)
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
                    fmt1.extend([r'^\s+(?P<val>ROA\d\s+ Fr= ?\d+-- \s*'
                                 + r'-?\d.*)\s*$',
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
    # Check that the DB has been correctly built with all arrays having the
    # same size.
    if not (len(lnk1) == len(key1) == len(sub1) == len(end1) == len(fmt1)
            == len(num1)):
        msg = 'Inconsistency in the construction of parsing data for GLog.'
        raise InternalError(msg)
    if len(lnk1) == 1:
        lnk = lnk1[0]
        key = key1[0]
        sub = sub1[0]
        end = end1[0]
        fmt = fmt1[0]
        num = num1[0]
    else:
        lnk = tuple(lnk1)
        key = tuple(key1)
        sub = tuple(sub1)
        end = tuple(end1)
        fmt = tuple(fmt1)
        num = tuple(num1)

    return lnk, key, sub, fmt, end, num
