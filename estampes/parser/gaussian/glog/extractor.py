"""Provide the routines to parse extracted data.

Provides the necessary parsing functions to process data extracted from
a Gaussian log file.
"""

import re
import typing as tp

from estampes.base import QData, \
    ParseKeyError, \
    TypeDGLog, TypeQData, TypeQInfo
from estampes.data.physics import PHYSFACT
from estampes.parser.gaussian.glog.parse_prp import parse_1xx_dat, \
    parse_3xx_dat, parse_13xx_dat, parse_en_dat
from estampes.parser.gaussian.glog.parse_qty import parse_fcdat, \
    parse_ramact_data, parse_trans_str, parse_vptdat, parse_vtrans_data

__ang2au = 1.0 / PHYSFACT.bohr2ang


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
            if not isinstance(datablocks[first], (list, tuple)):
                if (qlabel.rstate == 'c'
                        and 'Excited State' in ''.join(datablocks[first])):
                    realfirst = first + 1
                else:
                    realfirst = first
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
            # We extract the first block, stored in first
            data = [[]]
            for num, line in enumerate(datablocks[first]):
                ov, iops, links = line.split('/')
                if '(' in links:
                    links, jump = links[:-1].split('(')
                    toline = num + 1 + int(jump)
                else:
                    toline = 0
                iops = tuple(iops.split(','))
                for link in links.split(','):
                    fulllink = fmt.format(int(ov), int(link))
                    data[-1].append((fulllink, num+1, iops, toline))
            # Now let us check if more than one job
            # Since we look for Normal termination, we may have an extra blank
            # block last.  We check this:
            if len(datablocks[last]) > 0:
                if not datablocks[last][-1]:
                    del datablocks[last][-1]
            for block in datablocks[last]:
                data.append([])
                for num, line in enumerate(block):
                    ov, iops, links = line.split('/')
                    if '(' in links:
                        links, jump = links[:-1].split('(')
                        toline = num + 1 + int(jump)
                    else:
                        toline = 0
                    iops = tuple(iops.split(','))
                    for link in links.split(','):
                        fulllink = fmt.format(int(ov), int(link))
                        data[-1].append((fulllink, num+1, iops, toline))
            dobjs[qkey].set(data=data)
        elif qlabel.label == 'natoms':
            dobjs[qkey].set(data=int(datablocks[iref][0]))
        # Coordinates
        # -----------
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
                data = {}
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
            dobjs[qkey] = parse_fcdat(qlabel, datablocks[first:last+1],
                                      iref-first)
        # Vibrational transitions
        # -----------------------
        elif qlabel.label in ('vlevel', 'vtrans'):
            dobjs[qkey] = parse_vtrans_data(qlabel, datablocks[first:last+1],
                                            iref-first)
        # Anharmonic Information
        # ----------------------
        elif qlabel.label == 'vptdat':
            dobjs[qkey] = parse_vptdat(qlabel, datablocks[first:last+1])
        # Transition strengths
        # --------------------
        elif qlabel.label in ('dipstr', 'rotstr'):
            dobjs[qkey] = parse_trans_str(qlabel, datablocks[first:last+1],
                                          iref-first)
        # Raman and ROA activity
        # ----------------------
        elif qlabel.label == 'ramact':
            dobjs[qkey] = parse_ramact_data(qlabel, datablocks[first:last+1])
        elif qlabel.label == 'roaact':
            dobjs[qkey] = parse_ramact_data(qlabel, datablocks[first:last+1],
                                            ROA=True)
        # Energy
        # ------
        elif qlabel.label == 1:
            dobjs[qkey] = parse_en_dat(qlabel, datablocks[first:last+1],
                                       iref-first)
        # Electric-field properties
        # -------------------------
        elif qlabel.label in range(101, 200):
            dobjs[qkey] = parse_1xx_dat(qlabel, datablocks[first:last+1])
        # Dynamic (frequency-dependent) properties
        # ----------------------------------------
        elif qlabel.label in range(300, 400):
            dobjs[qkey] = parse_3xx_dat(qlabel, datablocks[first:last+1],
                                        iref-first)
        # Properties 13xx
        # ---------------
        elif qlabel.label in (1300, 1301, 1302, 1303, 1304, 1305):
            dobjs[qkey] = parse_13xx_dat(qlabel, datablocks[first:last+1],
                                         iref-first)
        # State(s)-dependent quantities
        # -----------------------------
        else:
            raise NotImplementedError(f'Quantity {qlabel.label} not supported')
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
