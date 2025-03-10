"""Provide functions to parse specific quantities.

Provides specialized functions to extract data for specific quantities.
"""

import re
import typing as tp

from estampes.base import QData, QLabel, \
    ParseKeyError, QuantityError
from estampes.data.physics import PHYSFACT

__ang2au = 1.0 / PHYSFACT.bohr2ang

xyz2id = {'X': 0, 'Y': 1, 'Z': 2}


def parse_fcdat(qlab: QLabel, dblock: tp.List[str], iref: int = 0) -> QData:
    """Parse extracted data related to Franck-Condon spectroscopy.

    This functions parses the different quantities produced through
    Franck-Condon calculations.

    Notes
    -----
    * This function does not check `qlab.label`.
    """
    dobj = QData(qlab)
    if qlab.kind == 'SimInf':
        def getval(s):
            return s.split(':', maxsplit=1)[-1].strip().title()
        for line in dblock[iref]:
            if 'Temperature' in line:
                if 'not' in line:
                    val = None
                else:
                    val = float(line.split()[-1][:-1])
                dobj.add_field('temp', value=val)
            elif 'framework' in line:
                dobj.add_field('frame', value=getval(line))
            elif 'Spectroscopy' in line:
                dobj.add_field('spec', value=getval(line))
            elif 'Model' in line:
                dobj.add_field('model', value=getval(line))
            elif 'electronic transition moment' in line:
                dobj.add_field('tmom',
                               value=line.split(':', maxsplit=1)[-1].strip())
    elif qlab.kind == 'ExcState':
        line = dblock[iref][0].strip(' .')
        if line.startswith('No'):
            dobj.set(data=0)
        else:
            try:
                dobj.set(data=int(line.split()[-1]))
            except ValueError as err:
                raise QuantityError(
                    'excstate', 'Could not parse index of excited state') \
                    from err
    elif qlab.kind == 'JMat':
        if dblock[-1]:
            i = -1
        else:
            i = 0
        if 'diagonal' in dblock[i]:
            dobj.set(data=[])
        else:
            data = []
            N = 0
            for line in dblock[i]:
                cols = line.split()
                irow = int(cols[0]) - 1
                if irow == N:
                    data.append([])
                    N += 1
                data[irow].extend([float(item.replace('D', 'e'))
                                   for item in cols[1:]])
            dobj.set(data=data)
    elif qlab.kind == 'JMatF':
        if 'diagonal' in dblock[iref]:
            dobj.set(data=[])
        else:
            data = []
            N = 0
            for line in dblock[iref]:
                cols = line.split()
                irow = int(cols[0]) - 1
                if irow == N:
                    data.append([])
                    N += 1
                data[irow].extend([float(item.replace('D', 'e'))
                                   for item in cols[1:]])
            dobj.set(data=data)
    elif qlab.kind == 'KVec':
        data = []
        for line in dblock[iref]:
            data.append(float(line.replace('D', 'e')))
        dobj.set(data=data)
    elif qlab.kind == 'SRAMat':
        data = []
        N = 0
        for line in dblock[iref]:
            cols = line.split()
            irow = int(cols[0]) - 1
            if irow == N:
                data.append([])
                N += 1
            data[irow].extend([float(item.replace('D', 'e'))
                               for item in cols[1:]])
        dobj.set(data=data)
    elif qlab.kind == 'SRBVec':
        data = []
        for line in dblock[iref]:
            data.append(float(line.replace('D', 'e')))
        dobj.set(data=data)
    elif qlab.kind == 'SRCMat':
        data = []
        N = 0
        for line in dblock[iref]:
            cols = line.split()
            irow = int(cols[0]) - 1
            if irow == N:
                data.append([])
                N += 1
            data[irow].extend([float(item.replace('D', 'e'))
                               for item in cols[1:]])
        dobj.set(data=data)
    elif qlab.kind == 'SRDVec':
        data = []
        for line in dblock[iref]:
            data.append(float(line.replace('D', 'e')))
        dobj.set(data=data)
    elif qlab.kind == 'SREMat':
        data = []
        N = 0
        for line in dblock[iref]:
            cols = line.split()
            irow = int(cols[0]) - 1
            if irow == N:
                data.append([])
                N += 1
            data[irow].extend([float(item.replace('D', 'e'))
                               for item in cols[1:]])
        dobj.set(data=data)
    elif qlab.kind == 'Spec':
        # Look for last blocks first, which should contain all blocks:
        # Note: For the vibronic spectrum, X is necessarily the same
        #       over multiple blocks.  We take this for granted to
        #       simplify the processing and the output.
        discard = []
        if len(dblock) > 1:
            nblocks = 0
            nyaxes = 0
            for i, bloc in enumerate(dblock[-1]):
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
            msg = 'Inconsistency in spectral data. This should not happen!'
            raise IndexError(msg)
        if discard:
            discard.reverse()
            for i in discard:
                del dblock[-1][i]
        data = {}
        # First block contains the full initial legend block
        # Necessary for the different parameters
        iref = 0
        data['x'] = []
        if nyaxes == 0:
            nyaxes = len(dblock[iref][0].split()) - 1
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
        iref = -1
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
            for line in dblock[iref][bloc]:
                cols = [float(item.replace('D', 'e'))
                        for item in line.split()]
                if bloc == 0:
                    xax.append(cols[0])
                for i, item in enumerate(cols[1:]):
                    yax[yfmt.format(idy=yoffset+i+1)].append(item)
            yoffset += len(cols) - 1
        for key, val in data.items():
            dobj.add_field(key, value=val)
    elif qlab.kind == 'SpcPar':
        # Look for last blocks first, which should contain all blocks:
        discard = []
        if len(dblock) > 1:
            nblocks = 0
            nyaxes = 0
            for i, bloc in enumerate(dblock[-1]):
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
                del dblock[-1][i]
        # First block contains the full initial legend block
        # Necessary for the different parameters
        iref = 0
        i = 0
        while dblock[iref][i].strip() != 'Legend:':
            i += 1
            if i >= len(dblock[iref]):
                raise IndexError('Legend not found')
        if 'No' in dblock[iref][i-1]:
            dobj.add_field('func', value='stick')
            dobj.add_field('hwhm', value=None)
        elif 'broadening' in dblock[iref][i-2]:
            func = dblock[iref][i-2].split()[-3].lower()
            hwhm = float(dblock[iref][i-1].split()[-2])
            dobj.add_field('func', value=func)
            dobj.add_field('hwhm', value=hwhm)
        else:
            raise IndexError('Unrecognized broadening function.')
        offset = i + 1
        # offset: num. lines for legend block + legend section title
        if nyaxes == 0:
            nyaxes = len(dblock[iref]) - offset - 2
        if nyaxes <= 1:
            # X, Y, Int, don't use numbering
            yfmt = 'y'
        else:
            yfmt = f'y{{idy:0{len(str(nyaxes))}d}}'
        yoffset = 0
        nyi = 0
        for line in dblock[iref][offset:]:
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
                dobj.add_field('unitx', value=txt)
            elif key.strip() == 'Intensity':
                _key = 'I'
                txt = title.rsplit('in ', maxsplit=1)[1].\
                    replace(')', '').strip()
                if dobj.get('func') == 'stick':
                    _desc = 'II:'
                    # Fix a stupidity in the unit in some versions
                    #   for the stick spectrum
                    if txt == 'dm^3.mol^-1.cm^-1':
                        txt = 'dm^3.mol^-1.cm^-2'
                else:
                    _desc = 'I:'
                dobj.add_field('unity', value=_desc + txt)
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
                dobj.add_field(_key, value=title.strip())
            else:
                if _key[0] == 'I':
                    dobj.add_field(_key, value=title.strip())
                else:
                    dobj.add_field(_key, value=title.strip())
        # More than 1 block, read the new ones
        for bloc in range(1, nblocks):
            nyi = 0
            for line in dblock[-1][bloc][1:]:
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
                        dobj.add_field(_key, value=title.strip())
                elif 'intensity' not in res[0]:
                    raise NotImplementedError('Unrecognized legend')
            yoffset += nyi
    elif qlab.kind == 'Conv':
        raise NotImplementedError()
    elif qlab.kind == 'Assign':
        data = {'T': [], 'E': [], 'I': []}
        if qlab.label == 'DipStr':
            qty = 'DS'
        elif qlab.label == 'RotStr':
            qty = 'RS'
        else:
            msg = 'Unrecognized spectroscopy-specific quantity'
            raise IndexError(msg)
        dobj.add_field('other', value=qty)
        data[qty] = []
        for l1, l2 in zip(dblock[0], dblock[-1]):
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
            dobj.add_field(key, value=val)
    elif qlab.kind == 'GeomIS':
        data = []
        # By default, we choose the standard orientation if present
        for line in dblock[iref]:
            data.append([float(item)*__ang2au for item in line.split()])
        dobj.set(data=data)
    elif qlab.kind == 'GeomFS':
        data = []
        # By default, we choose the standard orientation if present
        for line in dblock[iref]:
            data.append([float(item)*__ang2au for item in line.split()])
        dobj.set(data=data)
    elif qlab.kind == 'GeomMS':
        data = []
        # By default, we choose the standard orientation if present
        for line in dblock[iref]:
            data.append([float(item)*__ang2au for item in line.split()])
        dobj.set(data=data)
    elif qlab.kind == 'ExGeom':
        data = []
        # By default, we choose the standard orientation if present
        for line in dblock[iref]:
            data.append([float(item)*__ang2au for item in line.split()])
        dobj.set(data=data)
    elif qlab.kind == 'E(0-0)':
        if len(set(dblock[iref])) > 1:
            msg = 'Excessive information on 0-0 energy.'
            raise ParseKeyError(msg)
        val, unit = dblock[iref][0].split()
        if re.match(r'cm\^.?-1.?\s*', unit, re.I):
            dobj.set(unit='cm^-1')
        else:
            dobj.set(unit='???')
        dobj.set(data=float(val))
    elif qlab.kind == 'RedDim':
        if len(dblock[iref][0].split()) % 3 != 0:
            raise ParseKeyError('Expected reddim structure: num = num')
        nstates = len(dblock[iref][0].split())//3
        data = {}
        for i in range(nstates):
            data['state{}'.format(i+1)] = {}
        for line in dblock[iref]:
            cols = line.split()
            for i in range(nstates):
                data['state{}'.format(i+1)][int(cols[i*3])] = int(cols[i*3+2])
        for key, val in data.items():
            dobj.add_field(key, value=val)
    else:
        raise NotImplementedError('Unknown option for FCDat')

    return dobj


def parse_vptdat(qlab: QLabel, dblock: tp.List[str]) -> QData:
    """Parse extracted data related to VPTx.

    Parses data related to calculations within the vibrational
    perturbation theory.

    Notes
    -----
    * This function does not check `qlab.label`.
    """
    if isinstance(qlab.rstate, tuple):
        raise NotImplementedError('Vibronic VPT2 not yet available.')
    dobj = QData(qlab)
    if qlab.kind == 'CICoef':
        # First, check if reduced-dimensionality used.
        if dblock[-1]:
            eqv = {}
            for line in dblock[-1]:
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
        for line in dblock[0]:
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
        dobj.set(data=data)
    elif qlab.kind == 'XMat':
        maxcols = 5
        nbloc = 0
        dobj.set(unit='cm^{-1}')
        dobj.set(shape='LT')
        data = []
        # store in linear form
        for line in dblock[-1]:
            if '.' not in line:
                nbloc += 1
            else:
                cols = line.split()
                i = int(cols[0])
                if nbloc == 1:
                    data.extend([float(item.replace('D', 'e'))
                                 for item in cols[1:]])
                    if i > 5:
                        data.extend(0.0 for _ in range(maxcols, i))
                else:
                    ioff = i*(i-1)//2 + (nbloc-1)*maxcols
                    ncols = len(cols) - 1
                    data[ioff:ioff+ncols] = [float(item.replace('D', 'e'))
                                             for item in cols[1:]]
        dobj.set(data=data)
    else:
        raise NotImplementedError()

    return dobj


def parse_vtrans_data(qlab: QLabel, dblock: tp.List[str],
                      iref: int = 0) -> QData:
    """Parse extracted data related to vibrational transitions.

    Parses data related to vibrational transition information.
    """
    if isinstance(qlab.rstate, tuple):
        raise NotImplementedError(
            'Pure vibronic transitions not yet supported.')
    msg_noqty = f'Missing quantity "{qlab.label}" in file'
    dobj = QData(qlab)
    if qlab.label == 'vlevel':
        if qlab.kind == 'RR':
            if len(dblock[-1]) != len(dblock[-2]):
                msg = 'Incident frequencies data and transition energies' \
                    + ' do not match.'
                raise ParseKeyError(msg)
            dobj.set(unit='cm-1')
            data = {}
            # Let us check the structure of the incident frequency line.
            # Depending on the version of Gaussian, it can be of two forms:
            # Omega = value cm^-1, Sigma = value // old version
            # Omega = value cm^-1, Gamma = value cm^-1, Sigma = value // new
            # Since the states are the same for any value of Omega and Gamma
            #   we take the first block, so we just need to parse the relevant
            #   part.
            # We take the first element
            key = dblock[-1][0].split(', Sigma')[0]
            count = 0
            for setup_line, trans in zip(dblock[-1], dblock[-2]):
                if setup_line.startswith(key):
                    count += 1
                    data[count] = float(trans)
            dobj.set(data=data)
        else:
            if qlab.level == 'H':
                for i, item in enumerate(reversed(dblock)):
                    if item:
                        iref = -1 - i
                        break
                else:
                    raise ParseKeyError(msg_noqty)
                dobj.set(unit='cm-1')
                i = 0
                data = {}
                for line in dblock[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[i] = float(col)
                        except ValueError:
                            data[i] = float('inf')
                dobj.set(data=data)
            elif qlab.level == 'A':
                if dblock[-1]:
                    iref = -1
                else:
                    iref = 0
                dobj.set(unit='cm-1')
                i = 0
                data = {}
                for line in dblock[iref]:
                    i += 1
                    try:
                        data[i] = float(line)
                    except ValueError:
                        data[i] = float('inf')
                dobj.set(data=data)
            else:
                raise NotImplementedError()
    elif qlab.label == 'vtrans':
        if qlab.kind == 'RR':
            if len(dblock[-1]) != len(dblock[-2]):
                msg = 'Incident frequencies data and transition data ' \
                    + 'do not match.'
                raise ParseKeyError(msg)
            data = {}
            # Let us check the structure of the incident frequency line.
            # Depending on the version of Gaussian, it can be of two forms:
            # Omega = value cm^-1, Sigma = value // old version
            # Omega = value cm^-1, Gamma = value cm^-1, Sigma = value // new
            # Since the states are the same for any value of Omega and Gamma
            #   we take the first block, so we just need to parse the relevant
            #   part.
            # We take the first element
            key = dblock[-1][0].split(', Sigma')[0]
            count = 0
            for setup_line, trans in zip(dblock[-1], dblock[-2]):
                if setup_line.startswith(key):
                    count += 1
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
                    data[count] = tuple(svals)
            dobj.set(data=data)
        else:
            if qlab.level == 'H':
                for i, item in enumerate(reversed(dblock)):
                    if item:
                        iref = -1 - i
                        break
                else:
                    raise ParseKeyError(msg_noqty)
                i = 0
                data = {}
                for line in dblock[iref]:
                    cols = line.strip().split()
                    if cols[-1] in ('active', 'inactive', 'passive'):
                        del cols[-1]
                    for col in cols:
                        i += 1
                        res = col.split('(')
                        data[i] = (((0, 0), ), ((int(res[0]), 1), ))
                dobj.set(data=data)
            elif qlab.level == 'A':
                i = 0
                data = {}
                for line in dblock[iref]:
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
                dobj.set(data=data)
            else:
                raise NotImplementedError()

    return dobj


def parse_trans_str(qlab: QLabel, dblock: tp.List[str],
                    iref: int = 0) -> QData:
    """Parse extracted data related to transition strengths.

    Parses data related to transition strengths, currently:

    * dipole strengths
    * rotational strengths
    """
    msg_noqty = f'Missing quantity "{qlab.label}" in file'
    dobj = QData(qlab)
    if qlab.label == 'dipstr':
        if isinstance(qlab.rstate, tuple):
            Si, Sf = qlab.rstate
            if Si == 0:
                if isinstance(Sf, int):
                    dobj.set(data=float(dblock[iref]))
                    dobj.set(unit='DS:au')
                else:
                    dobj.set(unit='DS:au')
                    dobj.set(data=[float(item) for item in dblock[iref]])
            else:
                raise NotImplementedError()
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError()
            if qlab.level == 'H':
                for i, item in enumerate(reversed(dblock)):
                    if item:
                        iref = -1 - i
                        break
                else:
                    raise ParseKeyError(msg_noqty)
                dobj.set(unit='DS:esu^2.cm^2')
                i = 0
                data = {}
                for line in dblock[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[i] = float(col)*1.0e-40
                        except ValueError:
                            data[i] = float('inf')
                dobj.set(data=data)
            elif qlab.level == 'A':
                dobj.set(unit='DS:esu^2.cm^2')
                i = 0
                data = {}
                for line in dblock[iref]:
                    i += 1
                    try:
                        data[i] = float(line)*1.0e-40
                    except ValueError:
                        data[i] = float('inf')
                dobj.set(data=data)
            else:
                raise NotImplementedError()

    elif qlab.label == 'rotstr':
        if isinstance(qlab.rstate, tuple):
            Si, Sf = qlab.rstate
            if Si == 0:
                if isinstance(Sf, int):
                    dobj.set(
                        data=float(dblock[iref])*1.0e-40)
                    dobj.set(unit='RS:esu^2.cm^2')
                else:
                    dobj.set(unit='RS:esu^2.cm^2')
                    dobj.set(
                        data=[float(item) * 1.0e-40 for item in dblock[iref]])
            else:
                pass
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError()
            if qlab.level == 'H':
                for i, item in enumerate(reversed(dblock)):
                    if item:
                        iref = -1 - i
                        break
                else:
                    raise ParseKeyError(msg_noqty)
                dobj.set(unit='RS:esu^2.cm^2')
                i = 0
                data = {}
                for line in dblock[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[i] = float(col)*1.0e-44
                        except ValueError:
                            data[i] = float('inf')
                dobj.set(data=data)
            elif qlab.level == 'A':
                dobj.set(unit='RS:esu^2.cm^2')
                i = 0
                data = {}
                for line in dblock[iref]:
                    i += 1
                    try:
                        data[i] = float(line)*1.0e-44
                    except ValueError:
                        data[i] = float('inf')
                dobj.set(data=data)
            else:
                raise NotImplementedError()

    else:
        raise QuantityError('strength',
                            'Unsupported type of transition strength')

    return dobj


def parse_ramact_data(qlab: QLabel, dblock: tp.List[str],
                      ROA: bool = False) -> QData:
    """Parse extracted data related to Raman/ROA activity.

    Parses data related to the Raman/Raman optical activity.

    If ROA true, data are related to ROA.
    """
    if isinstance(qlab.rstate, tuple):
        raise QuantityError('Electronic Raman/ROA not available.')
    elif qlab.kind == 'RR':
        dobj = QData(qlab)
        dobj.set(unit='RS:cm^2.sr^-1.mol^-1')
        data = {}
        counts = {}
        line = dblock[-1][0]
        txt = r'Omega =\s*(?P<incfrq>\d+\.\d+) cm.-1\s*,\s+'
        if 'Gamma' in line:
            txt += r'Gamma =\s*(?P<gamma>\d+\.\d+) cm.-1\s*,\s+'
        txt += r'Sigma =\s*(?P<sigma>[-+]?\d\.\d+E?[+-]\d{2,3})'
        pattern = re.compile(txt)
        for line in dblock[-1]:
            res = pattern.match(line).groupdict()
            incfrq = res['incfrq']
            gamma = res.get('gamma')
            sigma = res['sigma']
            if incfrq not in data:
                if gamma is None:
                    data[incfrq] = {}
                    counts[incfrq] = 1
                    d_ptr = data[incfrq]
                    count = counts[incfrq]
                else:
                    data[incfrq] = {gamma: {}}
                    counts[incfrq] = {gamma: 1}
                    d_ptr = data[incfrq][gamma]
                    count = counts[incfrq][gamma]
            else:
                if gamma is not None:
                    if gamma not in data[incfrq]:
                        data[incfrq][gamma] = {}
                        counts[incfrq][gamma] = 1
                    else:
                        counts[incfrq][gamma] += 1
                    d_ptr = data[incfrq][gamma]
                    count = counts[incfrq][gamma]
                else:
                    d_ptr = data[incfrq]
                    counts[incfrq] += 1
                    count = counts[incfrq]
            d_ptr[count] = float(sigma)
        omegas = list(data)
        if gamma is not None:
            # We assume that same gammas used for each omega.
            gammas = list(data[omegas[0]])
        else:
            gammas = None
    else:
        gammas = None
        dobj = QData(qlab)
        data = {}
        for i, item in enumerate(reversed(dblock)):
            if item:
                iref = -1 - i
                break
        else:
            raise ParseKeyError('Missing data for Raman/ROA activity in file')
        # Data in Gaussian log are multiplied by 10^4 for ROA, so we need to
        #   take this into account
        yfactor = 1.0e-4 if ROA else 1.0
        if qlab.level == 'H':
            if iref == -1:
                dobj.set(unit='ROA:Ang^6' if ROA else 'RA:Ang^6')
                an_blk = True
            else:
                dobj.set(unit='ROA:amu.Ang^4' if ROA else 'RA:amu.Ang^4')
                an_blk = False
            # We create the database based on the quantity
            if qlab.kind == 'static':
                i = 0
                for line in dblock[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[i] = float(col)*yfactor
                        except ValueError:
                            data[i] = float('inf')
                omegas = None
            elif an_blk:
                incfrq = None
                omegas = []
                setups = []
                for item in dblock[2]:
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
                    if incfrq not in omegas:
                        omegas.append(incfrq)
                nlines = len(dblock[iref])/len(setups)
                block = 0
                iline = 0
                for line in dblock[iref]:
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
                incfreqs = [item for line in dblock[1]
                            for item in line.split()]
                if qlab.kind == 'dynamic':
                    omegas = incfreqs[:]
                    for freq in incfreqs:
                        data[freq] = {
                            'SCP(180)u': {}, 'SCP(90)z': {}, 'DCPI(180)': {}}
                    iline = 0
                    ifreq = 0
                    ioff = 1
                    i = 0
                    for line in dblock[iref]:
                        iline += 1
                        block = iline % 3  # 3: number of setups
                        if block == 1:
                            if ifreq == len(incfreqs):
                                ifreq = 0
                                ioff += i + 1
                            dfreq = data[incfreqs[ifreq]]
                            ifreq += 1
                        d = dfreq[
                            ('SCP(180)u', 'SCP(90)z', 'DCPI(180)')[block-1]]
                        for i, col in enumerate(line.strip().split()):
                            try:
                                d[ioff+i] = float(col)*yfactor
                            except ValueError:
                                d[ioff+i] = float('inf')
                else:
                    if qlab.kind in incfreqs:
                        ref_freq = incfreqs.index(qlab.kind)
                        omegas = [ref_freq]
                        data[incfreqs[ref_freq]] = {
                            'SCP(180)u': {}, 'SCP(90)z': {}, 'DCPI(180)': {}}
                        iline = 0
                        ifreq = 0
                        iref = 0
                        i = 0
                        for line in dblock[iref]:
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
        elif qlab.level == 'A':
            dobj.set(unit='ROA:Ang^6' if ROA else 'RA:Ang^6')
            if qlab.kind == 'static':
                i = 0
                for line in dblock[iref]:
                    for col in line.strip().split():
                        i += 1
                        try:
                            data[i] = float(col)*yfactor
                        except ValueError:
                            data[i] = float('inf')
                omegas = None
            else:
                incfrq = None
                setups = []
                omegas = []
                for item in dblock[1]:
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
                    if incfrq not in omegas:
                        omegas.append(incfrq)
                nlines = len(dblock[iref])/len(setups)
                block = 0
                iline = 0
                for line in dblock[iref]:
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
    if omegas is not None:
        dobj.add_field('omegas', value=omegas)
    if gammas is not None:
        dobj.add_field('gammas', value=gammas)
    return dobj
