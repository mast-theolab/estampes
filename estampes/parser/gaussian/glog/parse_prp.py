"""Provide functions to parse data for properties.

Provides specialized functions to extract data for specific properties.
"""

import re
import typing as tp

from estampes.base import QData, QLabel, \
    ParseDataError, ParseKeyError, QuantityError
from estampes.data.physics import phys_fact
from estampes.parser.gaussian import g_elquad_LT_to_2D, gfloat

xyz2id = {'X': 0, 'Y': 1, 'Z': 2}


def parse_en_dat(qlab: QLabel, dblock: tp.List[str], iref: int = 0) -> QData:
    """Parse extracted data related to the energy.

    Parses the energy from the extracted data.

    Notes
    -----
    * This function does not check `qlab.label`.
    """
    dobj = QData(qlab)
    if isinstance(qlab.rstate, tuple):
        Si, Sf = qlab.rstate
        if qlab.label == 1:
            if qlab.derord == 0:
                if Si == 0:
                    if isinstance(Sf, int):
                        dobj.set(data=float(dblock[iref]))
                        dobj.set(unit='eV')
                    else:
                        dobj.set(unit='eV')
                        dobj.set(data=[float(item) for item in dblock[iref]])
                else:
                    msg = 'Transition energies different from ground state ' \
                         + 'not implemented.'
                    raise NotImplementedError(msg)
            else:
                msg = 'Derivatives of the transition energies not supported.'
                raise NotImplementedError(msg)
    else:
        if qlab.rstate != 'c':
            # we should check if the state is the right one.
            raise NotImplementedError()

        if qlab.derord == 0:
            def conv(s):
                return gfloat(s.split('=')[1])

            fmt = re.compile(r'E\(?(.*?)\)?\s+')
            N = 2 if dblock[-1] else 1
            nblocks = len(dblock[iref])
            data = {}
            tag = None
            for i in range(iref, iref+N):
                if nblocks > 1:
                    txt = dblock[i][0][0]
                    val = [conv(item[0]) for item in dblock[i]]
                else:
                    txt = dblock[i][0]
                    val = conv(txt)
                res = fmt.search(txt)
                if res is None:
                    msg = 'Unsupported energy format'
                    raise ParseKeyError(msg)
                tag = res.group(1)
                data[tag] = val
            dobj.set(data=data[tag])
            dobj.set(unit='Eh')
            for key, val in data.items():
                dobj.add_field(key, value=val)
        elif qlab.derord == 2:
            maxcols = 5
            nbloc = 0
            dobj.set(unit='Eh.a0^{-2}')
            dobj.set(shape='LT')
            data = []
            # store in linear form
            for line in dblock[iref]:
                if '.' not in line:
                    nbloc += 1
                else:
                    cols = line.split()
                    i = int(cols[0])
                    if nbloc == 1:
                        data.extend([gfloat(item) for item in cols[1:]])
                        if i > 5:
                            data.extend(0.0 for _ in range(maxcols, i))
                    else:
                        ioff = i*(i-1)//2 + (nbloc-1)*maxcols
                        ncols = len(cols) - 1
                        data[ioff:ioff+ncols] \
                            = [gfloat(item) for item in cols[1:]]
            dobj.set(data=data)
        elif qlab.derord == 3:
            if qlab.dercrd == 'Q':
                dobj.set(unit='Eh.u^{-3/2}.a0^{-3}')
            elif qlab.dercrd == 'QRED':
                dobj.set(unit='cm^{-1}')
            else:
                dobj.set(unit='N/A')
            data = {}
            # store in linear form
            for line in dblock[iref]:
                cols = line.split()
                modes = [int(i) for i in cols[:3]]
                value = float(cols[-1])
                if modes[0] not in data:
                    data[modes[0]] = {modes[1]: {modes[2]: value}}
                elif modes[1] not in data[modes[0]]:
                    data[modes[0]][modes[1]] = {modes[2]: value}
                else:
                    data[modes[0]][modes[1]][modes[2]] = value
            dobj.set(data=data)
        elif qlab.derord == 4:
            if qlab.dercrd == 'Q':
                dobj.set(unit='Eh.u^{-3/2}.a0^{-3}')
            elif qlab.dercrd == 'QRED':
                dobj.set(unit='cm^{-1}')
            else:
                dobj.set(unit='N/A')
            data = {}
            # store in linear form
            for line in dblock[iref]:
                cols = line.split()
                modes = [int(i) for i in cols[:4]]
                value = float(cols[-1])
                if modes[0] not in data:
                    data[modes[0]] = {modes[1]: {modes[2]: {modes[3]: value}}}
                elif modes[1] not in data[modes[0]]:
                    data[modes[0]][modes[1]] = {modes[2]: {modes[3]: value}}
                elif modes[2] not in data[modes[0]][modes[1]]:
                    data[modes[0]][modes[1]][modes[2]] = {modes[3]: value}
                else:
                    data[modes[0]][modes[1]][modes[2]][modes[3]] = value
            dobj.set(data=data)
        else:
            raise NotImplementedError('Higher-order energy deriv. NYI.')

    return dobj


def parse_1xx_dat(qlab: QLabel, dblock: tp.List[str]) -> QData:
    """Parse extracted data related to 1xx properties/quantities."""
    dobj = QData(qlab)
    if qlab.label == 101:
        if isinstance(qlab.rstate, tuple):
            Si, Sf = qlab.rstate
            if qlab.derord == 0:
                dobj.set(unit='e.a0')
                if Si == 0:
                    if isinstance(Sf, int):
                        val = [float(item) for item in dblock[-1][0].split()]
                    elif Sf == 'a':
                        val = {}
                        for i, line in enumerate(dblock[-1]):
                            val[i+1] = [float(item) for item in line.split()]
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                elif isinstance(Si, int):
                    if isinstance(Sf, int):
                        val = [float(item)
                               for item in dblock[-1][0].split()[2:]]
                    elif Sf == 'a':
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            x, y, z = float(cols[2]), float(cols[3]), \
                                float(cols[4])
                            val[j] = [x, y, z]
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                elif Si == 'a':
                    if isinstance(Sf, int):
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            x, y, z = float(cols[2]), float(cols[3]), \
                                float(cols[4])
                            val[i] = [x, y, z]
                    elif Sf == 'a':
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            x, y, z = float(cols[2]), float(cols[3]), \
                                float(cols[4])
                            if i not in val:
                                val[i] = {}
                            val[i][j] = [x, y, z]
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                else:
                    raise NotImplementedError(
                        'Unsupported initial state spec.')
                dobj.set(data=val)
            else:
                raise NotImplementedError(
                    'Electronic transition dipole moments derivatives NYI')
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError()

            if qlab.derord == 0:
                if qlab.rstate == 'c':
                    val = [float(item)/phys_fact('au2Deb') for item in
                           dblock[-1][0].split()[1::2]]
                    dobj.set(data=val)
                    dobj.set(unit='e.a0')
                else:
                    raise NotImplementedError()
            else:
                raise NotImplementedError()
    elif qlab.label == 102:
        if isinstance(qlab.rstate, tuple):
            Si, Sf = qlab.rstate
            if qlab.derord == 0:
                dobj.set(unit='hbar.e/me')
                if Si == 0:
                    if isinstance(Sf, int):
                        val = [float(item) for item in dblock[-1][0].split()]
                    elif Sf == 'a':
                        val = {}
                        for i, line in enumerate(dblock[-1]):
                            val[i+1] = [float(item) for item in line.split()]
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                elif isinstance(Si, int):
                    if isinstance(Sf, int):
                        val = [float(item)
                               for item in dblock[-1][0].split()[2:]]
                    elif Sf == 'a':
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            x, y, z = float(cols[2]), float(cols[3]), \
                                float(cols[4])
                            val[j] = [x, y, z]
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                elif Si == 'a':
                    if isinstance(Sf, int):
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            x, y, z = float(cols[2]), float(cols[3]), \
                                float(cols[4])
                            val[i] = [x, y, z]
                    elif Sf == 'a':
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            x, y, z = float(cols[2]), float(cols[3]), \
                                float(cols[4])
                            if i not in val:
                                val[i] = {}
                            val[i][j] = [x, y, z]
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                else:
                    raise NotImplementedError(
                        'Unsupported initial state spec.')
                dobj.set(data=val)
            else:
                raise NotImplementedError(
                    'Electronic transition dipole moments derivatives NYI')
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError(f'Prp {qlab.label} NYI')

            if qlab.derord == 0:
                raise QuantityError('Magnetic dipole moment not available')
            else:
                raise NotImplementedError(f'Prp {qlab.label} NYI')
    elif qlab.label == 107:
        if isinstance(qlab.rstate, tuple):
            Si, Sf = qlab.rstate
            if qlab.derord == 0:
                dobj.set(unit='e.a0^2')
                if Si == 0:
                    if isinstance(Sf, int):
                        val = g_elquad_LT_to_2D(
                            [float(item) for item in dblock[-1][0].split()])
                    elif Sf == 'a':
                        val = {}
                        for i, line in enumerate(dblock[-1]):
                            val[i+1] = g_elquad_LT_to_2D(
                                [float(item) for item in line.split()])
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                elif isinstance(Si, int):
                    if isinstance(Sf, int):
                        val = g_elquad_LT_to_2D(
                            [float(item)
                             for item in dblock[-1][0].split()[2:]])
                    elif Sf == 'a':
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            val[j] = g_elquad_LT_to_2D(
                                [float(item) for item in cols[2:]])
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                elif Si == 'a':
                    if isinstance(Sf, int):
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            val[i] = g_elquad_LT_to_2D(
                                [float(item) for item in cols[2:]])
                    elif Sf == 'a':
                        val = {}
                        for line in dblock[-1]:
                            cols = line.split()
                            j, i = int(cols[0]), int(cols[1])
                            if i not in val:
                                val[i] = {}
                            val[i][j] = g_elquad_LT_to_2D(
                                [float(item) for item in cols[2:]])
                    else:
                        raise NotImplementedError(
                            'Unsupported final state spec.')
                else:
                    raise NotImplementedError(
                        'Unsupported initial state spec.')
                dobj.set(data=val)
            else:
                raise NotImplementedError(
                    'Electronic transition dipole moments derivatives NYI')
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError(f'Prp {qlab.label} NYI')

            if qlab.derord == 0:
                raise NotImplementedError(f'Prp {qlab.label} NYI')
            else:
                raise NotImplementedError(f'Prp {qlab.label} NYI')
    else:
        raise NotImplementedError(f'Quantity {qlab.label} not yet supported.')

    return dobj


def parse_3xx_dat(qlab: QLabel, dblock: tp.List[str], iref: int = 0) -> QData:
    """Parse extracted data related to 3xx properties/quantities.

    Notes
    -----
    * All values are computed, even if ~ 0.
    * There is no exclusion for active/inactive/frozen modes.
    """
    def list_incfrq(items):
        incfrqs = []
        # Incident frequencies can be given in items as:
        # 1 frequency at a time over multiple elements (anharm.)
        # multiple frequencies stored in one (or more) items.
        # We do a double loop on elements and content for this reason.
        for item in items:
            for element in item.strip().split():
                if element not in incfrqs:
                    incfrqs.append(element)
        return incfrqs

    def split3(block):
        return [gfloat(item) for item in block.split()]

    dobj = QData(qlab)
    if qlab.label == 300:
        data = {'data': [], 'keys': []}
        if qlab.level != 'A':
            unit = 'cm^-1'
            for line in dblock[iref]:
                for key in line.split():
                    if key not in data['keys']:
                        data['keys'].append(key)
                        data['data'].append(float(key))
        else:
            unit = dblock[iref][0].split()[-1]
            for line in dblock[iref]:
                key, unit = line.split()
                if key not in data['keys']:
                    data['keys'].append(key)
                    data['data'].append(float(key))
        dobj.set(unit=unit.lower())
        dobj.set(data=data['data'])
        dobj.add_field('keys', value=data['keys'])
    elif qlab.label in (301, 302, 303):
        if isinstance(qlab.rstate, tuple):
            raise NotImplementedError('Electronic transition dip. moment NYI.')
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError()

            # Extract unique incident frequencies if relevant
            if qlab.kind in ('static', 0):
                raise NotImplementedError('Static freq-dep. property NYI.')
            else:
                incfrqs = list_incfrq(dblock[iref])
                if qlab.derord == 0:
                    raise NotImplementedError(
                        'Reference freq-dep. property tensor NYI.')
                elif qlab.derord == 1:
                    if qlab.dercrd == 'X':
                        data = {freq: [] for freq in incfrqs}
                        n_incfrqs = len(incfrqs)
                        if len(dblock[-1]) % n_incfrqs != 0:
                            raise ParseDataError(
                                'd alpha/d X',
                                'Unknown data structure in Gaussian log file')
                        # NOTE: we assume the data stored on max. 5 cols.
                        for incfrq, block in zip(incfrqs,
                                                 dblock[-1][-n_incfrqs:]):
                            ioff = 0
                            dalp_dx = []
                            for iblk in range(0, len(block), 4):
                                dax = [gfloat(item) for item in
                                       block[iblk+1].strip().split()[1:]]
                                day = [gfloat(item) for item in
                                       block[iblk+2].strip().split()[1:]]
                                daz = [gfloat(item) for item in
                                       block[iblk+3].strip().split()[1:]]
                                ncols = len(block[iblk].strip().split())
                                if ioff == 0:
                                    if ncols < 3:
                                        raise ParseDataError(
                                            'd alpha/d X',
                                            'Incorrect number of columns in '
                                            + f'block {iblk+1}')
                                    data[incfrq].append([
                                        dax[0:3], day[0:3], daz[0:3]
                                    ])
                                    if ncols > 3:
                                        dalp_dx = [dax[3:], day[3:], daz[3:]]
                                        ioff = ncols - 3
                                else:
                                    if ncols < 3 - ioff:
                                        raise ParseDataError(
                                            'd alpha/d X',
                                            'Missing data to build tensor')
                                    dalp_dx[0].extend(dax[:3-ioff])
                                    dalp_dx[1].extend(day[:3-ioff])
                                    dalp_dx[2].extend(daz[:3-ioff])
                                    data[incfrq].append(dalp_dx)
                                    nrest = ncols - (3 - ioff)
                                    if nrest > 3:
                                        # only possible case: 1 + 3 + 1
                                        data[incfrq].append([
                                            dax[-nrest:-nrest+3],
                                            day[-nrest:-nrest+3],
                                            daz[-nrest:-nrest+3],
                                        ])
                                        dalp_dx = [
                                            dax[-nrest+3:],
                                            day[-nrest+3:],
                                            daz[-nrest+3:]]
                                        ioff = nrest - 3
                                    elif nrest == 3:
                                        # exactly 3 columns remaining
                                        data[incfrq].append([
                                            dax[-nrest:],
                                            day[-nrest:],
                                            daz[-nrest:],
                                        ])
                                        ioff = 0
                                    else:
                                        dalp_dx = [
                                            dax[-nrest:],
                                            day[-nrest:],
                                            daz[-nrest:]]
                                        ioff = nrest
                        dobj.set(data=data)
                        dobj.set(unit='a0^2')
                    elif qlab.dercrd == 'Q':
                        data = {freq: {} for freq in incfrqs}
                        ioff = -1
                        for ln, line in enumerate(dblock[-1]):
                            cols = line.split('|')
                            if cols[1].strip() == 'Y':
                                i = int(cols[0])
                                if i == 1:
                                    ioff += 1
                                data[incfrqs[ioff]][i] = [
                                    split3(dblock[-1][ln-1].split('|')[-1]),
                                    split3(cols[-1]),
                                    split3(dblock[-1][ln+1].split('|')[-1])
                                ]
                        if qlab.kind == 'dynamic':
                            dobj.set(data=data)
                        else:
                            if qlab.kind in data:
                                dobj.set(data=data[qlab.kind])
                            else:
                                raise ParseKeyError(
                                    'Missing incident frequency')
                    else:
                        raise NotImplementedError()
                elif qlab.derord == 2:
                    if qlab.dercrd == 'X':
                        raise NotImplementedError('d2 P/d X2 NYI.')
                    elif qlab.dercrd == 'Q':
                        data = {freq: {} for freq in incfrqs}
                        ioff = -1
                        for ln, line in enumerate(dblock[-1]):
                            cols = line.split('|')
                            if cols[1].strip() == 'Y':
                                i, j = (int(item) for item in cols[0].split())
                                if i == 1 and j == 1:
                                    ioff += 1
                                if i not in data[incfrqs[ioff]]:
                                    data[incfrqs[ioff]][i] = {}
                                data[incfrqs[ioff]][i][j] = [
                                    split3(dblock[-1][ln-1].split('|')[-1]),
                                    split3(cols[-1]),
                                    split3(dblock[-1][ln+1].split('|')[-1])
                                ]
                        if qlab.kind == 'dynamic':
                            dobj.set(data=data)
                        else:
                            if qlab.kind in data:
                                dobj.set(data=data[qlab.kind])
                            else:
                                raise ParseKeyError(
                                    'Missing incident frequency')
                    else:
                        raise NotImplementedError()
                elif qlab.derord == 3:
                    if qlab.dercrd == 'X':
                        raise NotImplementedError('d3 P/d X3 NYI.')
                    elif qlab.dercrd == 'Q':
                        data = {freq: {} for freq in incfrqs}
                        ioff = -1
                        for ln, line in enumerate(dblock[-1]):
                            cols = line.split('|')
                            if cols[1].strip() == 'Y':
                                i, j, k = \
                                    (int(item) for item in cols[0].split())
                                if i == 1 and j == 1 and k == 1:
                                    ioff += 1
                                if i not in data[incfrqs[ioff]]:
                                    data[incfrqs[ioff]][i] = {}
                                if j not in data[incfrqs[ioff]][i]:
                                    data[incfrqs[ioff]][i][j] = {}
                                data[incfrqs[ioff]][i][j][k] = [
                                    split3(dblock[-1][ln-1].split('|')[-1]),
                                    split3(cols[-1]),
                                    split3(dblock[-1][ln+1].split('|')[-1])
                                ]
                        if qlab.kind == 'dynamic':
                            dobj.set(data=data)
                        else:
                            if qlab.kind in data:
                                dobj.set(data=data[qlab.kind])
                            else:
                                raise ParseKeyError(
                                    'Missing incident frequency')
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()
    elif qlab.label == 304:
        if isinstance(qlab.rstate, tuple):
            raise NotImplementedError('Electronic transition dip. moment NYI.')
        else:
            if qlab.rstate != 'c':
                # we should check if the state is the right one.
                raise NotImplementedError()

            # Extract unique incident frequencies if relevant
            if qlab.kind in ('static', 0):
                raise NotImplementedError(
                    'Static dipole-quadrupole property NYI.')
            else:
                incfrqs = list_incfrq(dblock[iref])
                if qlab.derord == 0:
                    raise NotImplementedError(
                        'Reference freq-dep. property tensor NYI.')
                elif qlab.derord == 1:
                    if qlab.dercrd == 'X':
                        raise NotImplementedError('d A/d X NYI.')
                    elif qlab.dercrd == 'Q':
                        data = {freq: {} for freq in incfrqs}
                        ioff = -1
                        for ln, line in enumerate(dblock[-1]):
                            cols = line.split('|')
                            if cols[1].strip() == 'ZZ':
                                i = int(cols[0])
                                if i == 1:
                                    ioff += 1
                                xy = split3(dblock[-1][ln-5].split('|')[-1])
                                xz = split3(dblock[-1][ln-4].split('|')[-1])
                                yz = split3(dblock[-1][ln-3].split('|')[-1])
                                xx = split3(dblock[-1][ln-2].split('|')[-1])
                                yy = split3(dblock[-1][ln-1].split('|')[-1])
                                zz = split3(cols[-1])
                                data[incfrqs[ioff]][i] = [
                                    [  # x, ...
                                        [xx[0], xy[0], xz[0]],
                                        [xy[0], yy[0], yz[0]],
                                        [xz[0], xy[0], zz[0]]],
                                    [  # y, ...
                                        [xx[1], xy[1], xz[1]],
                                        [xy[1], yy[1], yz[1]],
                                        [xz[1], xy[1], zz[1]]],
                                    [  # z, ...
                                        [xx[2], xy[2], xz[2]],
                                        [xy[2], yy[2], yz[2]],
                                        [xz[2], xy[2], zz[2]]]
                                ]
                        if qlab.kind == 'dynamic':
                            dobj.set(data=data)
                        else:
                            if qlab.kind in data:
                                dobj.set(data=data[qlab.kind])
                            else:
                                raise ParseKeyError(
                                    'Missing incident frequency')
                    else:
                        raise NotImplementedError()
                elif qlab.derord == 2:
                    if qlab.dercrd == 'X':
                        raise NotImplementedError('d2 P/d X2 NYI.')
                    elif qlab.dercrd == 'Q':
                        data = {freq: {} for freq in incfrqs}
                        ioff = -1
                        for ln, line in enumerate(dblock[-1]):
                            cols = line.split('|')
                            if cols[1].strip() == 'ZZ':
                                i, j = (int(item) for item in cols[0].split())
                                if i == 1 and j == 1:
                                    ioff += 1
                                if i not in data[incfrqs[ioff]]:
                                    data[incfrqs[ioff]][i] = {}
                                xy = split3(dblock[-1][ln-5].split('|')[-1])
                                xz = split3(dblock[-1][ln-4].split('|')[-1])
                                yz = split3(dblock[-1][ln-3].split('|')[-1])
                                xx = split3(dblock[-1][ln-2].split('|')[-1])
                                yy = split3(dblock[-1][ln-1].split('|')[-1])
                                zz = split3(cols[-1])
                                data[incfrqs[ioff]][i][j] = [
                                    [  # x, ...
                                        [xx[0], xy[0], xz[0]],
                                        [xy[0], yy[0], yz[0]],
                                        [xz[0], xy[0], zz[0]]],
                                    [  # y, ...
                                        [xx[1], xy[1], xz[1]],
                                        [xy[1], yy[1], yz[1]],
                                        [xz[1], xy[1], zz[1]]],
                                    [  # z, ...
                                        [xx[2], xy[2], xz[2]],
                                        [xy[2], yy[2], yz[2]],
                                        [xz[2], xy[2], zz[2]]]
                                ]
                        if qlab.kind == 'dynamic':
                            dobj.set(data=data)
                        else:
                            if qlab.kind in data:
                                dobj.set(data=data[qlab.kind])
                            else:
                                raise ParseKeyError(
                                    'Missing incident frequency')
                    else:
                        raise NotImplementedError()
                elif qlab.derord == 3:
                    if qlab.dercrd == 'X':
                        raise NotImplementedError('d3 P/d X3 NYI.')
                    elif qlab.dercrd == 'Q':
                        data = {freq: {} for freq in incfrqs}
                        ioff = -1
                        for ln, line in enumerate(dblock[-1]):
                            cols = line.split('|')
                            if cols[1].strip() == 'ZZ':
                                i, j, k = \
                                    (int(item) for item in cols[0].split())
                                if i == 1 and j == 1 and k == 1:
                                    ioff += 1
                                xy = split3(dblock[-1][ln-5].split('|')[-1])
                                xz = split3(dblock[-1][ln-4].split('|')[-1])
                                yz = split3(dblock[-1][ln-3].split('|')[-1])
                                xx = split3(dblock[-1][ln-2].split('|')[-1])
                                yy = split3(dblock[-1][ln-1].split('|')[-1])
                                zz = split3(cols[-1])
                                if i not in data[incfrqs[ioff]]:
                                    data[incfrqs[ioff]][i] = {}
                                if j not in data[incfrqs[ioff]][i]:
                                    data[incfrqs[ioff]][i][j] = {}
                                data[incfrqs[ioff]][i][j][k] = [
                                    [  # x, ...
                                        [xx[0], xy[0], xz[0]],
                                        [xy[0], yy[0], yz[0]],
                                        [xz[0], xy[0], zz[0]]],
                                    [  # y, ...
                                        [xx[1], xy[1], xz[1]],
                                        [xy[1], yy[1], yz[1]],
                                        [xz[1], xy[1], zz[1]]],
                                    [  # z, ...
                                        [xx[2], xy[2], xz[2]],
                                        [xy[2], yy[2], yz[2]],
                                        [xz[2], xy[2], zz[2]]]
                                ]
                        if qlab.kind == 'dynamic':
                            dobj.set(data=data)
                        else:
                            if qlab.kind in data:
                                dobj.set(data=data[qlab.kind])
                            else:
                                raise ParseKeyError(
                                    'Missing incident frequency')
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()
    else:
        raise NotImplementedError(f'Quantity {qlab.label} not yet supported.')

    return dobj


def parse_13xx_dat(qlab: QLabel, dblock: tp.List[str], iref: int = 0) -> QData:
    """Parse extracted data related to 13xx properties/quantities."""
    if isinstance(qlab.rstate, tuple):
        raise NotImplementedError('Transition moment of 13xx NYI.')
    dobj = QData(qlab)
    if qlab.label == 1300:
        unit = dblock[iref][0].split()[-1]
        dobj.set(unit=unit.lower())
        data = {'data': [], 'keys': []}
        for line in dblock[iref]:
            key, unit = line.split()
            if key not in data['keys']:
                data['keys'].append(key)
                data['data'].append(float(key))
        dobj.set(data=data['data'])
        dobj.add_field('keys', value=data['keys'])
        return dobj

    if qlab.label in (1301, 1302, 1303, 1304, 1305):
        if len(dblock[-1]) != len(dblock[-2]):
            msg = 'Incident frequencies data and properties moments do not ' \
                + 'match.'
            raise ParseKeyError(msg)
    else:
        raise QuantityError('Unsupported property in parse_13xx.')
    counts = {}
    data = {}
    dosym = False
    for incfrq, dtens in zip(dblock[-1], dblock[-2]):
        if incfrq not in data:
            data[incfrq] = {}
            counts[incfrq] = 1
        else:
            counts[incfrq] += 1
        if qlab.label in (1301, 1302, 1303):
            tensor = [
                [None, None, None],
                [None, None, None],
                [None, None, None]
            ]
        elif qlab.label in (1304, 1305):
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
