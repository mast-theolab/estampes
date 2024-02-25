"""Provide the keys to extract data from log file.

Provides the keys necessary to delimit the data to extract from a
Gaussian log file.
"""

from estampes.base import QLabel
from estampes.parser.gaussian.glog.types import TypeQKwrd


FMT_F_DP = r'\s+-?\d\.\d+D?[+-]?\d{2,3}'


def keys_prp_3xx(qlab: QLabel) -> TypeQKwrd:
    """Provide extractor info for frequency-dependent properties."""
    # Get incident frequency information
    if qlab.kind != 'static':
        if qlab.level == 'H':
            lnk = [716]
            key = [' and normal coordinates:']
            sub = [1]
            end = [lambda s: not s.startswith(' Incident light:')]
            num = [0]
            fmt = [r'^\s+Incident light \(\S+\):\s+(?P<val>-?\d.*)\s*$']
        # -- Anharmonic level
        else:
            lnk = [-717]
            key = [' :    FREQUENCY-DEPENDENT PROPERTIES AND DERIVATIVES    :']
            end = [lambda s: s.startswith('     =====')]
            sub = [2]
            fmt = [r'^ ## INCIDENT FREQUENCY:\s+(?P<val>\S+).*\s+##+\s*$']
            num = [0]
    else:
        lnk = []
        key = []
        end = []
        sub = []
        fmt = []
        num = []
    # Extract properties
    if qlab.label == 301:
        # Polarizability tensor
        if qlab.kind == 'static':
            raise NotImplementedError('Static alpha not yet implemented')
        else:
            if qlab.derord == 0:
                lnk.extend([9999, 717])
                key.extend([
                    ' Dipole polarizability, Alpha (input',
                    ' Frequency-dependent Alpha(-w,w)'
                ])
                sub.extend([2, 8])
                end.extend([
                    lambda s: s.startswith(' --------'),
                    lambda s: s.startswith('     =====')
                    or s.startswith(' ----')
                ])
                num.extend([0, 0])
                fmt.extend([
                    r'^\s+[xyz]{2}\s+(?P<val>-?\d+\.\d+D?[+-]\d+)\s.*$',
                    r'^\s+(?P<p>P0)?\s*\|\s+(?P<val>[^|]+\|'
                    + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                ])
            elif qlab.derord == 1:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Alpha(-w,w)')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P1)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
            elif qlab.derord == 2:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Alpha(-w,w)')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P2)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
            elif qlab.derord == 3:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Alpha(-w,w)')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P3)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
    elif qlab.label == 302:
        # Optical rotation
        if qlab.kind == 'static':
            raise NotImplementedError('Static G not yet implemented')
        else:
            if qlab.derord == 0:
                lnk.append(717)
                key.append(' Frequency-dependent Optical rotations')
                sub.append(8)
                end.append(lambda s: s.startswith('     =====')
                           or s.startswith(' ----'))
                num.append(0)
                fmt.append(
                    r'^\s+(?P<p>P0)?\s*\|\s+(?P<val>[^|]+\|'
                    + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                )
            elif qlab.derord == 1:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Optical rotations')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P1)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
            elif qlab.derord == 2:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Optical rotations')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P2)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
            elif qlab.derord == 3:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Optical rotations')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P3)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)Y|[XZ]) +\|(?:' + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
    elif qlab.label == 304:
        # Electric dipole-electric quadrupole tensor
        if qlab.kind == 'static':
            raise NotImplementedError('Static A not yet implemented')
        else:
            if qlab.derord == 0:
                lnk.append(717)
                key.append(
                    ' Frequency-dependent Dipole-quadrupole polarizability')
                sub.append(8)
                end.append(lambda s: s.startswith('     =====')
                           or s.startswith(' ----'))
                num.append(0)
                fmt.append(
                    r'^\s+(?P<p>P0)?\s*\|\s+(?P<val>[^|]+\|'
                    + r' +(?(p)ZZ|(?:[XY][XYZ]|ZX|ZY)) +\|(?:'
                    + FMT_F_DP + r'){3})\s*$'
                )
            elif qlab.derord == 1:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Dipole-quadrupole '
                               + 'polarizability')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P1)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)ZZ|(?:[XY][XYZ]|ZX|ZY)) +\|(?:'
                        + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
            elif qlab.derord == 2:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Dipole-quadrupole '
                               + 'polarizability')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P2)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)ZZ|(?:[XY][XYZ]|ZX|ZY)) +\|(?:'
                        + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')
            elif qlab.derord == 3:
                if qlab.dercrd == 'Q':
                    lnk.append(717)
                    key.append(' Frequency-dependent Dipole-quadrupole '
                               + 'polarizability')
                    sub.append(8)
                    end.append(lambda s: s.startswith('     =====')
                               or s.startswith(' ----'))
                    num.append(0)
                    fmt.append(
                        r'^\s+(?P<p>P3)?\s*\|\s+(?P<val>[^|]+\|'
                        + r' +(?(p)ZZ|(?:[XY][XYZ]|ZX|ZY)) +\|(?:'
                        + FMT_F_DP + r'){3})\s*$'
                    )
                else:
                    raise NotImplementedError('Support of derivatives NYI.')

    return lnk, key, sub, fmt, end, num
