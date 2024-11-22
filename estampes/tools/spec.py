"""Toolbox for spectra

Module providing tools to manipulate data relative to spectra.

Developer notes:
The warning about eval is not fixed to offer some flexibility in
handling units with numerical factors for less standard conversions
which are not directly considered (e.g., 10^-40 esu^2 cm^2).
The calling functions should make sure that the string is compatible
with a unit.
"""

from math import ceil, log, pi
import typing as tp

from estampes.base import ArgumentError
from estampes.tools.char import convert_expr
from estampes.tools.math import f_gauss, f_lorentz
# from estampes.data.property import property_units as punits
from estampes.data.physics import PHYSFACT, PHYSCNST, phys_fact


# ==============
# Module Methods
# ==============

def broaden(xval: tp.Sequence[float],
            yval: tp.Sequence[float],
            xaxis: tp.Sequence[float],
            funcname: str,
            hwhm: float,
            yfactor: float = 1.0,
            xfunc: tp.Optional[tp.Callable[[float], float]] = None,
            ycorr: tp.Optional[tp.Callable[[float], float]] = None,
            ynorm: bool = False,
            truncate: bool = False) -> tp.List[float]:
    """Broadens a stick spectra with a given broadening function.

    Performs the broadening operation given a list of transitions stored in
    two separate arrays, `xval` and `yval`.
    Returns the broadened spectrum for each point in `xaxis`.

    `yaxis` is computed as:

    .. math:: y_a(i) = sum_j y_f*y_v(j)*f(x_a(i))*c(x_v(j)*g(x_a(i)-x_v(j))

    * :math:`y_a`: `yaxis`
    * :math:`y_f`: `yfactor`
    * :math:`y_v`: `yval`
    * :math:`x_a`: `xaxis`
    * :math:`x_v`: `xval`
    * :math:`f`: xfunc
    * :math:`c`: ycorr
    * :math:`g`: broadening function, related to `funcname`

    Parameters
    ----------
    xval
        List of X values.
    yval
        List of Y values.
    xaxis
        X axis of the broadened spectrum.
    funcname
        Broadening function: gaussian, lorentzian.
    hwhm
        Half-width at half-maximum (in unit of X).
    yfactor
        Scaling factor for the Y axis.
    xfunc
        Transformation of X to be included into the final Y value.
    ycorr
        Correction of Y values with respect to X values.
    ynorm
        Normalize Y (in this case, `yfactor` is ignored).
    truncate
        Truncate the broadening (for speed up).
        Note: `xaxis` must be linear and regular (delta x constant).

    Returns
    -------
    list
        Y axis with the broadened bands.

    Raises
    ------
    IndexError
        `xaxis` is empty.
    """
    supported_func = {'gaussian': f_gauss, 'lorentzian': f_lorentz}
    # Transform broadening function to be elemental
    func = supported_func[funcname]
    npoints = len(xaxis)
    if npoints == 0:
        raise IndexError('New X axis is empty')
    yaxis = [0.0 for _ in range(npoints)]
    if not truncate:
        if ycorr is None:
            for i in range(npoints):
                for x0, y0 in zip(xval, yval):
                    yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)
        else:
            for i in range(npoints):
                for x0, y0 in zip(xval, yval):
                    yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)*ycorr(x0)
    else:
        # hwtimes: Function is computer over an interval of 2*hwhm*hwtimes
        # 3.5*hwhm represents 99.7% of the Gaussian area
        # For lorentzian, the same coverage requires 212*hwhm, truncated to 15.
        hwtimes = {'gaussian': 3.5, 'lorentzian': 15}[funcname]
        xmin = xaxis[0]
        dx = abs(xaxis[-1] - xaxis[0])/(npoints-1)
        nhwpoints = int(ceil(hwhm*hwtimes/dx))
        if xmin > xaxis[-1]:
            # X axis is reversed
            xmin = xaxis[-1]
            for x0, y0 in zip(xval, yval):
                i0 = -int((x0 - xmin)/dx) - 1
                imin = max(i0 - nhwpoints, -npoints)
                imax = min(i0 + nhwpoints, -1)
                if ycorr is None:
                    for i in range(imin, imax+1):
                        yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)
                else:
                    for i in range(imin, imax+1):
                        yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor) \
                            * ycorr(x0)
        else:
            for x0, y0 in zip(xval, yval):
                i0 = int((x0 - xmin)/dx)
                imin = max(i0 - nhwpoints, 0)
                imax = min(i0 + nhwpoints, npoints-1)
                if ycorr is None:
                    for i in range(imin, imax+1):
                        yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)
                else:
                    for i in range(imin, imax+1):
                        yaxis[i] += func(xaxis[i], hwhm, x0, y0*yfactor)\
                            * ycorr(x0)
    if xfunc is not None:
        for i in range(npoints):
            yaxis[i] *= xfunc(xaxis[i])
    if ynorm:
        ymax = max(yaxis, key=abs)
        for i in range(npoints):
            yaxis[i] /= ymax
    return yaxis


def convert_x(unit_to: str,
              unit_from: str,
              xaxis: tp.Sequence[float],
              ) -> tp.List[float]:
    """Returns the converted X unit

    Returns the converted X axis from `unit_from` to `unit_to`.
    The units are strings.

    Parameters
    ----------
    unit_to
        Final X unit.
    unit_from
        Original X unit.
    xaxis
        Original X axis values.

    Returns
    -------
    list
        Converted values of X.

    Raises
    ------
    ZeroDivisionError
        The conversion requires an inversion, but `xaxis` contains
        null or near-null values.
    """
    raise NotImplementedError('X conversion not yet available')


def convert_y(specabbr: str,
              unit_to: str,
              unit_from: str,
              **subopts: tp.Dict[str, tp.Any]
              ) -> tp.Tuple[float,
                            tp.Optional[tp.Callable[[float], float]],
                            tp.Optional[tp.Callable[[float], float]]]:
    """Returns parameters for the conversion of Y units.

    Returns the scaling factor and the power of X needed to convert the
    Y unit from `unit_from` to `unit_to`.
    The units are strings composed of two parts, separated by colon:

    * a description: DS, RS, RA, I, II...
    * a unit: km/mol, cgs...

    A default unit may be available, in which case, the unit can be
    provided as: "desc:"

    Parameters
    ----------
    specabbr
        Spectroscopy name (abbreviated version).
    unit_to
        Final Y unit.
    unit_from
        Original Y unit.
    subopts
        Sub-options, dependent on quantity:

        'incfreq'
            incident frequency for Raman, ROA (in cm^-1)
            if not present or None, use def. 532 nm.

    Returns
    -------
    float
        Scaling factor on the Y unit.
    xfunc
        Function accounting for the contribution of X to final Y value.
    ycorr
        Function to correct/convert X values before broadening.
        For instance, to convert hybrid units used in Raman.

    Raises
    ------
    ArgumentError
        Error in sub-options.
    IndexError
        Unrecognized spectroscopy or unit.
    NotImplementedError
        Spectroscopy or unit not yet available.

    Notes
    -----
    * Recognized unit descriptions:

        * `I` - intensity
        * `II` - integrated intensity
        * `DS` - dipole strength
        * `RS` - rotatory strength
        * `RA` - Raman activity
        * `ROA` - Raman optical activity

    * Units are case-sensitive to avoid ambiguities.
      "/" and "." must be indicated correctly.

        * `M-1` or `/M` -> dm^3/mol
        * `M-1.cm-1` or `/M/cm` -> dm^3/mol/cm
        * `M-1.cm-2` or `/M/cm2` -> dm^3/mol/cm^2

    """
    msgNYI = f'For {specabbr}: conversion from "{unit_from}" to "{unit_to}"' +\
        ' not yet implemented.'
    msgNA = f'Conversion from "{unit_from}" to "{unit_to}" not recognized ' +\
        f'for {specabbr}'
    msgVal = 'Wrong value for option {subkey} in spectroscopy {specabbr}'
    UNIT_EPS = {  # Epsilon: Molar absorption coefficient
        '/M/cm': ('M-1.cm-1', '/M/cm', 'dm3.mol-1.cm-1', 'dm3/mol/cm',
                  'L.mol-1.cm-1', 'L/mol/cm', '')
    }
    UNIT_II = {  # Integrated intensity
        '/M/cm2': ('M-1.cm-2', '/M/cm2', 'dm3.mol-1.cm-2', 'dm3/mol/cm2',
                   'L.mol-1.cm-2', 'L/mol/cm2'),
        'km/M': ('km/mol', 'km/M', 'km.mol-1', 'km.M-1')
    }
    UNIT_DS = {  # Dipole strength
        'esu2.cm2': ('statC2.cm2', 'esu2.cm2'),
        'au': ('au', 'e2.a02')
    }
    UNIT_RS = {  # Rotatory strength
        'esu2.cm2': ('statA2.cm2', 'esu2.cm2')
    }
    UNIT_RDSG = {  # Raman differential sigma (cross section)
        'cm3/mol/sr': ('cm3/mol/sr', 'cm3.mol-1.sr-1'),
        'm2.cm/mol/sr': ('m2.cm/mol/sr', 'm2.cm.mol-1.sr-1'),
    }
    UNIT_RA = {  # Raman activity
        'Ang4': ('AA4', 'Ang4'),
        'Ang6': ('AA6', 'Ang6'),
        'bohr6': ('au6', 'bohr6'),
        'amu.Ang4': ('amu.AA4', 'amu.Ang4'),
    }
    UNIT_ROA = {  # Raman optical activity
        'Ang4': ('AA4', 'Ang4'),
        'Ang6': ('AA6', 'Ang6'),
        'bohr6': ('au6', 'bohr6'),
        'amu.Ang4': ('amu.AA4', 'amu.Ang4'),
    }
    # Initial setup
    _spec = specabbr.upper()
    res = unit_to.split(':')
    dfact = 1.0
    sfact = 1.0
    if len(res) == 1:
        _dest_unit = ''
    else:
        dunit = res[1].split()
        if len(dunit) == 1:
            _dest_unit = res[1].strip().replace('^', '')
            dfact = 1.0
        else:
            _dest_unit = dunit[1].strip().replace('^', '')
            dfact = eval(convert_expr(dunit[0], None))
    _dest_type = res[0].upper()
    res = unit_from.split(':')
    if len(res) == 1:
        _src_unit = ''
    else:
        sunit = res[1].split()
        if len(sunit) == 1:
            _src_unit = res[1].strip().replace('^', '')
            sfact = 1.0
        else:
            _src_unit = sunit[1].strip().replace('^', '')
            sfact = eval(convert_expr(sunit[0], None))
    _src_type = res[0].upper()
    yfactor = 1.0
    # Test
    if _spec == 'IR':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'DS':
                    if _src_unit in UNIT_DS['esu2.cm2']:
                        yfactor = 8*pi**3*PHYSCNST.avogadro * 1.0e-7 / \
                            (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                elif _src_type == 'II':
                    if _src_unit in UNIT_II['km/M']:
                        yfactor = 100/log(10)

                        xfunc = None
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'VCD':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'RS':
                    if _src_unit in UNIT_RS['esu2.cm2']:
                        yfactor = 32*pi**3*PHYSCNST.avogadro * 1.0e-7 / \
                            (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'RS':
        incfrq = subopts.get('incfreq', None)
        if incfrq is None:
            incfrq = 1.0e7/532.0
        elif not isinstance(incfrq, (int, float)):
            raise ArgumentError(msgVal.format(subkey='incfreq',
                                              specabbr=specabbr))
        if _dest_type == 'I':
            if _dest_unit in UNIT_RDSG['cm3/mol/sr']:
                if _src_type == 'RA':
                    # (AA->cm)^6 * Na * (2 pi * Wscat)^4
                    yfactor = 1.0e-48*PHYSCNST.avogadro*(2*pi)**4 / 45

                    def xfunc(x):
                        return (incfrq-x)**4
                    if _src_unit in UNIT_RA['amu.Ang4']:
                        def ycorr(x):  # (Q->q)^2
                            return phys_fact('mwq2q')**2*PHYSFACT.bohr2ang**2 \
                                / (2.0*x)
                    elif _src_unit in UNIT_RA['Ang6']:
                        ycorr = None
                    elif _src_unit in UNIT_RA['bohr6']:
                        # (bohr->AA)^6 * yfactor
                        yfactor *= PHYSFACT.bohr2ang**6
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            elif _dest_unit in UNIT_RDSG['m2.cm/mol/sr']:
                if _src_type == 'RA':
                    # 10^-4 * (AA->cm)^6 * Na * (2 pi * Wscat)^4
                    yfactor = 1.0e-52*PHYSCNST.avogadro*(2*pi)**4 / 45

                    def xfunc(x):
                        return (incfrq-x)**4
                    if _src_unit in UNIT_RA['amu.Ang4']:
                        def ycorr(x):  # (Q->q)^2
                            return phys_fact('mwq2q')**2*PHYSFACT.bohr2ang**2 \
                                / (2.0*x)
                    elif _src_unit in UNIT_RA['Ang6']:
                        ycorr = None
                    elif _src_unit in UNIT_RA['bohr6']:
                        # (bohr->AA)^6 * yfactor
                        yfactor *= PHYSFACT.bohr2ang**6
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'ROA':
        incfrq = subopts.get('incfreq', None)
        if incfrq is None:
            incfrq = 1.0e7/532.0
        elif not isinstance(incfrq, (int, float)):
            raise ArgumentError(msgVal.format(subkey='incfreq',
                                              specabbr=specabbr))
        if _dest_type == 'I':
            if _dest_unit in UNIT_RDSG['cm3/mol/sr']:
                if _src_type == 'ROA':
                    # (AA->cm)^6 * Na * (2 pi * Wscat)^4
                    yfactor = 1.0e-48*PHYSCNST.avogadro*(2*pi)**4 / 45

                    def xfunc(x):
                        return (incfrq-x)**4
                    if _src_unit in UNIT_ROA['amu.Ang4']:
                        def ycorr(x):  # (Q->q)^2
                            return phys_fact('mwq2q')**2*PHYSFACT.bohr2ang**2 \
                                / (2.0*x)
                    elif _src_unit in UNIT_ROA['Ang6']:
                        ycorr = None
                    elif _src_unit in UNIT_ROA['bohr6']:
                        # (bohr->AA)^6 * yfactor
                        yfactor *= PHYSFACT.bohr2ang**6
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            elif _dest_unit in UNIT_RDSG['m2.cm/mol/sr']:
                if _src_type == 'ROA':
                    # 10^-4 * (AA->cm)^6 * Na * (2 pi * Wscat)^4
                    yfactor = 1.0e-52*PHYSCNST.avogadro*(2*pi)**4 / 45

                    def xfunc(x):
                        return (incfrq-x)**4
                    if _src_unit in UNIT_ROA['amu.Ang4']:
                        def ycorr(x):  # (Q->q)^2
                            return phys_fact('mwq2q')**2*PHYSFACT.bohr2ang**2 \
                                / (2.0*x)
                    elif _src_unit in UNIT_ROA['Ang6']:
                        ycorr = None
                    elif _src_unit in UNIT_ROA['bohr6']:
                        # (bohr->AA)^6 * yfactor
                        yfactor *= PHYSFACT.bohr2ang**6
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'OPA':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'DS':
                    if _src_unit in UNIT_DS['esu2.cm2']:
                        yfactor = 8*pi**3*PHYSCNST.avogadro * 1.0e-7 \
                            / (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                        ycorr = None
                    elif _src_unit in UNIT_DS['au']:
                        yfactor = 8*pi**3*PHYSCNST.avogadro * 1.0e-7 \
                            * (phys_fact('au2esu')*PHYSFACT.bohr2ang
                               * 1.0e-8)**2 \
                            / (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                elif _src_type == 'II':
                    if _src_unit in UNIT_II['/M/cm2']:
                        yfactor = 1.0
                        xfunc = None
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    elif _spec == 'ECD':
        if _dest_type == 'I':
            if _dest_unit in UNIT_EPS['/M/cm']:
                if _src_type == 'RS':
                    if _src_unit in UNIT_RS['esu2.cm2']:
                        yfactor = 32*pi**3*PHYSCNST.avogadro * 1.0e-7 / \
                            (3000.*PHYSCNST.planck*PHYSCNST.slight*log(10))

                        def xfunc(x): return x
                        ycorr = None
                    else:
                        raise NotImplementedError(msgNYI)
                else:
                    raise NotImplementedError(msgNYI)
            else:
                raise IndexError(msgNA)
        else:
            raise IndexError(msgNA)
    else:
        raise NotImplementedError(msgNYI)

    return yfactor*sfact/dfact, xfunc, ycorr
